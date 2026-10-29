/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.function.Consumer;

/**
 * The set of fragment-pass grouping keys (library, read1 reference index, read1 coordinate, orientation) of paired
 * read ends. MarkDuplicates' fragment pass marks only unpaired ends; a paired end there only tells its chunk that a
 * pair shares the chunk's key, which this set answers without keeping the paired end as an object.
 *
 * Each key is packed into a non-negative long whose numeric order is the lexicographic order of the four fields.
 * A key whose library id or reference index does not fit is not added, and its end is kept by the caller.
 *
 * Keys are added into fixed-size chunks; a full chunk is sorted, de-duplicated and placed into an ordered list of
 * sorted runs on a helper thread while adding goes on, so that doneAdding has at most one partial chunk to sort.
 *
 * A run's array is large enough that G1 allocates it as a humongous object, which it never moves. Before an allocation
 * that needs a large contiguous block, {@link #splitRuns()} moves the keys into small arrays, and the overflow rule
 * moves them the same way before it hands them over.
 */
public final class PairedReadEndKeys {
    private static final int COORDINATE_AND_ORIENTATION_BITS = 33;
    private static final int INITIAL_CAPACITY = 1 << 16;
    private static final int CHUNK_SIZE = 1 << 20;
    // 256 KiB of keys: under half of G1's smallest region (1 MiB), so never a humongous object.
    static final int SPLIT_KEYS = 1 << 15;

    private final int referenceBits;
    private final int libraryBits;
    private final int maxKeysInRam;
    private final int chunkSize;
    private final Consumer<ReadEndsForMarkDuplicates> overflow;

    private long[] chunk;
    private int chunkLength = 0;
    /** The keys held as the overflow rule counts them: the last full sort's distinct keys, plus those added since. */
    private int held = 0;
    private long lastAdded;
    private boolean doneAdding = false;
    private long handedOver = 0;

    private ExecutorService sorter;
    private final List<Future<?>> sorting = new ArrayList<>();
    private final ConcurrentLinkedQueue<long[]> freeChunks = new ConcurrentLinkedQueue<>();

    /** Sorted, distinct and each above the one before; written by whichever thread places a chunk. */
    private final List<Run> runs = new ArrayList<>();
    /** Sorted, distinct runs that overlapped too many keys of the runs to merge in; sorted with them in gather. */
    private final List<Run> unplaced = new ArrayList<>();

    private long[][] runKeys;
    private int[] runFrom;
    private int[] runTo;
    private long[] runLast;
    private int cursorRun = 0;
    private int cursor = 0;
    private long lastQuery = Long.MIN_VALUE;

    private static final class Run {
        final long[] keys;
        final int from;
        int to;

        Run(final long[] keys, final int from, final int to) {
            this.keys = keys;
            this.from = from;
            this.to = to;
        }

        long last() {
            return keys[to - 1];
        }
    }

    /**
     * @param numberOfReferences the number of sequences in the input's header
     * @param maxKeysInRam       the number of keys held; when it is reached and half of them are distinct, every
     *                           key held is handed to overflow as a paired read end and dropped from this set
     * @param overflow           receives those read ends; one carries only the four key fields, with
     *                           read2ReferenceIndex set so that it is paired
     */
    public PairedReadEndKeys(final int numberOfReferences, final int maxKeysInRam,
                             final Consumer<ReadEndsForMarkDuplicates> overflow) {
        this(numberOfReferences, maxKeysInRam, overflow, CHUNK_SIZE);
    }

    PairedReadEndKeys(final int numberOfReferences, final int maxKeysInRam,
                      final Consumer<ReadEndsForMarkDuplicates> overflow, final int chunkSize) {
        this.referenceBits = numberOfReferences <= 1 ? 0 : 32 - Integer.numberOfLeadingZeros(numberOfReferences - 1);
        this.libraryBits = 63 - COORDINATE_AND_ORIENTATION_BITS - referenceBits;
        this.maxKeysInRam = Math.max(1, maxKeysInRam);
        this.chunkSize = Math.max(1, chunkSize);
        this.overflow = overflow;
        this.chunk = new long[initialCapacity()];
    }

    private int initialCapacity() {
        return Math.min(INITIAL_CAPACITY, Math.min(this.maxKeysInRam, this.chunkSize));
    }

    /** The packed key of the four fields, or -1 when they cannot be packed. */
    private long pack(final short libraryId, final int referenceIndex, final int coordinate, final byte orientation) {
        if (libraryBits < 0 || libraryId < 0 || libraryId >= (1L << libraryBits) ||
                referenceIndex < 0 || (referenceIndex >>> referenceBits) != 0 ||
                (orientation != ReadEnds.F && orientation != ReadEnds.R)) {
            return -1;
        }
        final long reference = ((long) libraryId << referenceBits) | referenceIndex;
        return (reference << COORDINATE_AND_ORIENTATION_BITS) |
                (((long) coordinate - Integer.MIN_VALUE) << 1) | orientation;
    }

    private long keyOf(final ReadEndsForMarkDuplicates end) {
        return pack(end.libraryId, end.read1ReferenceIndex, end.read1Coordinate, end.orientation);
    }

    /**
     * Adds the key of a paired read end.
     *
     * @return false, adding nothing, when the key cannot be packed
     */
    public boolean add(final ReadEndsForMarkDuplicates end) {
        if (doneAdding) {
            throw new IllegalStateException("Cannot add after doneAdding()");
        }
        final long key = keyOf(end);
        if (key < 0) {
            return false;
        }
        if (held > 0 && lastAdded == key) {
            return true;
        }
        if (held == maxKeysInRam) {
            gather();
            int distinct = 0;
            for (final Run run : runs) {
                distinct += run.to - run.from;
            }
            if (distinct > maxKeysInRam / 2) {
                // Each array is released once its keys are handed over, and those not yet handed over can be moved.
                final List<long[]> pieces = split(runs, SPLIT_KEYS, null);
                for (int p = 0; p < pieces.size(); p++) {
                    final long[] keys = pieces.set(p, null);
                    for (final long heldKey : keys) {
                        overflow.accept(unpack(heldKey));
                    }
                    handedOver += keys.length;
                }
                distinct = 0;
            }
            held = distinct;
            chunk = new long[initialCapacity()];
        } else if (chunkLength == chunk.length) {
            if (chunk.length < chunkSize) {
                chunk = Arrays.copyOf(chunk, (int) Math.min(2L * chunk.length, Math.min(chunkSize, maxKeysInRam)));
            } else {
                submit(chunk, chunkLength);
                final long[] free = freeChunks.poll();
                chunk = free != null ? free : new long[chunkSize];
                chunkLength = 0;
            }
        }
        chunk[chunkLength++] = key;
        lastAdded = key;
        held++;
        return true;
    }

    private ReadEndsForMarkDuplicates unpack(final long key) {
        final ReadEndsForMarkDuplicates end = new ReadEndsForMarkDuplicates();
        end.libraryId = (short) (key >>> (COORDINATE_AND_ORIENTATION_BITS + referenceBits));
        end.read1ReferenceIndex = (int) ((key >>> COORDINATE_AND_ORIENTATION_BITS) & ((1L << referenceBits) - 1));
        end.read1Coordinate = (int) (((key >>> 1) & 0xFFFFFFFFL) + Integer.MIN_VALUE);
        end.orientation = (byte) (key & 1);
        end.read2ReferenceIndex = end.read1ReferenceIndex;
        return end;
    }

    private void submit(final long[] keys, final int length) {
        if (sorter == null) {
            sorter = new ThreadPoolExecutor(0, 1, 1, TimeUnit.SECONDS, new LinkedBlockingQueue<>(), r -> {
                final Thread thread = new Thread(r, "PairedReadEndKeys sorter");
                thread.setDaemon(true);
                return thread;
            });
        }
        sorting.add(sorter.submit(() -> {
            place(keys, length, false);
            freeChunks.offer(keys);
        }));
    }

    /**
     * Sorts keys[0, length) and places its distinct keys into runs: after them when every key is above the runs'
     * last, otherwise merged with the runs' keys at or above its first key when those are few, otherwise into
     * unplaced. On the adding thread keys are sorted in parallel and become the runs' storage; on the sorter they
     * are copied, so that keys can be reused.
     */
    private void place(final long[] keys, final int length, final boolean onAddingThread) {
        if (onAddingThread) {
            Arrays.parallelSort(keys, 0, length);
        } else {
            Arrays.sort(keys, 0, length);
        }
        final int distinct = removeDuplicates(keys, 0, length);
        final long first = keys[0];

        // The runs' keys at or above first start at run startRun, index startIndex.
        int startRun = runs.size();
        int startIndex = 0;
        long overlapping = 0;
        while (startRun > 0 && runs.get(startRun - 1).last() >= first) {
            final Run run = runs.get(startRun - 1);
            final int index = lowerBound(run.keys, run.from, run.to, first);
            overlapping += run.to - index;
            if (overlapping > chunkSize) {
                unplaced.add(onAddingThread ? new Run(keys, 0, distinct) :
                        new Run(Arrays.copyOf(keys, distinct), 0, distinct));
                return;
            }
            startRun--;
            startIndex = index;
            if (index > run.from) {
                break;
            }
        }

        if (overlapping == 0) {
            runs.add(onAddingThread ? new Run(keys, 0, distinct) :
                    new Run(Arrays.copyOf(keys, distinct), 0, distinct));
            return;
        }

        // Merge the overlapping keys with keys[0, above), every key not above the runs' last.
        final long last = runs.get(runs.size() - 1).last();
        final int above = upperBound(keys, 0, distinct, last);
        final long[] old = new long[(int) overlapping];
        int n = 0;
        for (int r = startRun; r < runs.size(); r++) {
            final Run run = runs.get(r);
            final int from = r == startRun ? startIndex : run.from;
            System.arraycopy(run.keys, from, old, n, run.to - from);
            n += run.to - from;
        }
        final long[] merged = new long[n + above];
        int i = 0, j = 0, m = 0;
        while (i < n && j < above) {
            final long a = old[i];
            final long b = keys[j];
            if (a < b) {
                merged[m++] = a;
                i++;
            } else if (b < a) {
                merged[m++] = b;
                j++;
            } else {
                merged[m++] = a;
                i++;
                j++;
            }
        }
        while (i < n) {
            merged[m++] = old[i++];
        }
        while (j < above) {
            merged[m++] = keys[j++];
        }

        while (runs.size() > startRun + 1) {
            runs.remove(runs.size() - 1);
        }
        final Run startingRun = runs.get(startRun);
        if (startIndex == startingRun.from) {
            runs.remove(startRun);
        } else {
            startingRun.to = startIndex;
        }
        runs.add(new Run(merged, 0, m));
        if (above < distinct) {
            runs.add(onAddingThread ? new Run(keys, above, distinct) :
                    new Run(Arrays.copyOfRange(keys, above, distinct), 0, distinct - above));
        }
    }

    /** The first index in [from, to) whose key is at least key, or to. */
    private static int lowerBound(final long[] keys, int from, int to, final long key) {
        while (from < to) {
            final int middle = (from + to) >>> 1;
            if (keys[middle] < key) {
                from = middle + 1;
            } else {
                to = middle;
            }
        }
        return from;
    }

    /** The first index in [from, to) whose key is above key, or to. */
    private static int upperBound(final long[] keys, int from, int to, final long key) {
        while (from < to) {
            final int middle = (from + to) >>> 1;
            if (keys[middle] <= key) {
                from = middle + 1;
            } else {
                to = middle;
            }
        }
        return from;
    }

    /** Removes repeats from the sorted keys[from, to); returns the new end. */
    private static int removeDuplicates(final long[] keys, final int from, final int to) {
        int distinct = from;
        for (int i = from; i < to; i++) {
            if (distinct == from || keys[i] != keys[distinct - 1]) {
                keys[distinct++] = keys[i];
            }
        }
        return distinct;
    }

    /** Waits for the helper, then places the current chunk: runs hold every key held, and chunk is free. */
    private void gather() {
        boolean interrupted = false;
        for (final Future<?> future : sorting) {
            while (true) {
                try {
                    future.get();
                    break;
                } catch (final InterruptedException e) {
                    interrupted = true;
                } catch (final ExecutionException e) {
                    final Throwable cause = e.getCause();
                    if (cause instanceof Error) {
                        throw (Error) cause;
                    }
                    if (cause instanceof RuntimeException) {
                        throw (RuntimeException) cause;
                    }
                    throw new RuntimeException(cause);
                }
            }
        }
        if (interrupted) {
            Thread.currentThread().interrupt();
        }
        sorting.clear();
        if (chunkLength > 0) {
            place(chunk, chunkLength, true);
            chunkLength = 0;
        }
        chunk = null;
        if (!unplaced.isEmpty()) {
            long total = 0;
            for (final Run run : runs) {
                total += run.to - run.from;
            }
            for (final Run run : unplaced) {
                total += run.to - run.from;
            }
            final long[] all = new long[(int) total];
            int n = 0;
            for (final List<Run> list : Arrays.asList(runs, unplaced)) {
                for (final Run run : list) {
                    System.arraycopy(run.keys, run.from, all, n, run.to - run.from);
                    n += run.to - run.from;
                }
            }
            runs.clear();
            unplaced.clear();
            Arrays.parallelSort(all, 0, n);
            runs.add(new Run(all, 0, removeDuplicates(all, 0, n)));
        }
    }

    /**
     * The keys of the runs, in order, in arrays of at most maxKeys keys each, every array exactly full; a run whose array
     * is already that small and holds only its keys keeps it. Each run is dropped from runs once copied, so that its
     * array can be collected before the next is copied, and runs is left empty.
     *
     * @param firstPiece when not null, set to the index of each run's first array
     */
    private static List<long[]> split(final List<Run> runs, final int maxKeys, final int[] firstPiece) {
        final List<long[]> pieces = new ArrayList<>();
        for (int r = 0; r < runs.size(); r++) {
            final Run run = runs.set(r, null);
            if (firstPiece != null) {
                firstPiece[r] = pieces.size();
            }
            if (run.from == 0 && run.to == run.keys.length && run.keys.length <= maxKeys) {
                pieces.add(run.keys);
                continue;
            }
            for (int from = run.from; from < run.to; from += maxKeys) {
                pieces.add(Arrays.copyOfRange(run.keys, from, (int) Math.min(run.to, (long) from + maxKeys)));
            }
        }
        runs.clear();
        return pieces;
    }

    /**
     * Moves the keys into arrays of at most 32,768 keys each, small enough that G1 can move them, releasing each run's
     * array once its keys are copied. Call after doneAdding, before an allocation that needs a large contiguous block of
     * heap; contains answers as before.
     */
    public void splitRuns() {
        splitRuns(SPLIT_KEYS);
    }

    void splitRuns(final int maxKeys) {
        if (!doneAdding) {
            throw new IllegalStateException("Cannot split before doneAdding()");
        }
        if (maxKeys < 1) {
            throw new IllegalArgumentException("maxKeys must be positive: " + maxKeys);
        }
        final int n = runKeys.length;
        final List<Run> queried = new ArrayList<>(n);
        for (int r = 0; r < n; r++) {
            queried.add(new Run(runKeys[r], runFrom[r], runTo[r]));
            runKeys[r] = null;
        }
        runs.clear();
        final int[] firstPiece = new int[n];
        final List<long[]> pieces = split(queried, maxKeys, firstPiece);

        // The cursor moves to the same key in the array that now holds it.
        final int m = pieces.size();
        if (cursorRun < n) {
            final int offset = cursor - runFrom[cursorRun];
            cursorRun = firstPiece[cursorRun] + offset / maxKeys;
            cursor = offset % maxKeys;
        } else {
            cursorRun = m;
            cursor = 0;
        }
        runKeys = pieces.toArray(new long[m][]);
        runFrom = new int[m];
        runTo = new int[m];
        runLast = new long[m];
        for (int p = 0; p < m; p++) {
            runTo[p] = runKeys[p].length;
            runLast[p] = runKeys[p][runTo[p] - 1];
        }
    }

    /** The keys handed to overflow so far. */
    public long keysHandedOver() {
        return handedOver;
    }

    /** Sorts the keys; call once, after the last add and before contains. */
    public void doneAdding() {
        gather();
        if (sorter != null) {
            sorter.shutdown();
            sorter = null;
        }
        freeChunks.clear();
        final int n = runs.size();
        runKeys = new long[n][];
        runFrom = new int[n];
        runTo = new int[n];
        runLast = new long[n];
        for (int r = 0; r < n; r++) {
            final Run run = runs.get(r);
            runKeys[r] = run.keys;
            runFrom[r] = run.from;
            runTo[r] = run.to;
            runLast[r] = run.last();
        }
        cursorRun = 0;
        cursor = n > 0 ? runFrom[0] : 0;
        doneAdding = true;
    }

    /**
     * Whether a paired read end added here, and not handed to overflow, has the key of the given end. Exact in any
     * call order; fastest when calls come in nondecreasing key order, as in a pass over ReadEndsMDComparator's order.
     */
    public boolean contains(final ReadEndsForMarkDuplicates end) {
        if (!doneAdding) {
            throw new IllegalStateException("Cannot query before doneAdding()");
        }
        final long key = keyOf(end);
        if (key < 0) {
            return false;
        }
        final int n = runKeys.length;
        if (key < lastQuery) {
            final int r = lowerBound(runLast, 0, n, key);
            return r < n && Arrays.binarySearch(runKeys[r], runFrom[r], runTo[r], key) >= 0;
        }
        lastQuery = key;
        while (cursorRun < n && runLast[cursorRun] < key) {
            cursorRun++;
            cursor = cursorRun < n ? runFrom[cursorRun] : 0;
        }
        if (cursorRun == n) {
            return false;
        }
        final long[] keys = runKeys[cursorRun];
        final int size = runTo[cursorRun];
        if (keys[cursor] < key) {
            // Gallop: keys[low] < key, then search between low and the first step that reaches key.
            int low = cursor;
            long step = 1;
            long high = low + step;
            while (high < size && keys[(int) high] < key) {
                low = (int) high;
                step <<= 1;
                high = low + step;
            }
            final int found = Arrays.binarySearch(keys, low + 1, (int) Math.min(high, size), key);
            cursor = found >= 0 ? found : -found - 1;
        }
        return keys[cursor] == key;
    }
}
