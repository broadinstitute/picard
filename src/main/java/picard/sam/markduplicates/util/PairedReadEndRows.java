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

import java.util.Arrays;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Consumer;
import java.util.stream.IntStream;

/**
 * Paired read ends held as rows of primitive fields instead of as objects, for MarkDuplicates' pair pass, which
 * would otherwise keep every completed pair as an object from the first pass to the end of the pair pass.
 *
 * {@link #sort()} orders the rows as MarkDuplicates.ReadEndsMDComparator without barcodes orders the ends, and ends
 * it compares as equal in the order they were added, as a stable sort does: library id, read1 reference index and
 * coordinate, orientation, read2 reference index and coordinate, tile, x, y, read1 and read2 index in the file. That
 * comparator subtracts ints, which gives this order only while each int it subtracts lies in [-2^30, 2^30).
 * {@link #add} refuses an end outside that range, one past the capacity, and one it could not give back unchanged;
 * the caller then keeps the ends as objects, starting with the rows already added ({@link #drainTo}).
 */
public final class PairedReadEndRows {
    // 7 longs used. Rows start at index 6 of a chunk: past a 16-byte array header, at a 64-byte boundary when the
    // array is, so that a row is one cache line.
    private static final int LONGS_PER_ROW = 8;
    private static final int USED_LONGS = 7;
    private static final int ROW_BASE = 6;
    // Chunks of up to 2^20 rows, 64 MiB. G1 puts an array this large straight into the old generation, where young
    // collections do not copy it, and may start a concurrent marking cycle on each one, so chunks are few; but no
    // more than an eighth of the capacity, so that a small heap holds no more than its capacity's rows.
    private static final int MAX_CHUNK_SHIFT = 20;
    private static final int MIN_CHUNK_SHIFT = 10;
    // Flips the sign bit of the low int of a long holding two ints, so that signed long order is their order.
    private static final long LOW_SIGN = 0x80000000L;
    // After sort(), set in the order entry of each position that starts a duplicate set.
    private static final long SET_START = Long.MIN_VALUE;
    private static final int INSERTION_SORT_MAX = 16;
    private static final int ROWS_PER_SORT_TASK = 1 << 16;

    private final int capacity;
    private final int chunkShift;
    private final int chunkRows;
    private final int chunkMask;
    private long[][] chunks = new long[1][];
    private int size = 0;

    // The ranges of the four leading sort fields, for packing them into sort keys.
    private int minLibrary = Integer.MAX_VALUE, maxLibrary = Integer.MIN_VALUE;
    private int minReference = Integer.MAX_VALUE, maxReference = Integer.MIN_VALUE;
    private int minCoordinate = Integer.MAX_VALUE, maxCoordinate = Integer.MIN_VALUE;
    private int minOrientation = Integer.MAX_VALUE, maxOrientation = Integer.MIN_VALUE;

    // During sort(), a prefix of each position's leading sort fields above its row, in the low rowBits bits;
    // afterwards the row and SET_START.
    private long[] order;
    private int rowBits;
    private long rowMask;
    // Whether the prefix holds the four leading sort fields whole, so that rows with equal prefixes have them equal.
    private boolean wholePrefix;
    // After sort(), the ends of duplicate sets of more than one end in sorted order, copied out of the rows so that
    // getAt reads them in sequence: their positions, and their 7 fields each.
    private int[] setPositions = new int[0];
    private long[] setEnds = new long[0];
    private int setCursor = 0;

    /** @param capacity the number of rows held; an end added past it is refused */
    public PairedReadEndRows(final int capacity) {
        this.capacity = Math.max(0, capacity);
        final int eighthShift = 31 - Integer.numberOfLeadingZeros(Math.max(1, this.capacity / 8));
        this.chunkShift = Math.max(MIN_CHUNK_SHIFT, Math.min(MAX_CHUNK_SHIFT, eighthShift));
        this.chunkRows = 1 << chunkShift;
        this.chunkMask = chunkRows - 1;
    }

    public int size() {
        return size;
    }

    private static boolean outsideComparatorRange(final int value) {
        return value + (1 << 30) < 0;
    }

    private static long twoInts(final int high, final int low) {
        return ((long) high << 32) | (low & 0xFFFFFFFFL);
    }

    private int at(final int row) {
        return ROW_BASE + (row & chunkMask) * LONGS_PER_ROW;
    }

    private long[] chunk(final int row) {
        return chunks[row >>> chunkShift];
    }

    /**
     * Adds a paired read end as a row.
     *
     * @return false, adding nothing, when the end is past the capacity, is not a plain ReadEndsForMarkDuplicates with
     * duplicateSetSize -1 and isOpticalDuplicate false, or holds an int the comparator subtracts outside [-2^30, 2^30)
     */
    public boolean add(final ReadEndsForMarkDuplicates end) {
        if (order != null) {
            throw new IllegalStateException("Cannot add after sort()");
        }
        if (size >= capacity || end.getClass() != ReadEndsForMarkDuplicates.class ||
                end.duplicateSetSize != -1 || end.isOpticalDuplicate ||
                outsideComparatorRange(end.read1ReferenceIndex) || outsideComparatorRange(end.read1Coordinate) ||
                outsideComparatorRange(end.read2ReferenceIndex) || outsideComparatorRange(end.read2Coordinate) ||
                outsideComparatorRange(end.x) || outsideComparatorRange(end.y)) {
            return false;
        }
        final int row = size;
        final long[] chunk = chunkFor(row);
        final int at = at(row);
        chunk[at] = end.read1IndexInFile;
        chunk[at + 1] = end.read2IndexInFile;
        chunk[at + 2] = twoInts(end.read1ReferenceIndex, end.read1Coordinate);
        chunk[at + 3] = twoInts(end.read2ReferenceIndex, end.read2Coordinate);
        chunk[at + 4] = twoInts(end.x, end.y);
        chunk[at + 5] = ((long) end.libraryId << 48) | ((end.tile & 0xFFFFL) << 32) |
                ((end.readGroup & 0xFFFFL) << 16) | (end.score & 0xFFFFL);
        chunk[at + 6] = ((end.orientation & 0xFFL) << 8) | (end.orientationForOpticalDuplicates & 0xFFL);
        minLibrary = Math.min(minLibrary, end.libraryId);
        maxLibrary = Math.max(maxLibrary, end.libraryId);
        minReference = Math.min(minReference, end.read1ReferenceIndex);
        maxReference = Math.max(maxReference, end.read1ReferenceIndex);
        minCoordinate = Math.min(minCoordinate, end.read1Coordinate);
        maxCoordinate = Math.max(maxCoordinate, end.read1Coordinate);
        minOrientation = Math.min(minOrientation, end.orientation);
        maxOrientation = Math.max(maxOrientation, end.orientation);
        size = row + 1;
        return true;
    }

    private long[] chunkFor(final int row) {
        final int c = row >>> chunkShift;
        if (c == chunks.length) {
            chunks = Arrays.copyOf(chunks, 2 * c);
        }
        long[] chunk = chunks[c];
        if (chunk == null) {
            chunk = new long[ROW_BASE + (c == 0 ? 1 << MIN_CHUNK_SHIFT : chunkRows) * LONGS_PER_ROW];
            chunks[c] = chunk;
        } else if (chunk.length < at(row) + LONGS_PER_ROW) {
            chunk = Arrays.copyOf(chunk, ROW_BASE + Math.min(2 * (chunk.length - ROW_BASE), chunkRows * LONGS_PER_ROW));
            chunks[c] = chunk;
        }
        return chunk;
    }

    /** The end added as the given row (counting from 0), as a new object equal in every field to the one added. */
    public ReadEndsForMarkDuplicates get(final int row) {
        return toEnd(chunk(row), at(row));
    }

    private static ReadEndsForMarkDuplicates toEnd(final long[] chunk, final int at) {
        final ReadEndsForMarkDuplicates end = new ReadEndsForMarkDuplicates();
        end.read1IndexInFile = chunk[at];
        end.read2IndexInFile = chunk[at + 1];
        end.read1ReferenceIndex = (int) (chunk[at + 2] >> 32);
        end.read1Coordinate = (int) chunk[at + 2];
        end.read2ReferenceIndex = (int) (chunk[at + 3] >> 32);
        end.read2Coordinate = (int) chunk[at + 3];
        end.x = (int) (chunk[at + 4] >> 32);
        end.y = (int) chunk[at + 4];
        final long fields = chunk[at + 5];
        end.libraryId = (short) (fields >>> 48);
        end.tile = (short) (fields >>> 32);
        end.readGroup = (short) (fields >>> 16);
        end.score = (short) fields;
        end.orientation = (byte) (chunk[at + 6] >>> 8);
        end.orientationForOpticalDuplicates = (byte) chunk[at + 6];
        return end;
    }

    /** Hands every row to the consumer as an object, in the order added, and empties this; call before sort(). */
    public void drainTo(final Consumer<ReadEndsForMarkDuplicates> consumer) {
        if (order != null) {
            throw new IllegalStateException("Cannot drain after sort()");
        }
        final int n = size;
        size = 0;
        for (int row = 0; row < n; row++) {
            consumer.accept(get(row));
            if ((row & chunkMask) == chunkMask || row == n - 1) {
                chunks[row >>> chunkShift] = null;
            }
        }
        chunks = new long[1][];
    }

    private static int bitsFor(final long range) {
        return 64 - Long.numberOfLeadingZeros(range);
    }

    /** Sorts the rows and marks where duplicate sets start; call once, after the last add and before rowAt. */
    public void sort() {
        if (order != null) {
            throw new IllegalStateException("sort() was called already");
        }
        final int n = size;
        rowBits = n <= 1 ? 1 : 32 - Integer.numberOfLeadingZeros(n - 1);
        rowMask = (1L << rowBits) - 1;
        order = new long[n];
        if (n == 0) {
            return;
        }

        // The leading fields' offsets from their minimums, packed in order; bits dropped from the right end to fit
        // beside the row keep the keys in order, and rows whose prefixes tie are ordered by the full comparison.
        final int[] widths = {bitsFor((long) maxLibrary - minLibrary), bitsFor((long) maxReference - minReference),
                bitsFor((long) maxCoordinate - minCoordinate), bitsFor((long) maxOrientation - minOrientation)};
        final int[] dropped = new int[4];
        int excess = widths[0] + widths[1] + widths[2] + widths[3] - (63 - rowBits);
        wholePrefix = excess <= 0;
        for (int f = 3; f >= 0 && excess > 0; f--) {
            dropped[f] = Math.min(excess, widths[f]);
            widths[f] -= dropped[f];
            excess -= dropped[f];
        }
        final int chunkCount = (n - 1 >>> chunkShift) + 1;
        IntStream.range(0, chunkCount).parallel().forEach(c -> {
            final long[] chunk = chunks[c];
            final int first = c << chunkShift;
            final int last = Math.min(n, first + chunkRows);
            for (int row = first, at = ROW_BASE; row < last; row++, at += LONGS_PER_ROW) {
                final long read1 = chunk[at + 2];
                long prefix = (((short) (chunk[at + 5] >>> 48)) - (long) minLibrary) >>> dropped[0];
                prefix = (prefix << widths[1]) | ((((int) (read1 >> 32)) - (long) minReference) >>> dropped[1]);
                prefix = (prefix << widths[2]) | ((((int) read1) - (long) minCoordinate) >>> dropped[2]);
                prefix = (prefix << widths[3]) |
                        ((((byte) (chunk[at + 6] >>> 8)) - (long) minOrientation) >>> dropped[3]);
                order[row] = (prefix << rowBits) | row;
            }
        });
        Arrays.parallelSort(order);

        // Tasks start where a prefix does, found before any task rewrites its part of order.
        final int tasks = n < 2 * ROWS_PER_SORT_TASK ? 1 :
                Math.min(4 * ForkJoinPool.getCommonPoolParallelism(), n / ROWS_PER_SORT_TASK);
        final int[] starts = new int[tasks + 1];
        for (int t = 0; t <= tasks; t++) {
            starts[t] = prefixStart((int) ((long) t * n / tasks));
        }
        final int[] inSets = new int[tasks + 1];
        IntStream.range(0, tasks).parallel()
                .forEach(t -> inSets[t + 1] = sortTiesAndMarkSets(starts[t], starts[t + 1]));
        for (int t = 0; t < tasks; t++) {
            inSets[t + 1] += inSets[t];
        }
        // Without the copies getAt reads the rows: when they would take more than a sixteenth of the heap.
        if ((long) inSets[tasks] * (USED_LONGS * Long.BYTES + Integer.BYTES) > Runtime.getRuntime().maxMemory() / 16 ||
                (long) inSets[tasks] * USED_LONGS > Integer.MAX_VALUE - 8) {
            return;
        }
        setPositions = new int[inSets[tasks]];
        setEnds = new long[inSets[tasks] * USED_LONGS];
        IntStream.range(0, tasks).parallel().forEach(t -> copySetEnds(starts[t], starts[t + 1], inSets[t]));
    }

    /** Copies the ends of duplicate sets of more than one end in [from, to) out of the rows, from index on. */
    private void copySetEnds(final int from, final int to, int index) {
        for (int position = from; position < to; position++) {
            if (order[position] < 0 && (position + 1 == to || order[position + 1] < 0)) {
                continue;
            }
            final int row = rowAt(position);
            setPositions[index] = position;
            System.arraycopy(chunk(row), at(row), setEnds, index * USED_LONGS, USED_LONGS);
            index++;
        }
    }

    /** The first position at or after the given one that starts a run of equal prefixes. */
    private int prefixStart(int position) {
        while (position > 0 && position < order.length &&
                (order[position] >>> rowBits) == (order[position - 1] >>> rowBits)) {
            position++;
        }
        return position;
    }

    /**
     * Orders each run of equal prefixes in [from, to) by the full comparison and rewrites each entry as its row, with
     * SET_START where a duplicate set starts: at a new prefix, which differs in a leading field, or where a row is
     * not in its predecessor's set.
     *
     * @return the number of positions in duplicate sets of more than one end
     */
    private int sortTiesAndMarkSets(final int from, final int to) {
        int inSets = 0;
        long[] rows = new long[64], keys = new long[64], rowsScratch = new long[64], keysScratch = new long[64];
        for (int start = from, end; start < to; start = end) {
            final long prefix = order[start] >>> rowBits;
            end = start + 1;
            while (end < to && (order[end] >>> rowBits) == prefix) {
                end++;
            }
            final int m = end - start;
            if (m == 1) {
                order[start] = (order[start] & rowMask) | SET_START;
                continue;
            }
            if (rows.length < m) {
                final int length = Math.max(m, 2 * rows.length);
                rows = new long[length];
                keys = new long[length];
                rowsScratch = new long[length];
                keysScratch = new long[length];
            }
            // Read2 first, from each row once: with whole prefixes it decides most comparisons in the run.
            for (int i = 0; i < m; i++) {
                final int row = (int) (order[start + i] & rowMask);
                rows[i] = row;
                keys[i] = chunk(row)[at(row) + 3] ^ LOW_SIGN;
            }
            mergeSort(rows, keys, 0, m, rowsScratch, keysScratch);
            order[start] = rows[0] | SET_START;
            for (int i = 1, setStart = 0; i <= m; i++) {
                if (i == m || !sameDuplicateSet((int) rows[i - 1], (int) rows[i])) {
                    inSets += i - setStart > 1 ? i - setStart : 0;
                    setStart = i;
                    if (i < m) {
                        order[start + i] = rows[i] | SET_START;
                    }
                } else {
                    order[start + i] = rows[i];
                }
            }
        }
        return inSets;
    }

    private int compareInRun(final long rowA, final long keyA, final long rowB, final long keyB) {
        if (wholePrefix && keyA != keyB) {
            return keyA < keyB ? -1 : 1;
        }
        return compareRows((int) rowA, (int) rowB);
    }

    private void mergeSort(final long[] rows, final long[] keys, final int from, final int to,
                           final long[] rowsScratch, final long[] keysScratch) {
        if (to - from <= INSERTION_SORT_MAX) {
            for (int i = from + 1; i < to; i++) {
                final long row = rows[i], key = keys[i];
                int j = i - 1;
                while (j >= from && compareInRun(rows[j], keys[j], row, key) > 0) {
                    rows[j + 1] = rows[j];
                    keys[j + 1] = keys[j];
                    j--;
                }
                rows[j + 1] = row;
                keys[j + 1] = key;
            }
            return;
        }
        final int middle = (from + to) >>> 1;
        mergeSort(rows, keys, from, middle, rowsScratch, keysScratch);
        mergeSort(rows, keys, middle, to, rowsScratch, keysScratch);
        if (compareInRun(rows[middle - 1], keys[middle - 1], rows[middle], keys[middle]) <= 0) {
            return;
        }
        final int leftLength = middle - from;
        System.arraycopy(rows, from, rowsScratch, 0, leftLength);
        System.arraycopy(keys, from, keysScratch, 0, leftLength);
        int left = 0, right = middle, out = from;
        while (left < leftLength && right < to) {
            if (compareInRun(rowsScratch[left], keysScratch[left], rows[right], keys[right]) <= 0) {
                rows[out] = rowsScratch[left];
                keys[out++] = keysScratch[left++];
            } else {
                rows[out] = rows[right];
                keys[out++] = keys[right++];
            }
        }
        System.arraycopy(rowsScratch, left, rows, out, leftLength - left);
        System.arraycopy(keysScratch, left, keys, out, leftLength - left);
    }

    /** ReadEndsMDComparator without barcodes on two rows, then the order they were added. */
    private int compareRows(final int a, final int b) {
        final long[] ca = chunk(a);
        final long[] cb = chunk(b);
        final int ia = at(a);
        final int ib = at(b);
        int c = Short.compare((short) (ca[ia + 5] >>> 48), (short) (cb[ib + 5] >>> 48));
        if (c == 0) {
            c = Long.compare(ca[ia + 2] ^ LOW_SIGN, cb[ib + 2] ^ LOW_SIGN);
        }
        if (c == 0) {
            c = Byte.compare((byte) (ca[ia + 6] >>> 8), (byte) (cb[ib + 6] >>> 8));
        }
        if (c == 0) {
            c = Long.compare(ca[ia + 3] ^ LOW_SIGN, cb[ib + 3] ^ LOW_SIGN);
        }
        if (c == 0) {
            c = Short.compare((short) (ca[ia + 5] >>> 32), (short) (cb[ib + 5] >>> 32));
        }
        if (c == 0) {
            c = Long.compare(ca[ia + 4] ^ LOW_SIGN, cb[ib + 4] ^ LOW_SIGN);
        }
        if (c == 0) {
            c = Long.compare(ca[ia], cb[ib]);
        }
        if (c == 0) {
            c = Long.compare(ca[ia + 1], cb[ib + 1]);
        }
        return c != 0 ? c : Integer.compare(a, b);
    }

    /** The row at the given position of the sorted order. */
    public int rowAt(final int position) {
        return (int) (order[position] & rowMask);
    }

    /**
     * The end at the given position of the sorted order, as get(rowAt(position)) gives it; read in sequence when the
     * positions of the duplicate sets of more than one end are asked for in increasing order.
     */
    public ReadEndsForMarkDuplicates getAt(final int position) {
        while (setCursor < setPositions.length && setPositions[setCursor] < position) {
            setCursor++;
        }
        if (setCursor < setPositions.length && setPositions[setCursor] == position) {
            return toEnd(setEnds, setCursor++ * USED_LONGS);
        }
        return get(rowAt(position));
    }

    /** Whether the row at the given position of the sorted order is not in the duplicate set of the one before. */
    public boolean startsDuplicateSet(final int position) {
        return order[position] < 0;
    }

    /**
     * Whether two rows are comparable for duplicates as pairs, as MarkDuplicates.areComparableForDuplicates with
     * compareRead2 and without barcodes: equal library id, read1 reference index and coordinate, orientation, and
     * read2 reference index and coordinate.
     */
    public boolean sameDuplicateSet(final int a, final int b) {
        final long[] ca = chunk(a);
        final long[] cb = chunk(b);
        final int ia = at(a);
        final int ib = at(b);
        return ca[ia + 2] == cb[ib + 2] && ca[ia + 3] == cb[ib + 3] &&
                (ca[ia + 5] >>> 48) == (cb[ib + 5] >>> 48) &&
                ((ca[ia + 6] ^ cb[ib + 6]) & 0xFF00L) == 0;
    }
}
