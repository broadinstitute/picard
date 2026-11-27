package picard.sam.util;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockGunzipper;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * Iterates over the records of a BAM file in file order as htsjdk's BAMFileReader iterator does, with the records
 * decoded, validated and eagerly decoded on several threads. One thread inflates the file with htsjdk's
 * BlockCompressedInputStream, so inflation and its errors are htsjdk's, and cuts the stream at record boundaries into
 * batches, a pool decodes each batch as BAMRecordCodec does, and the consumer takes the batches in order. The stream is
 * not htsjdk's AsyncBlockCompressedInputStream, which never returns from a read after a corrupt deflate block. An empty
 * BGZF block ends the records where it ends BAMFileReader's: where BAMRecordCodec starts a read of the stream (cleanly
 * at a record's start, with an EOF exception inside one), or right after another empty block; elsewhere it is read past.
 *
 * Each record gets what BAMFileReader gives it: setValidationStringency, then unless SILENT isValid and
 * SAMUtils.processValidationErrors with its 1-based record number, then eagerDecode when asked. As in BAMFileReader the
 * iterator reads one record ahead, so a failure in record k is thrown by the next() that would return record k-1, or by
 * {@link #open} when k is 1, with the same exception. LENIENT validation warnings are printed in file order, when the
 * iterator reaches their record. Records carry no file source (SamReaderFactory's INCLUDE_SOURCE_IN_RECORDS). The
 * records are made by a SAMRecordFactory, as SamReaderFactory.samRecordFactory sets for BAMFileReader. An optional
 * annotator runs on the decoding threads too, once per record, and {@link #annotation()} returns its result.
 */
public final class ParallelBAMRecordIterator implements CloseableIterator<SAMRecord> {
    private static final int FIXED_BLOCK_SIZE = 32;

    // Package-private and not final so that a test can force small batches and the sequential tail. Larger batches
    // or deeper queues, which put more decoded records in flight, make pass 1 slower.
    static int batchBytes = 64 * 1024;
    // A record claiming more than this is read, with every record after it, on the consumer's thread.
    static int maxBufferedRecordBytes = 64 * 1024 * 1024;
    private static final int QUEUED_BATCHES_PER_THREAD = 1;

    private static final AtomicInteger POOL_NUMBER = new AtomicInteger();

    /** Bytes of whole records cut from the decompressed stream, and what ended the cut. */
    private static final class RawBatch {
        final byte[] bytes;
        final int length;
        final int records;
        final long firstIndex;
        final Throwable failure;
        final boolean last;
        final byte[] tail;

        RawBatch(final byte[] bytes, final int length, final int records, final long firstIndex,
                 final Throwable failure, final boolean last, final byte[] tail) {
            this.bytes = bytes;
            this.length = length;
            this.records = records;
            this.firstIndex = firstIndex;
            this.failure = failure;
            this.last = last;
            this.tail = tail;
        }
    }

    private static final class DecodedBatch {
        final SAMRecord[] records;
        // LENIENT validation errors per record, printed when the record is reached; null when there are none.
        final List<SAMValidationError>[] errors;
        // The annotator's result per record; null without an annotator.
        final Object[] annotations;
        final int size;
        final long firstIndex;
        // Thrown once the records are consumed, after failureErrors are printed.
        final Throwable failure;
        final List<SAMValidationError> failureErrors;
        final long failureIndex;
        final boolean last;
        // Not null when the bytes after this batch are read on the consumer's thread, starting with these.
        final byte[] tail;

        DecodedBatch(final SAMRecord[] records, final List<SAMValidationError>[] errors, final Object[] annotations,
                     final int size, final long firstIndex, final Throwable failure,
                     final List<SAMValidationError> failureErrors, final long failureIndex, final boolean last,
                     final byte[] tail) {
            this.records = records;
            this.errors = errors;
            this.annotations = annotations;
            this.size = size;
            this.firstIndex = firstIndex;
            this.failure = failure;
            this.failureErrors = failureErrors;
            this.failureIndex = failureIndex;
            this.last = last;
            this.tail = tail;
        }
    }

    /**
     * The bytes of one batch, then the stored failure, or the rest of another stream. Reads stop at the end of the
     * bytes, so a failure is thrown by the read after the last good byte, within the same BinaryCodec read as
     * BlockCompressedInputStream would throw it.
     */
    private static final class ReplayStream extends InputStream {
        private final byte[] bytes;
        private final int end;
        private final Throwable failure;
        private final InputStream next;
        private int position = 0;

        ReplayStream(final byte[] bytes, final int end, final Throwable failure, final InputStream next) {
            this.bytes = bytes;
            this.end = end;
            this.failure = failure;
            this.next = next;
        }

        @Override
        public int read() throws IOException {
            if (position < end) {
                return bytes[position++] & 0xff;
            }
            return next == null ? endOfBytes() : next.read();
        }

        @Override
        public int read(final byte[] buffer, final int offset, final int length) throws IOException {
            if (length == 0) {
                return 0;
            }
            if (position < end) {
                final int n = Math.min(length, end - position);
                System.arraycopy(bytes, position, buffer, offset, n);
                position += n;
                if (n == length || next == null) {
                    return n;
                }
                // One read across the hand-over, so an empty block there is read past as within one stream read.
                final int more = next.read(buffer, offset + n, length - n);
                return more < 0 ? n : n + more;
            }
            return next == null ? endOfBytes() : next.read(buffer, offset, length);
        }

        private int endOfBytes() throws IOException {
            if (failure instanceof IOException) {
                throw (IOException) failure;
            }
            if (failure != null) {
                throw rethrow(failure);
            }
            return -1;
        }
    }

    /**
     * The number of records, from the first, that another read of the same file validated at the same stringency
     * without a failure, while that read is still under way: it raises the count as it goes and gives the final count
     * when it ends. Decoding a record past the count waits until the record is known to be among those validated or
     * not. {@link #abort} ends every wait, and the records still waiting are then validated.
     */
    public static final class ValidatedRecords {
        private final Object lock = new Object();
        private volatile long known = 0;
        private volatile long total = -1;
        private volatile boolean aborted = false;
        // Guarded by lock.
        private int waiters = 0;

        /** At least count records are validated; called by one thread. */
        public void advance(final long count) {
            if (count > known) {
                known = count;
                wake();
            }
        }

        /** Exactly count records are validated. */
        public void finish(final long count) {
            total = count;
            wake();
        }

        /** Ends every wait; records not yet known to be validated are validated. */
        public void abort() {
            aborted = true;
            wake();
        }

        private void wake() {
            synchronized (lock) {
                if (waiters > 0) {
                    lock.notifyAll();
                }
            }
        }

        /** Whether the record with this 1-based number is among the validated ones, once that is known. */
        boolean covers(final long index) {
            if (index <= known) {
                return true;
            }
            synchronized (lock) {
                waiters++;
                try {
                    while (true) {
                        if (index <= known) {
                            return true;
                        }
                        final long count = total;
                        if (count >= 0) {
                            return index <= count;
                        }
                        if (aborted) {
                            return false;
                        }
                        try {
                            lock.wait();
                        } catch (final InterruptedException e) {
                            // Only release() interrupts a decoder, and its batch is dropped.
                            Thread.currentThread().interrupt();
                            return false;
                        }
                    }
                } finally {
                    waiters--;
                }
            }
        }
    }

    private final SAMFileHeader header;
    private final SAMRecordFactory recordFactory;
    private final ValidationStringency stringency;
    // Records 1 to validatedRecords are not validated again; validatedCount, when not null, gives that count instead.
    private final long validatedRecords;
    private final ValidatedRecords validatedCount;
    private final boolean eagerDecode;
    private final String fileName;
    private final BlockCompressedInputStream stream;
    private final ExecutorService decoders;
    private final BlockingQueue<Future<DecodedBatch>> batches;
    // Each decoding thread's annotator; null without one.
    private final ThreadLocal<Function<SAMRecord, Object>> annotators;
    private final Thread splitter;
    private volatile boolean stopping = false;

    private DecodedBatch current;
    private int position = 0;
    private BAMRecordCodec tailCodec = null;
    private long tailIndex = 0;
    private SAMRecord nextRecord = null;
    private Object pulledAnnotation = null;
    private Object nextAnnotation = null;
    private Object annotation = null;
    private boolean released = false;
    private boolean closed = false;

    /**
     * Opens an iterator over every record of a BAM file, from its first record.
     *
     * @param bam the BAM file, read through a stream of its own
     * @param header the header the records get, the one the file's SamReader returned
     * @param stringency the SamReader's validation stringency
     * @param eagerDecode whether each record is eagerly decoded, as SamReaderFactory.Option.EAGERLY_DECODE does
     * @param threads the number of decoding threads, at least 1
     * @param threadName the prefix of the threads' names
     */
    public static ParallelBAMRecordIterator open(final File bam, final SAMFileHeader header,
                                                 final ValidationStringency stringency, final boolean eagerDecode,
                                                 final int threads, final String threadName) {
        return open(bam, header, stringency, eagerDecode, threads, threadName, DefaultSAMRecordFactory.getInstance());
    }

    /**
     * As {@link #open(File, SAMFileHeader, ValidationStringency, boolean, int, String)}, with the records made by
     * recordFactory, which the decoding threads share.
     */
    public static ParallelBAMRecordIterator open(final File bam, final SAMFileHeader header,
                                                 final ValidationStringency stringency, final boolean eagerDecode,
                                                 final int threads, final String threadName,
                                                 final SAMRecordFactory recordFactory) {
        return open(bam, header, stringency, eagerDecode, threads, threadName, recordFactory, null);
    }

    /**
     * As {@link #open(File, SAMFileHeader, ValidationStringency, boolean, int, String, SAMRecordFactory)}, with each
     * record passed, on its decoding thread and after its validation and eager decoding, to an annotator that
     * annotatorFactory makes once per decoding thread; {@link #annotation()} returns the result. A record whose
     * annotator throws, and every record read on the consumer's thread, has a null annotation.
     *
     * @param annotatorFactory makes one decoding thread's annotator; null for none
     */
    public static ParallelBAMRecordIterator open(final File bam, final SAMFileHeader header,
                                                 final ValidationStringency stringency, final boolean eagerDecode,
                                                 final int threads, final String threadName,
                                                 final SAMRecordFactory recordFactory,
                                                 final Supplier<Function<SAMRecord, Object>> annotatorFactory) {
        return open(bam, header, stringency, eagerDecode, threads, threadName, recordFactory, annotatorFactory, 0);
    }

    /**
     * As {@link #open(File, SAMFileHeader, ValidationStringency, boolean, int, String, SAMRecordFactory, Supplier)},
     * for a file whose first validatedRecords records an earlier read of it validated and eagerly decoded at
     * stringency without a failure: those are not validated again. Their tags and cigar are decoded at SILENT on the
     * decoding threads, the values validation would decode, or for a RawBlockBAMRecord left for it to decode at
     * SILENT when first needed, and they are handed out at stringency. Every later record
     * is validated as BAMFileReader validates it.
     *
     * @param validatedRecords the number of records, from the first, already validated at stringency
     */
    public static ParallelBAMRecordIterator open(final File bam, final SAMFileHeader header,
                                                 final ValidationStringency stringency, final boolean eagerDecode,
                                                 final int threads, final String threadName,
                                                 final SAMRecordFactory recordFactory,
                                                 final Supplier<Function<SAMRecord, Object>> annotatorFactory,
                                                 final long validatedRecords) {
        return open(bam, header, stringency, eagerDecode, threads, threadName, recordFactory, annotatorFactory,
                validatedRecords, null);
    }

    /**
     * As {@link #open(File, SAMFileHeader, ValidationStringency, boolean, int, String, SAMRecordFactory, Supplier,
     * long)}, with the number of validated records given by validatedCount as another read of the file goes on.
     */
    public static ParallelBAMRecordIterator open(final File bam, final SAMFileHeader header,
                                                 final ValidationStringency stringency, final boolean eagerDecode,
                                                 final int threads, final String threadName,
                                                 final SAMRecordFactory recordFactory,
                                                 final Supplier<Function<SAMRecord, Object>> annotatorFactory,
                                                 final ValidatedRecords validatedCount) {
        return open(bam, header, stringency, eagerDecode, threads, threadName, recordFactory, annotatorFactory, 0,
                validatedCount);
    }

    private static ParallelBAMRecordIterator open(final File bam, final SAMFileHeader header,
                                                  final ValidationStringency stringency, final boolean eagerDecode,
                                                  final int threads, final String threadName,
                                                  final SAMRecordFactory recordFactory,
                                                  final Supplier<Function<SAMRecord, Object>> annotatorFactory,
                                                  final long validatedRecords, final ValidatedRecords validatedCount) {
        if (threads < 1) {
            throw new IllegalArgumentException("threads must be at least 1: " + threads);
        }
        final long firstRecord;
        try (SeekableStream headerStream = new SeekableFileStream(bam)) {
            firstRecord = SAMUtils.findVirtualOffsetOfFirstRecordInBam(headerStream);
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
        final BlockCompressedInputStream stream;
        try {
            stream = new BlockCompressedInputStream(bam, BlockGunzipper.getDefaultInflaterFactory());
            stream.seek(firstRecord);
        } catch (final IOException e) {
            // BAMFileReader.getIterator's message
            throw new RuntimeIOException(e.getMessage(), e);
        }
        final ParallelBAMRecordIterator iterator = new ParallelBAMRecordIterator(stream, header, stringency,
                eagerDecode, bam.getAbsolutePath(), threads, threadName, recordFactory, annotatorFactory,
                validatedRecords, validatedCount);
        try {
            iterator.advance();
        } catch (final RuntimeException | Error e) {
            iterator.release();
            throw e;
        }
        return iterator;
    }

    private ParallelBAMRecordIterator(final BlockCompressedInputStream stream, final SAMFileHeader header,
                                      final ValidationStringency stringency, final boolean eagerDecode,
                                      final String fileName, final int threads, final String threadName,
                                      final SAMRecordFactory recordFactory,
                                      final Supplier<Function<SAMRecord, Object>> annotatorFactory,
                                      final long validatedRecords, final ValidatedRecords validatedCount) {
        this.stream = stream;
        this.validatedRecords = validatedRecords;
        this.validatedCount = validatedCount;
        this.annotators = annotatorFactory == null ? null : ThreadLocal.withInitial(annotatorFactory);
        this.recordFactory = recordFactory;
        this.header = header;
        this.stringency = stringency;
        this.eagerDecode = eagerDecode;
        this.fileName = fileName;
        this.batches = new ArrayBlockingQueue<>(QUEUED_BATCHES_PER_THREAD * threads);
        final String poolName = threadName + "-" + POOL_NUMBER.incrementAndGet();
        final AtomicInteger decoderNumber = new AtomicInteger();
        this.decoders = Executors.newFixedThreadPool(threads, r -> {
            final Thread thread = new Thread(r, poolName + "-decoder-" + decoderNumber.incrementAndGet());
            thread.setDaemon(true);
            return thread;
        });
        this.current = new DecodedBatch(new SAMRecord[0], null, null, 0, 0, null, null, 0, false, null);
        this.splitter = new Thread(this::split, poolName + "-reader");
        this.splitter.setDaemon(true);
        this.splitter.start();
    }

    private static int intAt(final byte[] bytes, final int offset) {
        return (bytes[offset] & 0xff) | (bytes[offset + 1] & 0xff) << 8 | (bytes[offset + 2] & 0xff) << 16 |
                (bytes[offset + 3] & 0xff) << 24;
    }

    /** Whether BAMRecordCodec.decode starts a read of the stream at this offset into a record. */
    private static boolean startsRead(final int offset) {
        switch (offset) {
            case 0: case 4: case 8: case 12: case 13: case 14: case 16: case 18: case 20: case 24: case 28: case 32:
            case 36:
                return true;
            default:
                return false;
        }
    }

    /** Cuts the decompressed stream into batches of whole records and queues each one's decoding, in order. */
    private void split() {
        try {
            long index = 0;
            byte[] buffer = new byte[batchBytes];
            int length = 0;
            // Decompressed bytes read so far, and their count when an empty block was last read past.
            long consumed = 0;
            long emptyAt = -1;
            while (!stopping) {
                int scanned = 0;
                int records = 0;
                Throwable failure = null;
                boolean end = false;
                boolean invalid = false;
                boolean tail = false;
                while (true) {
                    while (length - scanned >= 4) {
                        final int recordLength = intAt(buffer, scanned);
                        if (recordLength < FIXED_BLOCK_SIZE) {
                            // BAMRecordCodec throws on this length; give it nothing after it.
                            invalid = true;
                            break;
                        }
                        if (recordLength > length - scanned - 4) {
                            break;
                        }
                        scanned += 4 + recordLength;
                        records++;
                    }
                    if (invalid) {
                        break;
                    }
                    if (length == buffer.length) {
                        if (scanned > 0) {
                            break;
                        }
                        final long needed = 4L + intAt(buffer, 0);
                        if (needed > maxBufferedRecordBytes) {
                            tail = true;
                            break;
                        }
                        buffer = Arrays.copyOf(buffer, (int) needed);
                    }
                    // Reads stay within one BGZF block, so a block that fails to read loses no bytes before it.
                    try {
                        final int available = stream.available();
                        if (available == 0) {
                            // An empty block, or the end of the file. BAMFileReader's stream returns -1 to a read that
                            // meets one before its first byte and reads on past it otherwise.
                            if (consumed == emptyAt || startsRead(length - scanned)) {
                                end = true;
                                break;
                            }
                            emptyAt = consumed;
                            continue;
                        }
                        final int n = stream.read(buffer, length, Math.min(available, buffer.length - length));
                        length += n;
                        consumed += n;
                    } catch (final Throwable t) {
                        failure = t;
                        end = true;
                        break;
                    }
                    if (stopping) {
                        return;
                    }
                }

                final RawBatch raw;
                byte[] next = null;
                if (invalid) {
                    raw = new RawBatch(buffer, scanned + 4, records, index, null, true, null);
                } else if (end) {
                    raw = new RawBatch(buffer, length, records, index, failure, true, null);
                } else if (tail) {
                    raw = new RawBatch(buffer, scanned, records, index, null, false,
                            Arrays.copyOfRange(buffer, scanned, length));
                } else {
                    raw = new RawBatch(buffer, scanned, records, index, null, false, null);
                    next = new byte[batchBytes];
                    System.arraycopy(buffer, scanned, next, 0, length - scanned);
                }
                if (!queue(decoders.submit(() -> decode(raw)))) {
                    return;
                }
                if (raw.last || raw.tail != null) {
                    return;
                }
                index += records;
                length -= scanned;
                buffer = next;
            }
        } catch (final Throwable t) {
            queue(CompletableFuture.completedFuture(
                    new DecodedBatch(new SAMRecord[0], null, null, 0, 0, t, null, 0, true, null)));
        }
    }

    /** Puts a batch in the queue; false if the iterator is being released. */
    private boolean queue(final Future<DecodedBatch> batch) {
        boolean interrupted = false;
        try {
            while (!stopping) {
                try {
                    batches.put(batch);
                    return true;
                } catch (final InterruptedException e) {
                    interrupted = true;
                }
            }
            return false;
        } finally {
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
        }
    }

    /**
     * BAMRecordCodec.decode of the whole record at offset: its fixed fields, their conversions and the factory call.
     * The splitter has checked that block_size is at least 32 and that the record ends within the bytes.
     */
    private SAMRecord createRecord(final byte[] bytes, final int offset) {
        return recordFactory.createBAMRecord(header, intAt(bytes, offset + 4), intAt(bytes, offset + 8) + 1,
                (short) (bytes[offset + 12] & 0xff), (short) (bytes[offset + 13] & 0xff), ushortAt(bytes, offset + 14),
                ushortAt(bytes, offset + 16), ushortAt(bytes, offset + 18), intAt(bytes, offset + 20),
                intAt(bytes, offset + 24), intAt(bytes, offset + 28) + 1, intAt(bytes, offset + 32),
                Arrays.copyOfRange(bytes, offset + 36, offset + 4 + intAt(bytes, offset)));
    }

    private static int ushortAt(final byte[] bytes, final int offset) {
        return (bytes[offset] & 0xff) | (bytes[offset + 1] & 0xff) << 8;
    }

    @SuppressWarnings("unchecked")
    private DecodedBatch decode(final RawBatch raw) {
        // A batch that is not the last holds raw.records whole records and no failure; the last one may end in a
        // partial record, an invalid length or a stored failure, which the codec reports as BAMFileReader's does.
        final BAMRecordCodec codec;
        if (raw.last) {
            codec = new BAMRecordCodec(header, recordFactory);
            codec.setInputStream(new ReplayStream(raw.bytes, raw.length, raw.failure, null), fileName);
        } else {
            codec = null;
        }
        int offset = 0;
        SAMRecord[] records = new SAMRecord[raw.records];
        final Function<SAMRecord, Object> annotator = annotators == null ? null : annotators.get();
        Object[] annotations = annotator == null ? null : new Object[raw.records];
        List<SAMValidationError>[] errors = null;
        List<SAMValidationError> recordErrors = null;
        long index = raw.firstIndex;
        int size = 0;
        Throwable failure = null;
        try {
            while (true) {
                recordErrors = null;
                final SAMRecord record;
                if (codec != null) {
                    record = codec.decode();
                    if (record == null) {
                        break;
                    }
                } else {
                    if (size == raw.records) {
                        break;
                    }
                    record = createRecord(raw.bytes, offset);
                    offset += 4 + intAt(raw.bytes, offset);
                }
                ++index;
                if (validatedCount == null ? index <= validatedRecords : validatedCount.covers(index)) {
                    decodeValidated(record);
                } else {
                    record.setValidationStringency(stringency);
                    if (stringency != ValidationStringency.SILENT) {
                        recordErrors = record.isValid(stringency == ValidationStringency.STRICT);
                        if (recordErrors != null && !recordErrors.isEmpty() && stringency == ValidationStringency.STRICT) {
                            SAMUtils.processValidationErrors(recordErrors, index, stringency);
                        }
                    }
                }
                if (eagerDecode) {
                    eagerDecode(record);
                }
                if (size == records.length) {
                    records = Arrays.copyOf(records, Math.max(1, 2 * size));
                    if (errors != null) {
                        errors = Arrays.copyOf(errors, records.length);
                    }
                    if (annotations != null) {
                        annotations = Arrays.copyOf(annotations, records.length);
                    }
                }
                if (annotator != null) {
                    annotations[size] = annotate(annotator, record);
                }
                if (recordErrors != null && !recordErrors.isEmpty()) {
                    if (errors == null) {
                        errors = new List[records.length];
                    }
                    errors[size] = recordErrors;
                }
                records[size++] = record;
            }
        } catch (final Throwable t) {
            failure = t;
        }
        // LENIENT warnings of a record whose eager decoding then failed; at STRICT the failure is the warning.
        final List<SAMValidationError> failureErrors = failure != null && stringency == ValidationStringency.LENIENT &&
                recordErrors != null && !recordErrors.isEmpty() ? recordErrors : null;
        return new DecodedBatch(records, errors, annotations, size, raw.firstIndex, failure, failureErrors, index,
                raw.last, raw.tail);
    }

    private static Object annotate(final Function<SAMRecord, Object> annotator, final SAMRecord record) {
        try {
            return annotator.apply(record);
        } catch (final Throwable t) {
            return null;
        }
    }

    private DecodedBatch take() {
        final Future<DecodedBatch> batch;
        try {
            batch = batches.take();
        } catch (final InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted waiting for a decoded batch", e);
        }
        boolean interrupted = false;
        try {
            while (true) {
                try {
                    return batch.get();
                } catch (final InterruptedException e) {
                    interrupted = true;
                } catch (final ExecutionException e) {
                    return new DecodedBatch(new SAMRecord[0], null, null, 0, 0, e.getCause(), null, 0, true, null);
                }
            }
        } finally {
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
        }
    }

    /** The next record in file order, null at the end; prints its LENIENT warnings and throws its failure. */
    private SAMRecord pull() {
        while (true) {
            if (tailCodec != null) {
                return pullTail();
            }
            if (position < current.size) {
                final SAMRecord record = current.records[position];
                current.records[position] = null;
                if (current.annotations != null) {
                    pulledAnnotation = current.annotations[position];
                    current.annotations[position] = null;
                }
                if (current.errors != null && current.errors[position] != null) {
                    SAMUtils.processValidationErrors(current.errors[position], current.firstIndex + position + 1,
                            stringency);
                }
                position++;
                return record;
            }
            if (current.failure != null) {
                final DecodedBatch failed = current;
                release();
                if (failed.failureErrors != null) {
                    SAMUtils.processValidationErrors(failed.failureErrors, failed.failureIndex, stringency);
                }
                throw rethrow(failed.failure);
            }
            if (current.tail != null) {
                tailCodec = new BAMRecordCodec(header, recordFactory);
                tailCodec.setInputStream(new ReplayStream(current.tail, current.tail.length, null, stream), fileName);
                tailIndex = current.firstIndex + current.size;
                continue;
            }
            if (current.last) {
                release();
                return null;
            }
            current = take();
            position = 0;
        }
    }

    /** What BAMFileReader's iterator does for one record. */
    private SAMRecord pullTail() {
        final SAMRecord record = tailCodec.decode();
        if (record == null) {
            release();
            return null;
        }
        ++tailIndex;
        if (validatedCount == null ? tailIndex <= validatedRecords : validatedCount.covers(tailIndex)) {
            decodeValidated(record);
        } else {
            record.setValidationStringency(stringency);
            if (stringency != ValidationStringency.SILENT) {
                SAMUtils.processValidationErrors(record.isValid(stringency == ValidationStringency.STRICT), tailIndex,
                        stringency);
            }
        }
        if (eagerDecode) {
            eagerDecode(record);
        }
        return record;
    }

    /**
     * For a record already validated at {@link #stringency}: decodes what isValid would, its tags and then its cigar,
     * in isValid's order, at SILENT, so that getCigar does not check the cigar again; then sets the stringency. A
     * record's stringency decides only whether decoding reports a problem, never the decoded value. A
     * RawBlockBAMRecord whose tags {@link RawBlockBAMRecord#deferTags} accepts decodes them when first needed instead.
     */
    private void decodeValidated(final SAMRecord record) {
        if (!(record instanceof RawBlockBAMRecord) || !((RawBlockBAMRecord) record).deferTags()) {
            record.setValidationStringency(ValidationStringency.SILENT);
            // Any tag: BAMRecord.getAttribute decodes them all first.
            record.getAttribute((short) 0);
            record.getCigar();
        }
        record.setValidationStringency(stringency);
    }

    /**
     * BAMRecord.eagerDecode, which is protected, through public getters in its order, so the same decoding fails first.
     * Unlike it this keeps the record's undecoded bytes, which no getter reads once everything is decoded.
     */
    private static void eagerDecode(final SAMRecord record) {
        record.getReadName();
        record.getCigar();
        record.getReadBases();
        record.getBaseQualities();
        // Any tag: BAMRecord.getAttribute decodes them all first.
        record.getAttribute((short) 0);
        record.getCigarString();
    }

    private void advance() {
        pulledAnnotation = null;
        nextRecord = pull();
        nextAnnotation = pulledAnnotation;
    }

    private static RuntimeException rethrow(final Throwable failure) {
        if (failure instanceof RuntimeException) {
            throw (RuntimeException) failure;
        }
        if (failure instanceof Error) {
            throw (Error) failure;
        }
        throw new RuntimeException(failure);
    }

    @Override
    public boolean hasNext() {
        if (closed) {
            throw new IllegalStateException("iterator has been closed");
        }
        return nextRecord != null;
    }

    @Override
    public SAMRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        final SAMRecord result = nextRecord;
        annotation = nextAnnotation;
        advance();
        return result;
    }

    /** The annotator's result for the record the last {@link #next()} returned; null when it has none. */
    public Object annotation() {
        return annotation;
    }

    /** Stops the reading thread and the decoders and closes the stream; the records already read stay readable. */
    private void release() {
        if (released) {
            return;
        }
        released = true;
        stopping = true;
        batches.clear();
        boolean interrupted = false;
        while (true) {
            try {
                splitter.join();
                break;
            } catch (final InterruptedException e) {
                interrupted = true;
            }
        }
        if (interrupted) {
            Thread.currentThread().interrupt();
        }
        batches.clear();
        decoders.shutdownNow();
        try {
            stream.close();
        } catch (final IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    @Override
    public void close() {
        if (closed) {
            return;
        }
        closed = true;
        nextRecord = null;
        nextAnnotation = null;
        annotation = null;
        current = null;
        release();
    }
}
