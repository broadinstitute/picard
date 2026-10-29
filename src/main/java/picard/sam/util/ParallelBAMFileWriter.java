package picard.sam.util;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.BinaryTagCodec;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.GenomicIndexUtil;
import htsjdk.samtools.SAMBinaryTagAndValue;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

/**
 * A BAM writer whose BGZF blocks are deflated on several threads by {@link ParallelBlockCompressedOutputStream}.
 * Its bytes are those of htsjdk's BAMFileWriter without an index: the header and the first records share the first
 * block, and records are encoded as {@link BAMRecordCodec} encodes them. A {@link RawBlockBAMRecord} that still holds
 * its block has its read name, cigar, bases and qualities copied from the block, and its deferred tags written from it.
 */
public class ParallelBAMFileWriter extends SAMFileWriterImpl {
    private static final byte[] BAM_MAGIC = {'B', 'A', 'M', 1};
    // A cap for machines with many processors: one thread encodes every record for all the deflate threads.
    private static final int MAX_DEFAULT_THREADS = 8;

    // BAM_FIXED_BLOCK_SIZE plus the block_size field itself.
    private static final int FIXED_FIELDS_LENGTH = 36;
    // Reference-consuming cigar operators by their BAM code: M, D, N, = and X.
    private static final int REFERENCE_CONSUMING_OPS = 1 << 0 | 1 << 2 | 1 << 3 | 1 << 7 | 1 << 8;

    private final BinaryCodec outputBinaryCodec;
    private final ParallelBlockCompressedOutputStream blockStream;
    private BAMRecordCodec bamRecordCodec = null;
    private final RecordBuffer recordBuffer = new RecordBuffer();
    private final BinaryTagCodec tagCodec = new BinaryTagCodec(new BinaryCodec(recordBuffer));
    // The last header and reference name binCompatible looked up, and the answer.
    private SAMFileHeader binHeader = null;
    private String binReferenceName = null;
    private boolean binCompatible = false;

    public ParallelBAMFileWriter(final OutputStream os, final String absoluteFilename, final int compressionLevel,
                                 final DeflaterFactory deflaterFactory, final int threads) {
        blockStream = new ParallelBlockCompressedOutputStream(os, compressionLevel, deflaterFactory, threads);
        outputBinaryCodec = new BinaryCodec(blockStream);
        outputBinaryCodec.setOutputFileName(absoluteFilename);
    }

    /**
     * The number of deflate threads {@link #makeAsyncWriter} uses: the available processors less a quarter of them,
     * rounded down, at most 8. The quarter is left to the threads that read, decode and encode the records.
     */
    public static int defaultThreads() {
        final int processors = Runtime.getRuntime().availableProcessors();
        return Math.max(1, Math.min(processors - processors / 4, MAX_DEFAULT_THREADS));
    }

    /**
     * Whether {@link #makeAsyncWriter} writes this output as {@link SAMFileWriterFactory#makeWriter} would: a file
     * named .bam, with the factory's defaults creating neither an index nor an MD5 file.
     */
    public static boolean canWrite(final File output) {
        if (output == null || SAMFileWriterFactory.getDefaultCreateIndexWhileWriting() || SAMFileWriterFactory.getDefaultCreateMd5File()) {
            return false;
        }
        final Path fileName = IOUtil.toPath(output).getFileName();
        return fileName != null && SamReader.Type.BAM_TYPE.hasValidFileExtension(fileName.toString());
    }

    /**
     * The writer {@code new SAMFileWriterFactory().setUseAsyncIo(true).makeWriter(header, presorted, output, null)}
     * returns for an output {@link #canWrite} accepts, with BGZF blocks deflated on {@link #defaultThreads()} threads.
     * Records are encoded on a writer thread, so none may be changed after addAlignment.
     */
    public static SAMFileWriter makeAsyncWriter(final SAMFileHeader header, final boolean presorted, final File output) {
        return makeAsyncWriter(header, presorted, output, false);
    }

    /**
     * As {@link #makeAsyncWriter(SAMFileHeader, boolean, File)}; with createNew the file is created and must not
     * exist. Opening an existing file truncates it, after which some file systems (ext4's auto_da_alloc) write the
     * whole file out to disk when it is closed.
     */
    public static SAMFileWriter makeAsyncWriter(final SAMFileHeader header, final boolean presorted, final File output,
                                                final boolean createNew) {
        final Path path = IOUtil.toPath(output);
        try {
            final OutputStream os = IOUtil.maybeBufferOutputStream(createNew ?
                    Files.newOutputStream(path, StandardOpenOption.CREATE_NEW, StandardOpenOption.WRITE) :
                    Files.newOutputStream(path), Defaults.BUFFER_SIZE);
            final ParallelBAMFileWriter writer = new ParallelBAMFileWriter(os, path.toUri().toString(),
                    BlockCompressedOutputStream.getDefaultCompressionLevel(),
                    BlockCompressedOutputStream.getDefaultDeflaterFactory(), defaultThreads());
            writer.setSortOrder(header.getSortOrder(), presorted);
            writer.setHeader(header);
            return new AsyncWriter(writer);
        } catch (final IOException e) {
            throw new RuntimeIOException("Error opening file: " + path.toUri(), e);
        }
    }

    /**
     * For a writer {@link #makeAsyncWriter} made: keeps at most this many BGZF blocks submitted for deflating and not
     * yet written, from the next block on, or with a limit below 1 as many as at first
     * (ParallelBlockCompressedOutputStream.setBlocksInFlightLimit). Does nothing to any other writer.
     */
    public static void setBlocksInFlightLimit(final SAMFileWriter writer, final int limit) {
        if (writer instanceof AsyncWriter && ((AsyncWriter) writer).underlyingWriter instanceof ParallelBAMFileWriter) {
            ((ParallelBAMFileWriter) ((AsyncWriter) writer).underlyingWriter).blockStream.setBlocksInFlightLimit(limit);
        }
    }

    @Override
    protected void writeAlignment(final SAMRecord alignment) {
        if (alignment instanceof RawBlockBAMRecord && writeFromBlock((RawBlockBAMRecord) alignment)) {
            return;
        }
        if (bamRecordCodec == null) {
            bamRecordCodec = new BAMRecordCodec(getFileHeader());
            bamRecordCodec.setOutputStream(outputBinaryCodec.getOutputStream(), getFilename());
        }
        bamRecordCodec.encode(alignment);
    }

    /**
     * Writes the bytes BAMRecordCodec.encode writes for this record, copying the read name, cigar, bases and
     * qualities from the record's block, in one write. Returns false, having written nothing, for a record
     * BAMRecordCodec must encode: one without its block, one BAMRecordCodec writes from its block unchanged, one whose
     * cigar came from a CG tag, one whose block or fields BAMRecordCodec would reject or treat specially.
     */
    private boolean writeFromBlock(final RawBlockBAMRecord record) {
        final byte[] block = record.getOriginalBlock();
        if (block == null || record.getVariableBinaryRepresentation() != null) {
            return false;
        }
        // BAMRecordCodec.encode decodes the cigar first, so a cigar that cannot be decoded throws here as it does there.
        // Deferred tags hold no CG tag and a cigar that decodes, so that cigar is the block's.
        final boolean deferredTags = record.hasDeferredTags();
        final int cigarLength = record.getCigarLength();
        if (!deferredTags && record.getCigar().numCigarElements() != cigarLength) {
            return false;
        }
        final int readNameLength = record.getReadNameLength() + 1;
        final int readLength = record.getReadLength();
        final int mappingQuality = record.getMappingQuality();
        final int flags = record.getFlags();
        if (readNameLength < 1 || readLength < 0 || (mappingQuality & ~0xFF) != 0 || (flags & ~0xFFFF) != 0) {
            return false;
        }
        final int cigarOffset = readNameLength;
        final long basesOffset = cigarOffset + 4L * cigarLength;
        final long qualitiesOffset = basesOffset + (readLength + 1L) / 2;
        final long tagsOffset = qualitiesOffset + readLength;
        if (tagsOffset > block.length) {
            return false;
        }

        // BAMRecordCodec.encode's bin: 0 when unplaced or on a reference too long for the BAI scheme, otherwise
        // SAMRecord.computeIndexingBin with the reference length summed from the block's cigar.
        final int alignmentStart = record.getAlignmentStart();
        int bin = 0;
        if (alignmentStart != 0) {
            if (!binCompatible(record)) {
                return false;
            }
            final int start = alignmentStart - 1;
            int end = record.getReadUnmappedFlag() ? 0 : alignmentStart + referenceLength(block, cigarOffset, cigarLength) - 1;
            if (end <= 0) {
                end = start + 1;
            }
            if (start > GenomicIndexUtil.BIN_GENOMIC_SPAN || end > GenomicIndexUtil.BIN_GENOMIC_SPAN) {
                return false;
            }
            bin = GenomicIndexUtil.regionToBin(start, end) & 0xFFFF;
        }

        final int variableLength = (int) tagsOffset;
        final byte[] bytes = recordBuffer.reserve(FIXED_FIELDS_LENGTH + variableLength);
        putInt(bytes, 4, record.getReferenceIndex());
        putInt(bytes, 8, alignmentStart - 1);
        bytes[12] = (byte) readNameLength;
        bytes[13] = (byte) mappingQuality;
        putShort(bytes, 14, bin);
        putShort(bytes, 16, cigarLength);
        putShort(bytes, 18, flags);
        putInt(bytes, 20, readLength);
        putInt(bytes, 24, record.getMateReferenceIndex());
        putInt(bytes, 28, record.getMateAlignmentStart() - 1);
        putInt(bytes, 32, record.getInferredInsertSize());

        // The name, then a NUL whatever the block's terminator byte is; then cigar, bases and qualities as they are.
        System.arraycopy(block, 0, bytes, FIXED_FIELDS_LENGTH, readNameLength - 1);
        bytes[FIXED_FIELDS_LENGTH + readNameLength - 1] = 0;
        System.arraycopy(block, cigarOffset, bytes, FIXED_FIELDS_LENGTH + cigarOffset, variableLength - cigarOffset);
        final int qualities = FIXED_FIELDS_LENGTH + (int) qualitiesOffset;
        if ((readLength & 1) != 0) {
            // Re-packing decoded bases leaves the last byte's padding nybble zero.
            bytes[qualities - 1] &= (byte) 0xF0;
        }
        if (readLength > 0 && block[(int) qualitiesOffset] == (byte) 0xFF) {
            // Qualities starting with 0xFF decode as missing, and missing qualities are written as all 0xFF.
            Arrays.fill(bytes, qualities, qualities + readLength, (byte) 0xFF);
        }

        recordBuffer.length = FIXED_FIELDS_LENGTH + variableLength;
        if (deferredTags) {
            final byte[] out = recordBuffer.reserve(recordBuffer.length + record.maxTagsLength());
            recordBuffer.length = record.writeTags(out, recordBuffer.length);
        } else {
            for (SAMBinaryTagAndValue tag = record.getTags(); tag != null; tag = tag.getNext()) {
                tagCodec.writeTag(tag.tag, tag.value, tag.isUnsignedArray());
            }
        }
        putInt(recordBuffer.bytes, 0, recordBuffer.length - 4);
        outputBinaryCodec.writeBytes(recordBuffer.bytes, 0, recordBuffer.length);
        return true;
    }

    /** Whether BAMRecordCodec computes this record's bin: false for a reference too long for the BAI scheme. */
    private boolean binCompatible(final SAMRecord record) {
        final SAMFileHeader header = record.getHeader();
        final String referenceName = record.getReferenceName();
        if (header != binHeader || referenceName != binReferenceName) {
            final SAMSequenceRecord sequence = header == null ? null : header.getSequence(referenceName);
            binCompatible = sequence == null || !SAMUtils.isReferenceSequenceIncompatibleWithBAI(sequence);
            binHeader = header;
            binReferenceName = referenceName;
        }
        return binCompatible;
    }

    /** Cigar.getReferenceLength of a binary cigar whose operators all decode. */
    private static int referenceLength(final byte[] block, final int offset, final int cigarLength) {
        int length = 0;
        for (int i = 0, position = offset; i < cigarLength; i++, position += 4) {
            final int element = (block[position] & 0xFF) | (block[position + 1] & 0xFF) << 8 |
                    (block[position + 2] & 0xFF) << 16 | (block[position + 3] & 0xFF) << 24;
            if ((REFERENCE_CONSUMING_OPS >>> (element & 0xF) & 1) != 0) {
                length += element >>> 4;
            }
        }
        return length;
    }

    private static void putShort(final byte[] bytes, final int offset, final int value) {
        bytes[offset] = (byte) value;
        bytes[offset + 1] = (byte) (value >>> 8);
    }

    private static void putInt(final byte[] bytes, final int offset, final int value) {
        bytes[offset] = (byte) value;
        bytes[offset + 1] = (byte) (value >>> 8);
        bytes[offset + 2] = (byte) (value >>> 16);
        bytes[offset + 3] = (byte) (value >>> 24);
    }

    /** BAMFileWriter.writeHeader(BinaryCodec, SAMFileHeader, String). */
    @Override
    protected void writeHeader(final String textHeader) {
        outputBinaryCodec.writeBytes(BAM_MAGIC);
        outputBinaryCodec.writeString(textHeader, true, false);
        outputBinaryCodec.writeInt(getFileHeader().getSequenceDictionary().size());
        for (final SAMSequenceRecord sequenceRecord : getFileHeader().getSequenceDictionary().getSequences()) {
            outputBinaryCodec.writeString(sequenceRecord.getSequenceName(), true, true);
            outputBinaryCodec.writeInt(sequenceRecord.getSequenceLength());
        }
    }

    @Override
    protected void finish() {
        outputBinaryCodec.close();
    }

    @Override
    protected String getFilename() {
        return outputBinaryCodec.getOutputFileName();
    }

    /** One record's bytes, collected so they reach the BGZF stream in one write. */
    private static final class RecordBuffer extends OutputStream {
        byte[] bytes = new byte[1024];
        int length = 0;

        /** The buffer, grown to hold at least capacity bytes; its contents up to length are kept. */
        byte[] reserve(final int capacity) {
            if (bytes.length < capacity) {
                bytes = Arrays.copyOf(bytes, Math.max(capacity, 2 * bytes.length));
            }
            return bytes;
        }

        @Override
        public void write(final int b) {
            reserve(length + 1)[length++] = (byte) b;
        }

        @Override
        public void write(final byte[] b, final int offset, final int count) {
            System.arraycopy(b, offset, reserve(length + count), length, count);
            length += count;
        }
    }

    /**
     * htsjdk 5.0.0's AsyncSAMFileWriter (an AbstractAsyncWriter) with records handed to the writer thread in batches,
     * so that the queue's lock and the writer thread's wake-up are paid once per batch rather than once per record.
     * Failures follow AbstractAsyncWriter's: a writer-thread failure is rethrown, as it is if an Error or a
     * RuntimeException and otherwise wrapped in a RuntimeException, by a later addAlignment or by close, and marks the
     * writer closed, so a close after that rethrow does nothing and leaves the underlying writer open.
     */
    private static final class AsyncWriter implements SAMFileWriter {
        private static final int BATCH_SIZE = 1024;
        private static final int QUEUED_BATCHES = 2;
        private static final AtomicInteger threadsCreated = new AtomicInteger();
        // Queued by close after the last records; the writer thread stops when it takes it.
        private static final SAMRecord[] END = new SAMRecord[0];

        private final SAMFileWriter underlyingWriter;
        private final BlockingQueue<SAMRecord[]> queue = new ArrayBlockingQueue<>(QUEUED_BATCHES);
        private final AtomicBoolean isClosed = new AtomicBoolean(false);
        private final AtomicReference<Throwable> ex = new AtomicReference<>(null);
        private final Thread writer;
        private SAMRecord[] batch = new SAMRecord[BATCH_SIZE];
        private int batchSize = 0;

        AsyncWriter(final SAMFileWriter out) {
            this.underlyingWriter = out;
            this.writer = new Thread(this::writeBatches, "SAMFileWriterThread-" + threadsCreated.getAndIncrement());
            this.writer.setDaemon(true);
            this.writer.start();
        }

        private void writeBatches() {
            try {
                while (true) {
                    final SAMRecord[] records;
                    try {
                        records = queue.take();
                    } catch (final InterruptedException ie) {
                        // As AbstractAsyncWriter's writer thread: keep writing.
                        continue;
                    }
                    if (records == END) {
                        return;
                    }
                    for (final SAMRecord record : records) {
                        underlyingWriter.addAlignment(record);
                    }
                }
            } catch (final Throwable t) {
                ex.compareAndSet(null, t);
                // Unblocks a put on the full queue, so that the caller sees the failure.
                queue.clear();
            }
        }

        @Override
        public void addAlignment(final SAMRecord alignment) {
            if (isClosed.get()) {
                throw new RuntimeIOException("Attempt to add record to closed writer.");
            }
            checkAndRethrow();
            batch[batchSize++] = alignment;
            if (batchSize == BATCH_SIZE) {
                final SAMRecord[] full = batch;
                batch = new SAMRecord[BATCH_SIZE];
                batchSize = 0;
                try {
                    queue.put(full);
                } catch (final InterruptedException ie) {
                    throw new RuntimeException("Interrupted queueing item for writing.", ie);
                }
                checkAndRethrow();
            }
        }

        @Override
        public void close() {
            checkAndRethrow();
            if (!isClosed.getAndSet(true)) {
                try {
                    // Two puts at most, so neither blocks for good on a writer thread that failed and cleared the queue.
                    if (batchSize > 0) {
                        queue.put(Arrays.copyOf(batch, batchSize));
                    }
                    queue.put(END);
                    writer.join();
                } catch (final InterruptedException ie) {
                    throw new RuntimeException("Interrupted waiting on writer thread.", ie);
                }
                batch = null;
                underlyingWriter.close();
                checkAndRethrow();
            }
        }

        private void checkAndRethrow() {
            if (ex.get() == null) {
                return;
            }
            final Throwable t = ex.getAndSet(null);
            if (t != null) {
                isClosed.set(true);
                if (t instanceof Error) {
                    throw (Error) t;
                }
                if (t instanceof RuntimeException) {
                    throw (RuntimeException) t;
                }
                throw new RuntimeException(t);
            }
        }

        @Override
        public void setProgressLogger(final ProgressLoggerInterface progress) {
            underlyingWriter.setProgressLogger(progress);
        }

        @Override
        public void setSortOrderChecking(final boolean check) {
            underlyingWriter.setSortOrderChecking(check);
        }

        @Override
        public SAMFileHeader getFileHeader() {
            return underlyingWriter.getFileHeader();
        }
    }
}
