package picard.sam.util;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.zip.DeflaterFactory;

import java.io.IOException;
import java.io.InterruptedIOException;
import java.io.OutputStream;
import java.util.ArrayDeque;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.CRC32;
import java.util.zip.Deflater;

/**
 * A BGZF output stream that deflates its blocks on several threads and writes them in order.
 *
 * Given the same writes, flushes, compression level and deflater factory, its bytes are those of htsjdk's
 * {@link htsjdk.samtools.util.BlockCompressedOutputStream}: a block is cut when the buffer is exactly full or on
 * {@link #flush()}, each block is deflated with the same calls on a deflater of its own, falling back to a level-0
 * JDK deflater when the output does not fit, and {@link #close()} writes the empty EOF block. It has no file pointer
 * and no index. One thread at a time may write to it; while it waits for the oldest block, that thread deflates
 * queued blocks itself.
 */
public class ParallelBlockCompressedOutputStream extends OutputStream {
    private static final int UNCOMPRESSED_BLOCK_SIZE = BlockCompressedStreamConstants.DEFAULT_UNCOMPRESSED_BLOCK_SIZE;
    private static final int HEADER_LENGTH = BlockCompressedStreamConstants.BLOCK_HEADER_LENGTH;
    private static final int FOOTER_LENGTH = BlockCompressedStreamConstants.BLOCK_FOOTER_LENGTH;
    // The size of BlockCompressedOutputStream's compressedBuffer, which decides when the level-0 fallback is taken.
    private static final int COMPRESSED_BUFFER_SIZE = BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE - HEADER_LENGTH;
    private static final AtomicInteger streamsCreated = new AtomicInteger();

    private final BinaryCodec codec;
    private final int compressionLevel;
    private final DeflaterFactory deflaterFactory;
    private final ThreadPoolExecutor executor;
    private final int maxBlocksInFlight;
    // The blocks submitted and not yet written, at most: maxBlocksInFlight unless set otherwise; set by any thread.
    private volatile int blocksInFlightLimit;
    private final ConcurrentLinkedQueue<Compressor> idleCompressors = new ConcurrentLinkedQueue<>();
    private final ArrayDeque<Future<Block>> inFlight = new ArrayDeque<>();
    private final ArrayDeque<Block> freeBlocks = new ArrayDeque<>();
    private Block current = new Block();
    private boolean closed = false;

    /**
     * @param os               where the compressed blocks go; closed by {@link #close()}
     * @param compressionLevel the deflate level, as given to BlockCompressedOutputStream
     * @param deflaterFactory  makes one deflater per worker thread
     * @param threads          the number of threads that deflate blocks
     */
    public ParallelBlockCompressedOutputStream(final OutputStream os, final int compressionLevel,
                                               final DeflaterFactory deflaterFactory, final int threads) {
        if (threads < 1) {
            throw new IllegalArgumentException("threads must be at least 1: " + threads);
        }
        this.codec = new BinaryCodec(os);
        this.compressionLevel = compressionLevel;
        this.deflaterFactory = deflaterFactory;
        this.maxBlocksInFlight = 4 * threads;
        this.blocksInFlightLimit = this.maxBlocksInFlight;

        final String namePrefix = "BGZFDeflater-" + streamsCreated.incrementAndGet() + "-";
        final AtomicInteger threadsCreated = new AtomicInteger();
        // Daemon threads that time out, so a writer abandoned after an error leaves no thread behind.
        this.executor = new ThreadPoolExecutor(threads, threads, 10, TimeUnit.SECONDS, new LinkedBlockingQueue<>(), r -> {
            final Thread thread = new Thread(r, namePrefix + threadsCreated.incrementAndGet());
            thread.setDaemon(true);
            return thread;
        });
        this.executor.allowCoreThreadTimeOut(true);
    }

    /**
     * Keeps at most this many blocks submitted and not yet written, from the next block on, in place of four per
     * thread; below 1, four per thread again. The output is the same at any limit.
     */
    public void setBlocksInFlightLimit(final int limit) {
        this.blocksInFlightLimit = limit < 1 ? maxBlocksInFlight : limit;
    }

    @Override
    public void write(final byte[] bytes) throws IOException {
        write(bytes, 0, bytes.length);
    }

    @Override
    public void write(final byte[] bytes, int startIndex, int numBytes) throws IOException {
        while (numBytes > 0) {
            final Block block = current;
            final int bytesToWrite = Math.min(UNCOMPRESSED_BLOCK_SIZE - block.uncompressedLength, numBytes);
            System.arraycopy(bytes, startIndex, block.uncompressed, block.uncompressedLength, bytesToWrite);
            block.uncompressedLength += bytesToWrite;
            startIndex += bytesToWrite;
            numBytes -= bytesToWrite;
            if (block.uncompressedLength == UNCOMPRESSED_BLOCK_SIZE) {
                submitCurrent();
            }
        }
    }

    @Override
    public void write(final int b) throws IOException {
        final Block block = current;
        block.uncompressed[block.uncompressedLength++] = (byte) b;
        if (block.uncompressedLength == UNCOMPRESSED_BLOCK_SIZE) {
            submitCurrent();
        }
    }

    /** Ends the current block, if it holds anything, and writes every block deflated or in progress. */
    @Override
    public void flush() throws IOException {
        if (current.uncompressedLength > 0) {
            submitCurrent();
        }
        while (!inFlight.isEmpty()) {
            writeOldest();
        }
        codec.getOutputStream().flush();
    }

    /** Flushes, writes the empty EOF block and closes the underlying stream. */
    @Override
    public void close() throws IOException {
        if (closed) {
            return;
        }
        closed = true;
        try {
            flush();
            codec.writeBytes(BlockCompressedStreamConstants.EMPTY_GZIP_BLOCK);
            codec.close();
        } finally {
            executor.shutdown();
            try {
                // A deflater may be ended only once no task can still be using it.
                if (executor.awaitTermination(1, TimeUnit.MINUTES)) {
                    Compressor compressor;
                    while ((compressor = idleCompressors.poll()) != null) {
                        compressor.end();
                    }
                }
            } catch (final InterruptedException e) {
                Thread.currentThread().interrupt();
            }
        }
    }

    private void submitCurrent() throws IOException {
        final Block block = current;
        inFlight.addLast(executor.submit(() -> deflate(block)));
        final int limit = blocksInFlightLimit;
        while (!inFlight.isEmpty() && (inFlight.size() >= limit || inFlight.peekFirst().isDone())) {
            writeOldest();
        }
        current = freeBlocks.isEmpty() ? new Block() : freeBlocks.pop();
    }

    private void writeOldest() throws IOException {
        final Future<Block> oldest = inFlight.peekFirst();
        // Deflate blocks still waiting for a worker on this thread rather than wait idle for the oldest one.
        while (!oldest.isDone() && deflateQueuedBlock()) {
        }
        final Block block = await(inFlight.removeFirst());
        codec.writeBytes(block.gzip, 0, block.gzipLength);
        block.uncompressedLength = 0;
        freeBlocks.push(block);
    }

    /**
     * Runs one deflate task that no worker has started, on the calling thread, and returns whether there was one.
     * Its future is in {@link #inFlight}, so the block is still written in order.
     */
    private boolean deflateQueuedBlock() {
        final Runnable task = executor.getQueue().poll();
        if (task == null) {
            return false;
        }
        task.run();
        return true;
    }

    private Block deflate(final Block block) {
        Compressor compressor = idleCompressors.poll();
        if (compressor == null) {
            compressor = new Compressor();
        }
        try {
            compressor.compress(block);
        } finally {
            idleCompressors.add(compressor);
        }
        return block;
    }

    private static Block await(final Future<Block> future) throws IOException {
        try {
            return future.get();
        } catch (final InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new InterruptedIOException("Interrupted while waiting for a BGZF block to be deflated");
        } catch (final ExecutionException e) {
            final Throwable cause = e.getCause();
            if (cause instanceof RuntimeException) {
                throw (RuntimeException) cause;
            }
            if (cause instanceof Error) {
                throw (Error) cause;
            }
            throw new IOException(cause);
        }
    }

    private static final class Block {
        final byte[] uncompressed = new byte[UNCOMPRESSED_BLOCK_SIZE];
        int uncompressedLength = 0;
        final byte[] gzip = new byte[HEADER_LENGTH + COMPRESSED_BUFFER_SIZE + FOOTER_LENGTH];
        int gzipLength = 0;
    }

    /** The deflaters, checksum and buffer of BlockCompressedOutputStream, one set per worker thread. */
    private final class Compressor {
        private final Deflater deflater = deflaterFactory.makeDeflater(compressionLevel, true);
        private final Deflater noCompressionDeflater = new Deflater(Deflater.NO_COMPRESSION, true);
        private final CRC32 crc32 = new CRC32();
        // GKL's IntelDeflater accepts only offset 0, so the block is deflated here and then copied into place.
        private final byte[] compressedBuffer = new byte[COMPRESSED_BUFFER_SIZE];

        /** BlockCompressedOutputStream.deflateBlock and writeGzipBlock, into block.gzip. */
        void compress(final Block block) {
            final int bytesToCompress = block.uncompressedLength;
            deflater.reset();
            deflater.setInput(block.uncompressed, 0, bytesToCompress);
            deflater.finish();
            int compressedSize = deflater.deflate(compressedBuffer, 0, compressedBuffer.length);
            if (!deflater.finished()) {
                noCompressionDeflater.reset();
                noCompressionDeflater.setInput(block.uncompressed, 0, bytesToCompress);
                noCompressionDeflater.finish();
                compressedSize = noCompressionDeflater.deflate(compressedBuffer, 0, compressedBuffer.length);
                if (!noCompressionDeflater.finished()) {
                    throw new IllegalStateException("unpossible");
                }
            }
            crc32.reset();
            crc32.update(block.uncompressed, 0, bytesToCompress);

            final byte[] gzip = block.gzip;
            final int totalBlockSize = compressedSize + HEADER_LENGTH + FOOTER_LENGTH;
            gzip[0] = BlockCompressedStreamConstants.GZIP_ID1;
            gzip[1] = (byte) BlockCompressedStreamConstants.GZIP_ID2;
            gzip[2] = BlockCompressedStreamConstants.GZIP_CM_DEFLATE;
            gzip[3] = (byte) BlockCompressedStreamConstants.GZIP_FLG;
            putInt(gzip, 4, 0);
            gzip[8] = (byte) BlockCompressedStreamConstants.GZIP_XFL;
            gzip[9] = (byte) BlockCompressedStreamConstants.GZIP_OS_UNKNOWN;
            putShort(gzip, 10, BlockCompressedStreamConstants.GZIP_XLEN);
            gzip[12] = BlockCompressedStreamConstants.BGZF_ID1;
            gzip[13] = BlockCompressedStreamConstants.BGZF_ID2;
            putShort(gzip, 14, BlockCompressedStreamConstants.BGZF_LEN);
            putShort(gzip, 16, (short) (totalBlockSize - 1));
            System.arraycopy(compressedBuffer, 0, gzip, HEADER_LENGTH, compressedSize);
            putInt(gzip, HEADER_LENGTH + compressedSize, (int) crc32.getValue());
            putInt(gzip, HEADER_LENGTH + compressedSize + 4, bytesToCompress);
            block.gzipLength = totalBlockSize;
        }

        void end() {
            deflater.end();
            noCompressionDeflater.end();
        }
    }

    private static void putShort(final byte[] bytes, final int offset, final short value) {
        bytes[offset] = (byte) value;
        bytes[offset + 1] = (byte) (value >>> 8);
    }

    private static void putInt(final byte[] bytes, final int offset, final int value) {
        bytes[offset] = (byte) value;
        bytes[offset + 1] = (byte) (value >>> 8);
        bytes[offset + 2] = (byte) (value >>> 16);
        bytes[offset + 3] = (byte) (value >>> 24);
    }
}
