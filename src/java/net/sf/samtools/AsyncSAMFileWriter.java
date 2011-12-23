package net.sf.samtools;

import net.sf.samtools.util.AbstractAsyncWriter;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

/**
 * SAMFileWriter that can be wrapped around an underlying SAMFileWriter to provide asynchronous output. Records
 * added are placed into a queue, the queue is then drained into the underlying SAMFileWriter by a thread owned
 * by the instance.
 *
 * Exceptions experienced by the writer thread will be emitted back to the caller in subsequent calls to either
 * addAlignment() or close().
 *
 * @author Tim Fennell
 */
class AsyncSAMFileWriter extends AbstractAsyncWriter<SAMRecord> implements SAMFileWriter {
    private final SAMFileWriter underlyingWriter;

    /**
     * Creates a new AsyncSAMFileWriter wrapping the provided SAMFileWriter.
     */
    public AsyncSAMFileWriter(final SAMFileWriter out) {
        this(out, DEFAULT_QUEUE_SIZE);
    }

    /**
     * Creates an AsyncSAMFileWriter wrapping the provided SAMFileWriter and using the specified
     * queue size for buffer SAMRecords.
     */
    public AsyncSAMFileWriter(final SAMFileWriter out, final int queueSize) {
        super(queueSize);
        this.underlyingWriter = out;
    }

    @Override protected void synchronouslyWrite(final SAMRecord item) { this.underlyingWriter.addAlignment(item); }

    @Override protected void synchronouslyClose() { this.underlyingWriter.close();  }

    @Override protected final String getThreadNamePrefix() { return "SAMFileWriterThread-"; }

    /**
     * Adds an alignment to the queue to be written.  Will re-throw any exception that was received when
     * writing prior record(s) to the underlying SAMFileWriter.
     */
    public void addAlignment(final SAMRecord alignment) {
        write(alignment);
    }

    /** Returns the SAMFileHeader from the underlying SAMFileWriter. */
    public SAMFileHeader getFileHeader() {
        return this.underlyingWriter.getFileHeader();
    }
}
