package net.sf.samtools;

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
class AsyncSAMFileWriter implements SAMFileWriter {
    private static volatile int threadsCreated = 0; // Just used for thread naming.
    public static final int DEFAULT_QUEUE_SIZE = 2000;

    private final SAMFileWriter out;
    private final AtomicBoolean isClosed = new AtomicBoolean(false);
    private final BlockingQueue<SAMRecord> queue;
    private final Thread writer;
    private final AtomicReference<Throwable> ex = new AtomicReference<Throwable>(null);

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
        this.out = out;
        this.queue = new ArrayBlockingQueue<SAMRecord>(queueSize);
        this.writer = new Thread(new WriterRunnable(), "SAMFileWriterThread-" + threadsCreated++);
        this.writer.setDaemon(true);
        this.writer.start();
    }

    /**
     * Adds an alignment to the queue to be written.  Will re-throw any exception that was received when
     * writing prior record(s) to the underlying SAMFileWriter.
     * @param alignment
     */
    public void addAlignment(final SAMRecord alignment) {
        if (isClosed.get()) throw new RuntimeException("Attempt to add record to closed SAMFileWriter.");

        checkAndRethrow();
        try { this.queue.put(alignment); }
        catch (InterruptedException ie) { throw new RuntimeException("Interrupted queueing SAMRecord for writing.", ie); }
        checkAndRethrow();
    }

    /** Returns the SAMFileHeader from the underlying SAMFileWriter. */
    public SAMFileHeader getFileHeader() {
        return this.out.getFileHeader();
    }

    /**
     * Attempts to finishing draining the queue and then close the underlying SAMFileWriter.
     */
    public void close() {
        checkAndRethrow();

        if (this.isClosed.get()) {
            throw new RuntimeException("SAMFileWriter already closed.");
        }
        else {
            this.isClosed.set(true);

            try { this.writer.join(); }
            catch (InterruptedException ie) { throw new RuntimeException("Interrupted waiting on writer thread.", ie); }

            // Assert that the queue is empty
            if (!this.queue.isEmpty()) {
                throw new RuntimeException("Queue should be empty but is size: " + this.queue.size());
            }

            out.close();
            checkAndRethrow();
        }
    }

    /**
     * Checks to see if an exception has been raised in the writer thread and if so rethrows it as an Error
     * or RuntimeException as appropriate.
     */
    private final void checkAndRethrow() {
        if (this.ex.get() != null) {
            final Throwable t = this.ex.get();

            if (t instanceof Error) throw (Error) t;
            if (t instanceof RuntimeException) throw (RuntimeException) t;
            else throw new RuntimeException(t);
        }
    }

    /**
     * Small Runnable implementation that simply reads from the blocking queue and writes to the
     * output SAMFileWriter.
     */
    private class WriterRunnable implements Runnable {
        public void run() {
            try {
                while (!queue.isEmpty() || !isClosed.get()) {
                    try {
                        final SAMRecord rec = queue.poll(2, TimeUnit.SECONDS);
                        if (rec != null) out.addAlignment(rec);
                    }
                    catch (InterruptedException ie) {
                        /* Do Nothing */
                    }
                }
            }
            catch (Throwable t) {
                ex.compareAndSet(null, t);
            }
        }
    }
}
