package net.sf.samtools.util;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Abstract class that is designed to be extended and specialized to provide an asynchronous
 * wrapper around any kind of Writer class that takes an object and writes it out somehow.
 *
 * @author Tim Fennell
 */
public abstract class AbstractAsyncWriter<T> {
    private static volatile int threadsCreated = 0; // Just used for thread naming.
    public static final int DEFAULT_QUEUE_SIZE = 2000;

    private final AtomicBoolean isClosed = new AtomicBoolean(false);
    private final BlockingQueue<T> queue;
    private final Thread writer;
    private final WriterRunnable writerRunnable;
    private final AtomicReference<Throwable> ex = new AtomicReference<Throwable>(null);

    /** Returns the prefix to use when naming threads. */
    protected abstract String getThreadNamePrefix();

    protected abstract void synchronouslyWrite(final T item);

    protected abstract void synchronouslyClose();

    /**
     * Creates an AbstractAsyncWriter that will use the provided WriterRunnable to consume from the
     * internal queue and write records into the synchronous writer.
     */
    protected AbstractAsyncWriter(final int queueSize) {
        this.queue = new ArrayBlockingQueue<T>(queueSize);
        this.writerRunnable = new WriterRunnable();
        this.writer = new Thread(writerRunnable, getThreadNamePrefix() + threadsCreated++);
        this.writer.setDaemon(true);
        this.writer.start();
    }

    /**
     * Public method for sub-classes or ultimately consumers to put an item into the queue
     * to be written out.
     */
    public void write(final T item) {
        if (this.isClosed.get()) throw new RuntimeException("Attempt to add record to closed writer.");

        checkAndRethrow();
        try { this.queue.put(item); }
        catch (InterruptedException ie) { throw new RuntimeException("Interrupted queueing item for writing.", ie); }
        checkAndRethrow();
    }

    /**
     * Attempts to finishing draining the queue and then calls synchronoslyClose() to allow implementation
     * to do any one time clean up.
     */
    public void close() {
        checkAndRethrow();

        if (this.isClosed.get()) {
            throw new RuntimeException("AbstractAsyncWriter already closed.");
        }
        else {
            this.isClosed.set(true);

            try { this.writer.join(); }
            catch (InterruptedException ie) { throw new RuntimeException("Interrupted waiting on writer thread.", ie); }

            // Assert that the queue is empty
            if (!this.queue.isEmpty()) {
                throw new RuntimeException("Queue should be empty but is size: " + this.queue.size());
            }

            synchronouslyClose();
            checkAndRethrow();
        }
    }

    /**
     * Checks to see if an exception has been raised in the writer thread and if so rethrows it as an Error
     * or RuntimeException as appropriate.
     */
    private final void checkAndRethrow() {
        final Throwable t = this.ex.get();
        if (t != null) {
            if (t instanceof Error) throw (Error) t;
            if (t instanceof RuntimeException) throw (RuntimeException) t;
            else throw new RuntimeException(t);
        }
    }

    /**
     * Small Runnable implementation that simply reads from the blocking queue and writes to the
     * synchronous writer.
     */
    private class WriterRunnable implements Runnable {
        public void run() {
            try {
                while (!queue.isEmpty() || !isClosed.get()) {
                    try {
                        final T item = queue.poll(2, TimeUnit.SECONDS);
                        if (item != null) synchronouslyWrite(item);
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
