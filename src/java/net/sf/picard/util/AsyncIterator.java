/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.picard.util;

import net.sf.samtools.util.CloseableIterator;

import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicReference;

/**
 * Wrapper around a CloseableIterator that reads in a separate thread, for cases in which that might be
 * efficient.
 */
public class AsyncIterator<T> implements CloseableIterator<T> {
    private static volatile int threadsCreated = 0; // Just used for thread naming.
    public static final int DEFAULT_QUEUE_SIZE = 2000;

    private final AtomicBoolean isClosed = new AtomicBoolean(false);
    private final BlockingQueue<T> queue;
    private final Thread reader;
    private final ReaderRunnable readerRunnable;
    private final AtomicReference<Throwable> ex = new AtomicReference<Throwable>(null);
    private T theNext = null;
    private final CloseableIterator<T> underlyingIterator;


    public AsyncIterator(final CloseableIterator<T> underlyingIterator,
                            final int queueSize,
                            final String threadNamePrefix) {
        this.underlyingIterator = underlyingIterator;
        this.queue = new ArrayBlockingQueue<T>(queueSize);
        this.readerRunnable = new ReaderRunnable();
        this.reader = new Thread(readerRunnable, threadNamePrefix + threadsCreated++);
        this.reader.setDaemon(true);
        this.reader.start();
        getNext();
    }


    /**
     * Set theNext to the next item to be returned, or null if there are no more items.
     */
    private void getNext() {
        assertOpen();

        checkAndRethrow();
        try {
            theNext = null;
            while (!this.queue.isEmpty() || !this.readerRunnable.isDone()) {
                theNext = this.queue.poll(5, TimeUnit.SECONDS);
                checkAndRethrow();
                if (theNext != null) break;
            }
        } catch (InterruptedException ie) { throw new RuntimeException("Interrupted queueing item for writing.", ie); }
        checkAndRethrow();
    }

    public boolean hasNext() {
        assertOpen();
        return theNext != null;
    }

    public T next() {
        assertOpen();
        if (!hasNext()) throw new NoSuchElementException();
        final T ret = theNext;
        getNext();
        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Stops the thread and then calls synchronouslyClose() to allow implementation to do any one time clean up.
     */
    public void close() {
        checkAndRethrow();

        assertOpen();
        this.isClosed.set(true);

        try { this.reader.join(); }
        catch (InterruptedException ie) { throw new RuntimeException("Interrupted waiting on reader thread.", ie); }

        underlyingIterator.close();
        checkAndRethrow();
        this.queue.clear();
    }

    private void assertOpen() {
        if (this.isClosed.get()) {
            throw new RuntimeException("AsyncIterator already closed.");
        }
    }

    /**
     * Checks to see if an exception has been raised in the reader thread and if so rethrows it as an Error
     * or RuntimeException as appropriate.
     */
    private void checkAndRethrow() {
        final Throwable t = this.ex.get();
        if (t != null) {
            if (t instanceof Error) throw (Error) t;
            if (t instanceof RuntimeException) throw (RuntimeException) t;
            else throw new RuntimeException(t);
        }
    }

    /**
     * Small Runnable implementation that simply reads from underlying iterator and stores on the blocking queue.
     */
    private class ReaderRunnable implements Runnable {
        private final AtomicBoolean readerDone = new AtomicBoolean(false);

        public boolean isDone() { return readerDone.get(); }

        public void run() {
            try {
                boolean isEof = false;
                while (!isClosed.get() && !isEof) {
                    try {
                        if (!underlyingIterator.hasNext()) {
                            isEof = true;
                        } else {
                            final T item = underlyingIterator.next();
                            // Keep trying to put item on the queue unless close() has been called.
                            while (!isClosed.get() && !queue.offer(item, 2, TimeUnit.SECONDS)) {
                            }
                        }
                    }
                    catch (InterruptedException ie) {
                        /* Do Nothing */
                    }
                }
            }
            catch (Throwable t) {
                ex.compareAndSet(null, t);
            } finally {
                readerDone.set(true);
            }
        }
    }
}
