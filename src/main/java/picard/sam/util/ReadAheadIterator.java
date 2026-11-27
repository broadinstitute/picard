package picard.sam.util;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/**
 * Reads an iterator ahead on a background thread, in batches, so that the work done inside the wrapped iterator's
 * hasNext and next (for a SamReader's iterator: decoding and validating each record) runs beside the consumer's.
 * Elements come out in the wrapped iterator's order. An exception the wrapped iterator throws is rethrown unchanged,
 * by hasNext or next, once the elements before it have been returned.
 *
 * Unlike htsjdk's AsyncBufferedIterator, close never interrupts the background thread. An interrupt that reaches a
 * reader with asynchronous BGZF reading while it waits for a block becomes a RuntimeException and is cleared, after
 * which AsyncBufferedIterator's thread goes on reading, blocks on its full queue, and its close never returns.
 */
public final class ReadAheadIterator<T> implements CloseableIterator<T> {
    private static final int QUEUED_BATCHES = 2;

    private static final class Batch {
        final Object[] elements;
        final int size;
        // Thrown once the elements are consumed.
        final Throwable failure;
        final boolean last;

        Batch(final Object[] elements, final int size, final Throwable failure, final boolean last) {
            this.elements = elements;
            this.size = size;
            this.failure = failure;
            this.last = last;
        }
    }

    private final CloseableIterator<T> wrapped;
    private final int batchSize;
    private final BlockingQueue<Batch> batches = new ArrayBlockingQueue<>(QUEUED_BATCHES);
    private final Thread thread;
    private volatile boolean stopping = false;

    private Batch current = new Batch(new Object[0], 0, null, false);
    private int position = 0;
    private boolean closed = false;

    public ReadAheadIterator(final CloseableIterator<T> wrapped, final int batchSize, final String threadName) {
        if (batchSize < 1) {
            throw new IllegalArgumentException("batchSize must be at least 1: " + batchSize);
        }
        this.wrapped = wrapped;
        this.batchSize = batchSize;
        this.thread = new Thread(this::readAhead, threadName);
        this.thread.setDaemon(true);
        this.thread.start();
    }

    private void readAhead() {
        boolean last = false;
        while (!last) {
            final Object[] elements = new Object[batchSize];
            int size = 0;
            Throwable failure = null;
            try {
                while (size < batchSize) {
                    if (stopping) {
                        return;
                    }
                    if (!wrapped.hasNext()) {
                        last = true;
                        break;
                    }
                    // Not elements[size++] = wrapped.next(): Java increments size before calling next, so a throw
                    // would leave a null element ahead of the exception.
                    final Object element = wrapped.next();
                    elements[size++] = element;
                }
            } catch (final Throwable t) {
                failure = t;
                last = true;
            }
            // close() empties the queue after setting stopping, so this put cannot block once close has begun.
            boolean interrupted = false;
            while (true) {
                try {
                    batches.put(new Batch(elements, size, failure, last));
                    break;
                } catch (final InterruptedException e) {
                    interrupted = true;
                }
            }
            if (interrupted) {
                Thread.currentThread().interrupt();
            }
        }
    }

    @Override
    public boolean hasNext() {
        if (closed) {
            throw new IllegalStateException("iterator has been closed");
        }
        while (position == current.size) {
            if (current.failure != null) {
                rethrow(current.failure);
            }
            if (current.last) {
                return false;
            }
            try {
                current = batches.take();
            } catch (final InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new RuntimeException("Interrupted waiting for the read-ahead thread", e);
            }
            position = 0;
        }
        return true;
    }

    @Override
    @SuppressWarnings("unchecked")
    public T next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        final T element = (T) current.elements[position];
        current.elements[position++] = null;
        return element;
    }

    private static void rethrow(final Throwable failure) {
        if (failure instanceof RuntimeException) {
            throw (RuntimeException) failure;
        }
        if (failure instanceof Error) {
            throw (Error) failure;
        }
        throw new RuntimeException(failure);
    }

    /** Stops the background thread, waits for it, and closes the wrapped iterator on the calling thread. */
    @Override
    public void close() {
        if (closed) {
            return;
        }
        closed = true;
        stopping = true;
        batches.clear();
        boolean interrupted = false;
        while (true) {
            try {
                thread.join();
                break;
            } catch (final InterruptedException e) {
                interrupted = true;
            }
        }
        if (interrupted) {
            Thread.currentThread().interrupt();
        }
        current = null;
        CloserUtil.close(wrapped);
    }
}
