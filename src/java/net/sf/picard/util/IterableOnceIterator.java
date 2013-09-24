package net.sf.picard.util;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;

/**
 * Abstract implementation of an iterator that also implements Iterable (to return itself)
 * so that it can be used if for() loops.  Only supports calling iterator() once since new
 * iterators are not manufactured but the same object returned.
 *
 * @author Tim Fennell
 */
public abstract class IterableOnceIterator<T> implements Iterable<T>, Iterator<T>, Closeable {
    private boolean iterated = false;

    /**
     * On the first call returns this object which is also an iterator.  On subsequent calls throws
     * an exception since new iterators cannot be generated.
     */
    @Override
    public Iterator<T> iterator() {
        if (iterated) {
            throw new IllegalStateException("May not call iterator() more than once on IterableOnceIterator.");
        }
        else {
            iterated = true;
            return this;
        }
    }

    /** Operation not supported. */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported");
    }

    /** Does nothing, intended to be overridden when needed. */
    @Override public void close() throws IOException {
        // Default do nothing implementation
    }
}
