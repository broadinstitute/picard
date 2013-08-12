package net.sf.picard.util;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;

/**
 * Abstract implementation of an iterator that also implements Iterable (to return itself)
 * so that it can be used if for() loops.
 *
 * @author Tim Fennell
 */
public abstract class IterableIterator<T> implements Iterable<T>, Iterator<T>, Closeable {
    @Override
    public Iterator<T> iterator() { return this; }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("remove() not supported");
    }

    @Override public void close() throws IOException {
        // Default do nothing implementation
    }
}
