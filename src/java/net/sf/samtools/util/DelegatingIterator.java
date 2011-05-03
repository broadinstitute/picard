package net.sf.samtools.util;

import net.sf.samtools.util.CloseableIterator;

import java.util.Iterator;

/**
 * Simple iterator class that delegates all method calls to an underlying iterator. Useful
 * for in-line subclassing to add behaviour to one or more methods.
 *
 * @author Tim Fennell
 */
public class DelegatingIterator<T> implements CloseableIterator<T> {
    private final Iterator<T> iterator;

    public DelegatingIterator(final Iterator<T> iterator) {
        this.iterator = iterator;
    }

    public void close() {
        if (iterator instanceof CloseableIterator) {
            ((CloseableIterator) this.iterator).close();
        }
    }

    public boolean hasNext() {
        return this.iterator.hasNext();
    }

    public T next() {
        return this.iterator.next();
    }

    public void remove() {
        this.iterator.remove();
    }
}
