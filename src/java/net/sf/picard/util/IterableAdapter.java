package net.sf.picard.util;

import java.util.ConcurrentModificationException;
import java.util.Iterator;

/**
 * Provides an adapter to wrap an Iterator with an Iterable, allowing it to be run through a foreach loop. Will only
 * allow iterator() to be called a single time - this is intended to be called inline.
 *
 * @author jgentry@broadinstitute.org
 */
public class IterableAdapter<T> implements Iterable<T> {
    private boolean isIteratorCalled = false;
    private final Iterator<T> theIterator;

    public IterableAdapter(final Iterator<T> theIterator) {
        this.theIterator = theIterator;
    }

    @Override
    public Iterator<T> iterator() {
        if (isIteratorCalled) {
            throw new ConcurrentModificationException("iterator() can only be called once!");
        }

        isIteratorCalled = true;
        return theIterator;
    }
}
