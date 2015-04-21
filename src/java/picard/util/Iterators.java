package picard.util;

import com.google.common.base.Optional;

import java.util.Iterator;

/**
 * @author mccowan
 */
public class Iterators {
    public static <E> AtomicIterator<E> atomicIteratorOf(final Iterable<E> backingIterable) {
        return atomicIteratorOf(backingIterable.iterator());
    }

    public static <E> AtomicIterator<E> atomicIteratorOf(final Iterator<E> backingIterator) {
        final Object monitor = new Object();
        return new AtomicIterator<E>() {
            @Override
            public Optional<E> next() {
                synchronized (monitor) {
                    return backingIterator.hasNext() ? Optional.fromNullable(backingIterator.next()) : Optional.<E>absent();
                }
            }
        };
    }
}
