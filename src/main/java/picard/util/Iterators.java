package picard.util;


import java.util.Iterator;
import java.util.Optional;

/**
 * @author mccowan
 */
public class Iterators {
    public static <E> AtomicIterator<E> atomicIteratorOf(final Iterable<E> backingIterable) {
        return atomicIteratorOf(backingIterable.iterator());
    }

    public static <E> AtomicIterator<E> atomicIteratorOf(final Iterator<E> backingIterator) {
        final Object monitor = new Object();
        return () -> {
            synchronized (monitor) {
                return backingIterator.hasNext() ? Optional.ofNullable(backingIterator.next()) : Optional.empty();
            }
        };
    }
}
