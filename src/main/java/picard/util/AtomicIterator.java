package picard.util;

import java.util.Optional;

/**
 * Describes 
 * @author mccowan
 */
public interface AtomicIterator<T> {
    /** Produces the next element from the iterator, if there is one; otherwise, produces {@link Optional#EMPTY} */
    Optional<T> next();
}
