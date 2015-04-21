package picard.util;

import com.google.common.base.Optional;

/**
 * Describes 
 * @author mccowan
 */
public interface AtomicIterator<T> {
    /** Produces the next element from the iterator, if there is one; otherwise, produces {@link com.google.common.base.Optional.Absent} */
    Optional<T> next();
}
