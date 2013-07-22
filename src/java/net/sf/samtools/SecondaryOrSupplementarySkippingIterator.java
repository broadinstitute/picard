package net.sf.samtools;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.PeekIterator;

/**
 * Wrapper around SAMRecord iterator that skips over secondary and supplementary elements.
 * This iterator conflates a filtering iterator and a peekable iterator.  It would be cleaner to
 * handle those concerns separately. This class should be viewed as a replacement for NotPrimarySkippingIterator,
 * in that we did not want to change the functionality of NPSI to no longer match its name
 */
public class SecondaryOrSupplementarySkippingIterator {
    private final PeekIterator<SAMRecord> it;

    public SecondaryOrSupplementarySkippingIterator(final CloseableIterator<SAMRecord> underlyingIt) {
        it = new PeekIterator<SAMRecord>(underlyingIt);
        skipAnyNotprimary();
    }

    public boolean hasCurrent() {
        return it.hasNext();
    }

    public SAMRecord getCurrent() {
        assert(hasCurrent());
        return it.peek();
    }

    public boolean advance() {
        it.next();
        skipAnyNotprimary();
        return hasCurrent();
    }

    private void skipAnyNotprimary() {
        while (it.hasNext() && it.peek().isSecondaryOrSupplementary()) {
            it.next();
        }
    }

}
