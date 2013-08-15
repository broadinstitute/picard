package org.broad.tribble.readers;

import java.util.Iterator;

/**
 * A very simple descriptor for line-iterables.
 * @author mccowan
 */
public interface LineIterator extends Iterator<String> {
    /** Peeks at the next line, without expending any elements in the underlying iterator. */
    public String peek();
}
