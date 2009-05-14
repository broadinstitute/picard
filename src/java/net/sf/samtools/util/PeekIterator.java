/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools.util;

import java.util.Iterator;

/**
 * Wrapper around an iterator that enables non-destructive peeking at the next element that would
 * be returned by next()
 */
public class PeekIterator<T> implements Iterator<T> {
    Iterator<T> underlyingIterator;
    T peekedElement = null;

    public PeekIterator(final Iterator<T> underlyingIterator) {
        this.underlyingIterator = underlyingIterator;
    }

    /**
     * @return true if the iteration has more elements. (In other words, returns true if next would return an element 
     * rather than throwing an exception.)
     */
    public boolean hasNext() {
        return peekedElement != null || underlyingIterator.hasNext();  
    }

    /**
     * @return the next element in the iteration. Calling this method repeatedly until the hasNext() method returns
     * false will return each element in the underlying collection exactly once.
     */
    public T next() {
        if (peekedElement != null) {
            final T ret = peekedElement;
            peekedElement = null;
            return ret;
        }
        return underlyingIterator.next();
    }

    /**
     * @return the next element in the iteration, but without removing it, so the next call to next() or peek()
     * will return the same element as returned by the current call to peek().
     */
    public T peek() {
        if (peekedElement == null) {
            peekedElement = underlyingIterator.next();
        }
        return peekedElement;
    }

    /**
     * Unsupported
     */
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * @return the iterator wrapped by this object.
     */
    public Iterator<T> getUnderlyingIterator() {
        return underlyingIterator;
    }
}
