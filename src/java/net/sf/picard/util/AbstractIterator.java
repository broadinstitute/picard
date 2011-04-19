/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

package net.sf.picard.util;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Base class of implementing iterators. All you have to do is implement advance which gets
 * the next element.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public abstract class AbstractIterator<E> implements Iterator<E> {
    private E next;
    private boolean iterating = false;

    @Override
    public boolean hasNext() {
        // If this is the start of iteration, queue up the first item
        if(!iterating) {
            next = advance();
            iterating = true;
        }
        return next != null;
    }

    @Override
    public E next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        
        E ret = next;
        next = advance();
        return ret;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("Remove() not supported.");
    }

    /**
     * @return the next element or null if the iterator is at the end
     */
    protected abstract E advance();

    /**
     * @return true after the first time hasNext() or next() have been called
     */
    protected boolean isIterating() {
        return iterating;
    }
}
