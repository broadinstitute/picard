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
package net.sf.picard.util;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;

import java.util.Iterator;

/**
 * Generic Closable Iterator that allows you to peek at the next value before calling next
 */
public class PeekableIterator<Object> implements CloseableIterator<Object> {
    private Iterator<Object> iterator;
    private Object nextObject;

    /** Constructs a new iterator that wraps the supplied iterator. */
    public PeekableIterator(Iterator<Object> iterator) {
        this.iterator = iterator;
        advance();
    }

    /** Closes the underlying iterator. */
    public void close() {
        CloserUtil.close(iterator);
    }

    /** True if there are more items, in which case both next() and peek() will return a value. */
    public boolean hasNext() {
        return this.nextObject != null;
    }

    /** Returns the next object and advances the iterator. */
    public Object next() {
        Object retval = this.nextObject;
        advance();
        return retval;
    }

    /**
     * Returns the next object but does not advance the iterator. Subsequent calls to peek()
     * and next() will return the same object.
     */
    public Object peek(){
        return this.nextObject;
    }

    private void advance(){
        if (this.iterator.hasNext()) {
            this.nextObject = iterator.next();
        }
        else {
            this.nextObject = null;
        }
    }

    /** Unsupported Operation. */
    public void remove() {
        throw new UnsupportedOperationException("Not supported: remove");
    }
}
