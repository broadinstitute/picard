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
 * This interface is used by iterators that use releasable resources during iteration.
 * 
 * The consumer of a CloseableIterator should ensure that the close() method is always called,
 * for example by putting such a call in a finally block.  Two conventions should be followed
 * by all implementors of CloseableIterator:
 * 1) The close() method should be idempotent.  Calling close() twice should have no effect.
 * 2) When hasNext() returns false, the iterator implementation should automatically close itself.
 *    The latter makes it somewhat safer for consumers to use the for loop syntax for iteration:
 *    for (Type obj : getCloseableIterator()) { ... }
 * 
 * We do not inherit from java.io.Closeable because IOExceptions are a pain to deal with.
 */
public interface CloseableIterator<T>
    extends Iterator<T> {

    public void close();
}
