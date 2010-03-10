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
package net.sf.samtools;

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.PeekIterator;

/**
 * Wrapper around SAMRecord iterator that skips over non-primary elements.
 * This iterator conflates a filtering iterator and a peekable iterator.  It would be cleaner to
 * handle those concerns separately.
 */
public class NotPrimarySkippingIterator {
    private final PeekIterator<SAMRecord> it;

    public NotPrimarySkippingIterator(final CloseableIterator<SAMRecord> underlyingIt) {
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
        while (it.hasNext() && it.peek().getNotPrimaryAlignmentFlag()) {
            it.next();
        }
    }
}
