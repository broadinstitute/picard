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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

/**
 * While structurally identical to CompositeIndex, this class is maintained as it makes code more readable when the two are used together (see QSeqParser)
 * @author jburke@broadinstitute.org
 */
class Range {
    public final int start;
    public final int end;
    public final int length;
    public Range(final int start, final int end) {
        if(end < start) {
            throw new PicardException("Nonsensical Range(" + start + ", " + end + ")");
        }

        this.start = start;
        this.end   = end;
        this.length = end - start + 1;
    }

    @Override
    public boolean equals(final Object object) {
        if(object == null || !(object instanceof Range)) {
            return false;
        }

        final Range that = (Range) object;
        return that.start == this.start && that.end == this.end;
    }

    @Override
    public int hashCode() {
        return (int)Math.pow(start, end);
    }

    @Override
    public String toString() {
        return "Range(" + start + ", " + end + ", " + length +")";
    }
}