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

/**
 * Represents one set of cycles in an ReadStructure (e.g. if the ReadStructure is 36TB836T then
 * 36T, 8B, and 36T are invidually represented internally as a ReadDescriptor).
 */
public class ReadDescriptor {
    public final int length;
    public final ReadType type;

    public ReadDescriptor(final int length, final ReadType type) {
        this.length = length;
        this.type = type;
    }

    @Override
    public String toString() {
        return this.length + this.type.name();
    }

    public boolean equals(final Object other) {
        if(this == other) return true;
        if(other.getClass() != this.getClass()) return false;

        final ReadDescriptor that = (ReadDescriptor) other;
        return this.length == that.length && this.type == that.type;
    }

    @Override
    public int hashCode() {
        return 31 * this.type.ordinal() + 379 * length;
    }
}
