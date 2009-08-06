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
package net.sf.picard.sam;

/** Little struct-like class to hold read pair (and fragment) end data for MarkDuplicates. */
class ReadEnds {
    public static final int SIZE_OF = (1*1) + (2*1) + (4*4) + (8*2) + 8 + // last 8 == reference overhead
                                                                 13; // This is determined experimentally with JProfiler
    public static final byte F=0, R=1, FF=2, FR=3, RR=4, RF=5;

    short libraryId;
    short score = 0;
    byte orientation;
    int read1Sequence     = -1;
    int read1Coordinate   = -1;
    long read1IndexInFile = -1;
    int read2Sequence     = -1;
    int read2Coordinate   = -1;
    long read2IndexInFile = -1;

    boolean isPaired() { return this.read2Sequence != -1; }
}
