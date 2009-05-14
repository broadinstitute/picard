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

package net.sf.picard.reference;

/**
 * Wrapper around a reference sequence that has been read from a reference file.
 *
 * @author Tim Fennell
 */
public class ReferenceSequence {
    private final String name;
    private final byte[] bases;
    private final int contigIndex;
    private final int length;

    /**
     * creates a fully formed ReferenceSequence
     *
     * @param name the name of the sequence from the source file
     * @param index the zero based index of this contig in the source file
     * @param bases the bases themselves stored as one-byte characters
     */
    public ReferenceSequence(String name, int index, byte[] bases) {
        this.name = name;
        this.contigIndex = index;
        this.bases = bases;
        this.length = bases.length;
    }

    /** Gets the set of names given to this sequence in the source file. */
    public String getName() { return name; }

    /**
     * Gets the array of bases that define this sequence. The bases can include any
     * letter and possibly include masking information in the form of lower case
     * letters.  This array is mutable (obviously!) and it NOT a clone of the array
     * held interally.  Do not modify it!!!
     */
    public byte[] getBases() { return bases; }

    /** Gets the 0-based index of this contig in the source file from which it came. */
    public int getContigIndex() { return contigIndex; }

    /** Gets the length of this reference sequence in bases. */
    public int length() { return length; }
    
    public String toString() {
        return "ReferenceSequence " + getName();
    }
}
