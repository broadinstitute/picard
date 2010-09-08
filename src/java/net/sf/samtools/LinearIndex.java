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

/**
 * The linear index associated with a given reference in a BAM index.
 *
 * @author mhanna
 * @version 0.1
 */
public class LinearIndex {

    public static final int MAX_LINEAR_INDEX_SIZE = AbstractBAMFileIndex.MAX_LINEAR_INDEX_SIZE;

    public static final int BAM_LIDX_SHIFT = 14;

    /**
     * The reference sequence number for this linear index.
     */
    private final int mReferenceSequence;

    /**
     * Dictates the first stored element of the index.
     */
    private final int mIndexStart;

    /**
     * The linear index entries within this bin.
     */
    private final long[] mIndexEntries;

    public LinearIndex(final int referenceSequence, final int indexStart, final long[] indexEntries) {
        this.mReferenceSequence = referenceSequence;
        this.mIndexStart = indexStart;
        this.mIndexEntries = indexEntries;
    }

    public int getReferenceSequence() {
        return mReferenceSequence;
    }

    public int size() {
        return mIndexEntries.length;
    }

    public long get(final int index) {
        return mIndexEntries[index-mIndexStart];
    }

    public static int convertToLinearIndexOffset(final int contigPos) {
        final int indexPos = (contigPos <= 0) ? 0 : contigPos-1;
        return indexPos >> BAM_LIDX_SHIFT;
    }

    /**
     * Gets the minimum offset of any alignment start appearing in this index, according to the linear index. 
     * @param startPos Starting position for this query.
     * @return The minimum offset, in chunk format, of any read appearing in this position.
     */
    public long getMinimumOffset(final int startPos) {
        final int start = (startPos <= 0) ? 0 : startPos-1;
        final int regionLinearBin = start >> BAM_LIDX_SHIFT;
        // System.out.println("# regionLinearBin: " + regionLinearBin);
        long minimumOffset = 0;
        if (regionLinearBin-mIndexStart < mIndexEntries.length)
            minimumOffset = mIndexEntries[regionLinearBin-mIndexStart];
        return minimumOffset;
    }

    /**
     * Direct access to the array.  Be careful!
     * @return The elements of the linear index.
     */
    protected long[] getIndexEntries() {
        return mIndexEntries;
    }

     protected int getIndexStart() {
        return mIndexStart;
    }
}
