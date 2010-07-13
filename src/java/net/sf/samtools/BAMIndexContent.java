/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import java.util.*;

/** Represents the contents of a bam index file for one reference.
 * A BAM index (.bai) file contains information for all references in the bam file.
 * This class describes the data present in the index file for one of these references;
 * including the bins, chunks, and linear index.
 */
class BAMIndexContent {
    /**
     * The reference sequence for the data currently loaded.
     */
    private final int mReferenceSequence;

    /**
     * A list of all bins in the above reference sequence.
     */
    private final List<Bin> mBins;

    /**
     * The linear index for the reference sequence above.
     */
    private final LinearIndex mLinearIndex;

    public BAMIndexContent(final int referenceSequence, final List<Bin> bins, final LinearIndex linearIndex) {
        this.mReferenceSequence = referenceSequence;
        this.mBins = bins;
        this.mLinearIndex = linearIndex;
    }

    public int getReferenceSequence() {
        return mReferenceSequence;
    }

    public boolean containsBin(final Bin bin) {
        return Collections.binarySearch(mBins,bin) >= 0;
    }

    public List<Bin> getBins() {
        return Collections.unmodifiableList(mBins);
    }

    public List<Chunk> getChunksForBin(final Bin bin) {
        return Collections.unmodifiableList(bin.getChunkList());
    }

    public List<Chunk> getAllChunks() {
        List<Chunk> allChunks = new ArrayList<Chunk>();
        for(Bin bin: getBins()){
            if (getChunksForBin(bin) != null){
                allChunks.addAll(getChunksForBin(bin));
            }
        }
        return Collections.unmodifiableList(allChunks);
    }

    public LinearIndex getLinearIndex() {
        return mLinearIndex;
    }
}
