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
    private List<Bin> mBins;

    /**
     * Alternative representation of the bins;
     * Avoids having to copy the array to the mBins list when constructing bam index
     */
    private final Bin[] mBinArray;
    private final int mNumberOfBins;

    /**
     * Chunks containing metaData for the reference, e.g. number of aligned and unaligned records
     */
    private final List<Chunk> metaDataChunks;

    /**
     * The linear index for the reference sequence above.
     */
    private final LinearIndex mLinearIndex;

    /**
     * @param referenceSequence   Content corresponds to this reference.
     * @param bins                List of bins represented by this content
     * @param metaDataChunks      Chunks representing metaData
     * @param linearIndex         Additional index used to optimize queries
     */
    public BAMIndexContent(final int referenceSequence, final List<Bin> bins, final List<Chunk> metaDataChunks, final LinearIndex linearIndex) {
        this.mReferenceSequence = referenceSequence;
        this.mBins = bins;
        this.mBinArray = null;
        this.mNumberOfBins = 0;
        this.metaDataChunks = metaDataChunks;
        this.mLinearIndex = linearIndex;
    }

    /** Alternate constructor used when building an index.
     * Avoids copying bin array to bin list.
     * @param referenceSequence   Content corresponds to this reference.
     * @param bins                Array of bins represented by this content
     * @param numberOfBins        Maximum index in the bins array that has data
     * @param metaDataChunks      Chunks representing metaData
     * @param linearIndex         Additional index used to optimize queries
     */
    BAMIndexContent(final int referenceSequence, final Bin[] bins, final int numberOfBins, final List<Chunk> metaDataChunks, final LinearIndex linearIndex) {
        this.mReferenceSequence = referenceSequence;
        this.mBins = null;
        this.mBinArray = bins;
        this.mNumberOfBins = numberOfBins;
        this.metaDataChunks = metaDataChunks;
        this.mLinearIndex = linearIndex;
    }

    /** Reference for this Content */
    public int getReferenceSequence() {
        return mReferenceSequence;
    }

    /**
     * Does this content have anything in this bin?
    */
    public boolean containsBin(final Bin bin) {
        return Collections.binarySearch(mBins,bin) >= 0;
    }

    /**
    * @return list of bins represented by this content
    */
    public List<Bin> getBins() {
        if (mBins == null && mBinArray != null && mNumberOfBins > 0) {
            // copy the mBinArray to an unmodifiable list
            mBins = new ArrayList<Bin>();
            for (Bin bin : mBinArray) {
                if (bin != null) mBins.add(bin);
            }
        }
        return Collections.unmodifiableList(mBins);
    }

    /**
     * @return array of bins specified in the alternate constructor, if any; otherwise null
     */
    Bin[] getOriginalBins(){
        return mBinArray;
    }

    /**
     * @return  the number of bins represented by this content
     */
    int getNumberOfBins(){
        if (mBinArray != null)
            return mNumberOfBins;
        else
            return mBins.size();
    }

    /**
     * @return  the chunks associated with the specified bin
     */
    public List<Chunk> getChunksForBin(final Bin bin) {
        return Collections.unmodifiableList(bin.getChunkList());
    }

    /**
     * @return  the meta data chunks for this content
     */
    public List<Chunk> getMetaDataChunks() {
        return Collections.unmodifiableList(metaDataChunks);
    }

    /**
     * @return  all chunks associated with all bins in this content
     */
    public List<Chunk> getAllChunks() {
        List<Chunk> allChunks = new ArrayList<Chunk>();
        for(Bin b: mBins)
        if (b.getChunkList() != null){
            allChunks.addAll(b.getChunkList());
        }
        return Collections.unmodifiableList(allChunks);
    }

    /**
     * @return  the linear index represented by this content
     */
    public LinearIndex getLinearIndex() {
        return mLinearIndex;
    }
}
