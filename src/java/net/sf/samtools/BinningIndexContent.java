/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

/**
 * In-memory representation of the binning index for a single reference.  BAM and Tabix are both binning indices
 * with slightly different disk formats but identical in-memory representations.
 */
public class BinningIndexContent {
    /**
     * The reference sequence for the data currently loaded.
     */
    private final int mReferenceSequence;

    /**
     * A list of all bins in the above reference sequence.
     */
    private final BinList mBinList;

    /**
     * The linear index for the reference sequence above.
     */
    private final LinearIndex mLinearIndex;


    /**
     * @param referenceSequence Content corresponds to this reference.
     * @param binList           Array of bins represented by this content, possibly sparse
     * @param linearIndex       Additional index used to optimize queries
     */
    public BinningIndexContent(final int referenceSequence, final BinList binList, final LinearIndex linearIndex) {
        this.mReferenceSequence = referenceSequence;
        this.mBinList = binList;
        this.mLinearIndex = linearIndex;
    }

    /**
     * Reference for this Content
     */
    public int getReferenceSequence() {
        return mReferenceSequence;
    }

    /**
     * Does this content have anything in this bin?
     */
    public boolean containsBin(final Bin bin) {
        return mBinList.getBin(bin.getBinNumber()) != null;
    }

    /**
     * @return iterable list of bins represented by this content
     */
    public BinList getBins() {
        return mBinList;
    }

    /**
     * @return the number of non-null bins represented by this content
     */
    int getNumberOfNonNullBins() {
        return mBinList.getNumberOfNonNullBins();
    }

    /**
     * @return all chunks associated with all bins in this content
     */
    public List<Chunk> getAllChunks() {
        final List<Chunk> allChunks = new ArrayList<Chunk>();
        for (final Bin b : mBinList)
            if (b.getChunkList() != null) {
                allChunks.addAll(b.getChunkList());
            }
        return Collections.unmodifiableList(allChunks);
    }

    /**
     * @return the linear index represented by this content
     */
    public LinearIndex getLinearIndex() {
        return mLinearIndex;
    }


    /**
     *
     * @param startPos 1-based, inclusive
     * @param endPos 1-based, inclusive
     * @return List of Chunks overlapping the given region.  May return null if there are none.
     */
    public List<Chunk> getChunksOverlapping(final int startPos, final int endPos) {
        final BitSet overlappingBins = GenomicIndexUtil.regionToBins(startPos,endPos);
        if (overlappingBins == null) return null;

        // System.out.println("# Sequence target TID: " + referenceIndex);
        final List<Bin> bins = new ArrayList<Bin>();
        for(final Bin bin: this.getBins()) {
            if (overlappingBins.get(bin.getBinNumber()))
                bins.add(bin);
        }

        if (bins.isEmpty()) {
            return null;
        }

        final List<Chunk> chunkList = new ArrayList<Chunk>();
        for(final Bin bin: bins) {
            for(final Chunk chunk: bin.getChunkList())
                chunkList.add(chunk.clone());
        }

        if (chunkList.isEmpty()) {
            return null;
        }

        return Chunk.optimizeChunkList(chunkList,this.getLinearIndex().getMinimumOffset(startPos));
    }
    /**
     * This class is used to encapsulate the list of Bins store in the BAMIndexContent
     * While it is currently represented as an array, we may decide to change it to an ArrayList or other structure
     */
    public static class BinList implements Iterable<Bin> {

        private final Bin[] mBinArray;
        public final int numberOfNonNullBins;
        public final int maxBinNumber;  // invariant: maxBinNumber = mBinArray.length -1 since array is 0 based

        /**
         * @param binArray            a sparse array representation of the bins. The index into the array is the bin number.
         * @param numberOfNonNullBins
         */
        public BinList(final Bin[] binArray, final int numberOfNonNullBins) {
            this.mBinArray = binArray;
            this.numberOfNonNullBins = numberOfNonNullBins;
            this.maxBinNumber = mBinArray.length - 1;
        }

        Bin getBin(final int binNumber) {
            if (binNumber > maxBinNumber) return null;
            return mBinArray[binNumber];
        }

        int getNumberOfNonNullBins() {
            return numberOfNonNullBins;
        }

        /**
         * @return An iterator over all non-empty bins.
         */
        public Iterator<Bin> iterator() {
            return new BinIterator();
        }

        private class BinIterator implements Iterator<Bin> {
            /**
             * Stores the bin # of the Bin currently in use.
             */
            private int nextBin;

            public BinIterator() {
                nextBin = 0;
            }

            /**
             * Are there more bins in this set, waiting to be returned?
             *
             * @return True if more bins are remaining.
             */
            public boolean hasNext() {
                while (nextBin <= maxBinNumber) {
                    if (getBin(nextBin) != null) return true;
                    nextBin++;
                }
                return false;
            }

            /**
             * Gets the next bin in the provided BinList.
             *
             * @return the next available bin in the BinList.
             */
            public Bin next() {
                if (!hasNext())
                    throw new NoSuchElementException("This BinIterator is currently empty");
                final Bin result = getBin(nextBin);
                nextBin++;
                return result;
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from a bin iterator");
            }
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final BinList bins = (BinList) o;

            if (maxBinNumber != bins.maxBinNumber) return false;
            if (numberOfNonNullBins != bins.numberOfNonNullBins) return false;
            if (!Arrays.equals(mBinArray, bins.mBinArray)) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(mBinArray);
            result = 31 * result + numberOfNonNullBins;
            result = 31 * result + maxBinNumber;
            return result;
        }
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final BinningIndexContent that = (BinningIndexContent) o;

        if (mReferenceSequence != that.mReferenceSequence) return false;
        if (!mBinList.equals(that.mBinList)) return false;
        if (!mLinearIndex.equals(that.mLinearIndex)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = mReferenceSequence;
        result = 31 * result + mBinList.hashCode();
        result = 31 * result + mLinearIndex.hashCode();
        return result;
    }
}
