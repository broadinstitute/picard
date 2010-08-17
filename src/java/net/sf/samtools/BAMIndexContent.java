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

/**
 * Represents the contents of a bam index file for one reference.
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
    private final BinList mBinList;

    /**
     * Chunks containing metaData for the reference, e.g. number of aligned and unaligned records
     */
    private final BAMIndexMetaData mMetaData;

    /**
     * The linear index for the reference sequence above.
     */
    private final LinearIndex mLinearIndex;


    /**
     * @param referenceSequence Content corresponds to this reference.
     * @param bins              Array of bins represented by this content, possibly sparse
     * @param numberOfBins      Number of non-null bins
     * @param metaData          Extra information about the reference in this index
     * @param linearIndex       Additional index used to optimize queries
     */
    BAMIndexContent(final int referenceSequence, final Bin[] bins, final int numberOfBins, final BAMIndexMetaData metaData, final LinearIndex linearIndex) {
        this.mReferenceSequence = referenceSequence;
        this.mBinList = new BinList(bins, numberOfBins);
        this.mMetaData = metaData;
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
     * @return the meta data chunks for this content
     */
    public BAMIndexMetaData getMetaData() {
        return mMetaData;
    }

    /**
     * @return all chunks associated with all bins in this content
     */
    public List<Chunk> getAllChunks() {
        List<Chunk> allChunks = new ArrayList<Chunk>();
        for (Bin b : mBinList)
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
     * This class is used to encapsulate the list of Bins store in the BAMIndexContent
     * While it is currently represented as an array, we may decide to change it to an ArrayList or other structure
     */
    class BinList implements Iterable<Bin> {

        private final Bin[] mBinArray;
        public final int numberOfNonNullBins;
        public final int maxBinNumber;  // invariant: maxBinNumber = mBinArray.length -1 since array is 0 based

        /**
         * @param binArray            a sparse array representation of the bins. The index into the array is the bin number.
         * @param numberOfNonNullBins
         */
        BinList(Bin[] binArray, int numberOfNonNullBins) {
            this.mBinArray = binArray;
            this.numberOfNonNullBins = numberOfNonNullBins;
            this.maxBinNumber = mBinArray.length - 1;
        }

        Bin getBin(int binNumber) {
            if (binNumber > maxBinNumber) return null;
            return mBinArray[binNumber];
        }

        int getNumberOfNonNullBins() {
            return numberOfNonNullBins;
        }

        /**
         * Gets an iterator over all non-null bins.
         *
         * @return An iterator over all bins.
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
                Bin result = getBin(nextBin);
                nextBin++;
                return result;
            }

            public void remove() {
                throw new UnsupportedOperationException("Unable to remove from a bin iterator");
            }
        }
    }
}
