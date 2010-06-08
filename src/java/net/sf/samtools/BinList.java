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

import java.util.Iterator;
import java.util.BitSet;
import java.util.NoSuchElementException;

/**
 * Provides a list of all bins which could exist in the BAM file.
 * Allows the user to iterate over all bins, selecting ones of interest
 * for later inspection.
 *
 * @author mhanna
 * @version 0.1
 */
public class BinList implements Iterable<Bin> {
    /**
     * The reference sequence relevant to this bin list.
     */
    private final int referenceSequence;

    /**
     * For each sequence, which bins should be included in the BitSet.
     */
    private final BitSet bins;

    /**
     * Create a new BinList over sequenceCount sequences, consisting of the given bins.
     * @param referenceSequence Reference sequence to which these bins are relevant.
     * @param bins The given bins to include.
     */
    protected BinList(final int referenceSequence, final BitSet bins) {
        this.referenceSequence = referenceSequence;
        this.bins = bins;
    }

    /**
     * Gets an iterator over all selected bins.
     * @return An iterator over all selected bins.
     */
    public Iterator<Bin> iterator() {
        return new BinIterator();
    }

    /**
     * Get the reference sequence to which this bin belongs.
     * @return Integer representing the reference sequence.
     */
    protected int getReferenceSequence() {
        return referenceSequence;
    }

    /**
     * Retrieves the bins stored in this list.
     * @return A bitset where a bin is present in the list if the bit is true.
     */
    protected BitSet getBins() {
        return bins;
    }

    private class BinIterator implements Iterator<Bin> {
        /**
         * Stores the bin currently in use.  Will be -1 if no more bins remain in the set.
         */
        private int nextBin;

        public BinIterator() {
            // Initialize the bin iterator to just before the first bin.
            nextBin = bins.nextSetBit(0);
        }

        /**
         * Are there more bins in this set, waiting to be returned?
         * @return True if more bins are remaining.
         */
        public boolean hasNext() {
            return nextBin >= 0;
        }

        /**
         * Gets the next bin in the provided BinList.
         * @return the next available bin in the BinList.
         */
        public Bin next() {
            if(!hasNext())
                throw new NoSuchElementException("This BinIterator is currently empty");
            int currentBin = nextBin;
            nextBin = bins.nextSetBit(nextBin+1);
            return new Bin(referenceSequence,currentBin);
        }

        public void remove() {
            throw new UnsupportedOperationException("Unable to remove from a bin iterator");
        }
    }
}

