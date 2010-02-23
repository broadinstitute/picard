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
package net.sf.picard.util;

import net.sf.samtools.SAMFileHeader;

import java.util.BitSet;
import java.util.List;

/**
 * Serve up loci of interest based on an interval list.
 *
 * @author alecw at broadinstitute dot oh are gee
 */
public class IntervalListReferenceSequenceMask implements ReferenceSequenceMask {

    private final SAMFileHeader header;
    // if memory usage becomes a problem... this could be changed to a SparseBitSet
    // http://java.sun.com/developer/onlineTraining/collections/magercises/BitSet/index.html
    private final BitSet currentBitSet = new BitSet();
    private int currentSequenceIndex = -1;
    private final PeekableIterator<Interval> intervalIterator;
    private final int lastSequenceIndex;
    private final int lastPosition;

    public IntervalListReferenceSequenceMask(final IntervalList intervalList) {
        this.header = intervalList.getHeader();
        if (intervalList.getHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            intervalList.sort();
        }
        final List<Interval> uniqueIntervals = intervalList.getUniqueIntervals();
        if (uniqueIntervals.isEmpty()) {
            lastSequenceIndex = -1;
            lastPosition = 0;
        } else {
            final Interval lastInterval = uniqueIntervals.get(uniqueIntervals.size() - 1);
            lastSequenceIndex = header.getSequenceIndex((lastInterval.getSequence()));
            lastPosition = lastInterval.getEnd();
        }
        intervalIterator = new PeekableIterator<Interval>(uniqueIntervals.iterator());
    }

    /**
     * It is required that sequenceIndex is >= any previous sequenceIndex passed to this class.
     * @return true if the mask is set for the given sequence and position
     */
    public boolean get(final int sequenceIndex, final int position) {
        ensureSequenceLoaded(sequenceIndex);
        return currentBitSet.get(position);
    }

    /**
     * It is required that sequenceIndex is >= any previous sequenceIndex passed to this class.
     * @return the next pos on the given sequence >= position that is set, or -1 if there are no more set positions
     */
    public int nextPosition(final int sequenceIndex, final int position) {
        ensureSequenceLoaded(sequenceIndex);
        return currentBitSet.nextSetBit(position);
    }

    private void ensureSequenceLoaded(final int sequenceIndex) {
        if (sequenceIndex < this.currentSequenceIndex) {
            throw new IllegalArgumentException("Cannot look at an earlier sequence.  Current: " +
                    this.currentSequenceIndex + "; requested: " + sequenceIndex);
        }
        if (sequenceIndex > currentSequenceIndex) {
            currentBitSet.clear();
            while (intervalIterator.hasNext()) {
                final Interval interval = intervalIterator.peek();
                final int nextSequenceIndex = header.getSequenceIndex(interval.getSequence());
                if (nextSequenceIndex < sequenceIndex) {
                    intervalIterator.next();
                } else if (nextSequenceIndex == sequenceIndex) {
                    currentBitSet.set(interval.getStart(), interval.getEnd() + 1);
                    intervalIterator.next();
                } else {
                    break;
                }
            }
            currentSequenceIndex = sequenceIndex;
        }
    }

    /**
     * @return Largest sequence index for which there are set bits.
     */
    public int getMaxSequenceIndex() {
        return lastSequenceIndex;
    }

    /**
     * @return the largest position on the last sequence index
     */
    public int getMaxPosition() {
        return lastPosition;
    }

    public SAMFileHeader getHeader() {
        return header;
    }
}
