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
package net.sf.picard.liftover;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Interval;
import net.sf.picard.util.OverlapDetector;

import java.io.*;
import java.util.List;

/**
 * Java port of UCSC liftOver.  Only the most basic liftOver functionality is implemented.
 * Internally coordinates are 0-based, half-open. The API is standard Picard 1-based, inclusive.
 *
 * @author alecw@broadinstitute.org
 */
public class LiftOver {
    public static final double DEFAULT_LIFTOVER_MINMATCH = 0.95;

    private double liftOverMinMatch = DEFAULT_LIFTOVER_MINMATCH;
    private final OverlapDetector<Chain> chains;

    /**
     * Load UCSC chain file in order to lift over Intervals.
     */
    public LiftOver(File chainFile) {
        IoUtil.assertFileIsReadable(chainFile);
        chains = Chain.loadChains(chainFile);
    }

    /**
     * Lift over the given interval to the new genome build.
     * @param interval Interval to be lifted over.
     * @return Interval in the output build coordinates, or null if it cannot be lifted over.
     */
    public Interval liftOver(final Interval interval) {
        if (interval.length() == 0) {
            throw new IllegalArgumentException("Zero-length interval cannot be lifted over.  Interval: " +
                    interval.getName());
        }
        Chain chainHit = null;
        TargetIntersection targetIntersection = null;
        // Number of bases in interval that can be lifted over must be >= this.
        double minMatchSize = liftOverMinMatch * interval.length();

        // Find the appropriate Chain, and the part of the chain corresponding to the interval to be lifted over.
        for (final Chain chain : chains.getOverlaps(interval)) {
            final TargetIntersection candidateIntersection = targetIntersection(chain, interval);
            if (candidateIntersection != null && candidateIntersection.intersectionLength >= minMatchSize) {
                if (chainHit != null) {
                    // In basic liftOver, multiple hits are not allowed.
                    return null;
                }
                chainHit = chain;
                targetIntersection = candidateIntersection;
            }
        }
        if (chainHit == null) {
            // Can't be lifted over.
            return null;
        }

        int toIntersectionLength = 0;
        for (int i = targetIntersection.firstBlockIndex; i <= targetIntersection.lastBlockIndex; ++i) {
            final Chain.ContinuousBlock block = chainHit.getBlock(i);
            toIntersectionLength += (block.getToEnd() - block.toStart);
        }
        toIntersectionLength -= (targetIntersection.startOffset + targetIntersection.offsetFromEnd);
        if (toIntersectionLength < 0) {
            throw new PicardException("Something strange lifting over interval " + interval.getName());
        }
        if (toIntersectionLength < minMatchSize) {
            // This probably won't happen, because the targetIntersection won't be big enough.
            // If this exception doesn't happen, the computation of toIntersectionLength can be eliminated.
            throw new PicardException("Something strange lifting over interval " + interval.getName());
        }
        // Compute the query interval given the offsets of the target interval start and end into the first and
        // last ContinuousBlocks.
        int toStart = chainHit.getBlock(targetIntersection.firstBlockIndex).toStart + targetIntersection.startOffset;
        int toEnd = chainHit.getBlock(targetIntersection.lastBlockIndex).getToEnd() - targetIntersection.offsetFromEnd;
        if (toEnd <= toStart || toStart < 0) {
            throw new PicardException("Something strange lifting over interval " + interval.getName());
        }

        if (chainHit.toNegativeStrand) {
            // Flip if query is negative.
            int negativeStart = chainHit.toSequenceSize - toEnd;
            int negativeEnd = chainHit.toSequenceSize - toStart;
            toStart = negativeStart;
            toEnd = negativeEnd;
        }
        // Convert to 1-based, inclusive.
        return new Interval(chainHit.toSequenceName, toStart+1, toEnd, chainHit.toNegativeStrand, interval.getName());
    }

    /**
     * Add up overlap btw the blocks in this chain and the given interval.
     * @return Length of overlap, offsets into first and last ContinuousBlocks, and indices of first and
     * last ContinuousBlocks.
     */
    private static TargetIntersection targetIntersection(final Chain chain, final Interval interval) {
        int intersectionLength = 0;
        // Convert interval to 0-based, half-open
        int start = interval.getStart() - 1;
        int end = interval.getEnd();
        int firstBlockIndex = -1;
        int lastBlockIndex = -1;
        int startOffset = -1;
        int offsetFromEnd = -1;
        List<Chain.ContinuousBlock> blockList = chain.getBlocks();
        for (int i = 0; i < blockList.size(); ++i) {
            final Chain.ContinuousBlock block = blockList.get(i);
            if (block.fromStart >= end) {
                break;
            } else if (block.getFromEnd() <= start) {
                continue;
            }
            if (firstBlockIndex == -1) {
                firstBlockIndex = i;
                if (start > block.fromStart) {
                    startOffset = start - block.fromStart;
                } else {
                    startOffset = 0;
                }
            }
            lastBlockIndex = i;
            if (block.getFromEnd() > end) {
                offsetFromEnd = block.getFromEnd() - end;
            } else {
                offsetFromEnd = 0;
            }
            int thisIntersection = Math.min(end, block.getFromEnd()) - Math.max(start, block.fromStart);
            if (thisIntersection <= 0) {
                throw new PicardException("Should have been some intersection.");
            }
            intersectionLength += thisIntersection;
        }
        if (intersectionLength == 0) {
            return null;
        }
        return new TargetIntersection(intersectionLength, startOffset, offsetFromEnd, firstBlockIndex, lastBlockIndex);
    }

    /**
     * Get minimum fraction of bases that must remap.
     */
    public double getLiftOverMinMatch() {
        return liftOverMinMatch;
    }

    /**
     * Set minimum fraction of bases that must remap.
     */
    public void setLiftOverMinMatch(final double liftOverMinMatch) {
        this.liftOverMinMatch = liftOverMinMatch;
    }

    /**
    * Value class returned by targetIntersection()
    */
    private static class TargetIntersection {
        /** Total intersectionLength length */
        final int intersectionLength;
        /** Offset of target interval start in first block. */
        final int startOffset;
        /** Distance from target interval end to end of last block. */
        final int offsetFromEnd;
        /** Index of first ContinuousBlock matching interval. */
        final int firstBlockIndex;
        /** Index of last ContinuousBlock matching interval. */
        final int lastBlockIndex;

        TargetIntersection(final int intersectionLength, final int startOffset, final int offsetFromEnd,
                           final int firstBlockIndex, final int lastBlockIndex) {
            this.intersectionLength = intersectionLength;
            this.startOffset = startOffset;
            this.offsetFromEnd = offsetFromEnd;
            this.firstBlockIndex = firstBlockIndex;
            this.lastBlockIndex = lastBlockIndex;
        }
    }
}
