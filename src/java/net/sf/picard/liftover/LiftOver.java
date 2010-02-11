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
import net.sf.picard.util.Log;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Java port of UCSC liftOver.  Only the most basic liftOver functionality is implemented.
 * Internally coordinates are 0-based, half-open. The API is standard Picard 1-based, inclusive.
 *
 * @author alecw@broadinstitute.org
 */
public class LiftOver {
    private static final Log LOG = Log.getInstance(LiftOver.class);
    
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
     * Lift over the given interval to the new genome build using the liftOverMinMatch set for this
     * LiftOver object.
     * @param interval Interval to be lifted over.
     * @return Interval in the output build coordinates, or null if it cannot be lifted over.
     */
    public Interval liftOver(final Interval interval) {
        return liftOver(interval, liftOverMinMatch);
    }

    /**
     * Lift over the given interval to the new genome build.
     * @param interval Interval to be lifted over.
     * @param liftOverMinMatch Minimum fraction of bases that must remap.
     * @return Interval in the output build coordinates, or null if it cannot be lifted over.
     */
    public Interval liftOver(final Interval interval, final double liftOverMinMatch) {
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
            } else if (candidateIntersection != null) {
                LOG.info("Interval " + interval.getName() + " failed to match chain " + chain.id +
                " because intersection length " + candidateIntersection.intersectionLength + " < minMatchSize "
                + minMatchSize +
                " (" + (candidateIntersection.intersectionLength/(float)interval.length()) + " < " + liftOverMinMatch + ")");
            }
        }
        if (chainHit == null) {
            // Can't be lifted over.
            return null;
        }

        return createToInterval(interval.getName(), targetIntersection);
    }

    public List<PartialLiftover> diagnosticLiftover(final Interval interval) {
        final List<PartialLiftover> ret = new ArrayList<PartialLiftover>();
        if (interval.length() == 0) {
            throw new IllegalArgumentException("Zero-length interval cannot be lifted over.  Interval: " +
                    interval.getName());
        }
        for (final Chain chain : chains.getOverlaps(interval)) {
            Interval intersectingChain = interval.intersect(chain.interval);
            final TargetIntersection targetIntersection = targetIntersection(chain, intersectingChain);
            if (targetIntersection == null) {
                ret.add(new PartialLiftover(intersectingChain, chain.id));
            } else {
                Interval toInterval = createToInterval(interval.getName(), targetIntersection);
                float percentLiftedOver = targetIntersection.intersectionLength/(float)interval.length();
                ret.add(new PartialLiftover(intersectingChain, toInterval, targetIntersection.chain.id, percentLiftedOver));
            }
        }
        return ret;
    }

    private static Interval createToInterval(final String intervalName, final TargetIntersection targetIntersection) {
        // Compute the query interval given the offsets of the target interval start and end into the first and
        // last ContinuousBlocks.
        int toStart = targetIntersection.chain.getBlock(targetIntersection.firstBlockIndex).toStart + targetIntersection.startOffset;
        int toEnd = targetIntersection.chain.getBlock(targetIntersection.lastBlockIndex).getToEnd() - targetIntersection.offsetFromEnd;
        if (toEnd <= toStart || toStart < 0) {
            throw new PicardException("Something strange lifting over interval " + intervalName);
        }

        if (targetIntersection.chain.toNegativeStrand) {
            // Flip if query is negative.
            int negativeStart = targetIntersection.chain.toSequenceSize - toEnd;
            int negativeEnd = targetIntersection.chain.toSequenceSize - toStart;
            toStart = negativeStart;
            toEnd = negativeEnd;
        }
        // Convert to 1-based, inclusive.
        return new Interval(targetIntersection.chain.toSequenceName, toStart+1, toEnd, targetIntersection.chain.toNegativeStrand,
                intervalName);
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
        return new TargetIntersection(chain, intersectionLength, startOffset, offsetFromEnd, firstBlockIndex, lastBlockIndex);
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
        /** Chain used for this intersection */
        final Chain chain;
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

        TargetIntersection(final Chain chain,final int intersectionLength, final int startOffset,
                           final int offsetFromEnd, final int firstBlockIndex, final int lastBlockIndex) {
            this.chain = chain;
            this.intersectionLength = intersectionLength;
            this.startOffset = startOffset;
            this.offsetFromEnd = offsetFromEnd;
            this.firstBlockIndex = firstBlockIndex;
            this.lastBlockIndex = lastBlockIndex;
        }
    }

    /**
     * Represents a portion of a liftover operation, for use in diagnosing liftover failures.
     */
    public static class PartialLiftover {
        /** Intersection between "from" interval and "from" region of a chain. */
        final Interval fromInterval;
        /**
         * Result of lifting over fromInterval (with no percentage mapped requirement).  This is null
         * if fromInterval falls entirely with a gap of the chain. */
        final Interval toInterval;
        /** id of chain used for this liftover */
        final int chainId;
        /** Percentage of bases in fromInterval that lifted over.  0 if fromInterval is not covered by any chain. */
        final float percentLiftedOver;

        PartialLiftover(final Interval fromInterval, final Interval toInterval, final int chainId, final float percentLiftedOver) {
            this.fromInterval = fromInterval;
            this.toInterval = toInterval;
            this.chainId = chainId;
            this.percentLiftedOver = percentLiftedOver;
        }

        PartialLiftover(final Interval fromInterval, final int chainId) {
            this.fromInterval = fromInterval;
            this.toInterval = null;
            this.chainId = chainId;
            this.percentLiftedOver = 0.0f;
        }

        public String toString() {
            if (toInterval == null) {
                // Matched a chain, but entirely within a gap.
                return fromInterval.toString() + " (len " + fromInterval.length() + ")=>null using chain " + chainId;
            }
            final String strand = toInterval.isNegativeStrand()? "-": "+";
            return fromInterval.toString() + " (len " + fromInterval.length() + ")=>" + toInterval + "(" + strand
                    + ") using chain " + chainId + " ; pct matched " + percentLiftedOver;
        }
    }
}
