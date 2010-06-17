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

package net.sf.picard.analysis.directed;

import net.sf.picard.sam.DuplicationMetrics;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.Log;

import java.util.*;
import java.io.*;

import net.sf.picard.util.OverlapDetector;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.util.CoordMath;

/**
 * Calculates HS metrics for a given SAM or BAM file. Requires the input of a list of
 * target intervals and a list of bait intervals. Can be invoked either on an entire
 * iterator of SAMRecords or be passed SAMRecords one at a time.
 *
 * @author Tim Fennell
 */
public class HsMetricsCalculator {
    // What is considered "near" to the bait
    private static final int NEAR_BAIT_DISTANCE = 250;
    private static final Log log = Log.getInstance(HsMetricsCalculator.class);

    // Holds file names and other parameter related junk
    private SAMFileReader sam;
    private final File baitFile;
    private final File targetFile;
    private final IntervalList baits;
    private final IntervalList targets;

    // Overlap detector for finding overlaps between reads and the experimental targets
    private final OverlapDetector<Interval> targetDetector = new OverlapDetector<Interval>(0,0);

	// Overlap detector for finding overlaps between the reads and the baits (and the near bait space)
    private final OverlapDetector<Interval> baitDetector = new OverlapDetector<Interval>(-NEAR_BAIT_DISTANCE,0);

    // A Map to accumulate per-bait-region (i.e. merge of overlapping baits) coverage. */
    private final Map<Interval, Coverage> coverageByTarget;

    private final HsMetrics metrics = new HsMetrics();
    private double PF_BASES = 0; // Not reported in HsMetrics but needed for calculations
    private long PF_SELECTED_PAIRS; // Tracks the number of read pairs that we see that are PF (used to calc hs library size)
    private long PF_SELECTED_UNIQUE_PAIRS; // Tracks the number of unique PF reads pairs we see (used to calc hs library size)
    private long ON_TARGET_FROM_PAIR_BASES = 0;

    /**
	 * Constructor that parses the squashed reference to genome reference file and stores the
	 * information in a map for later use.
	 */
    public HsMetricsCalculator(final File baits, final File targets) {
        this.baitFile     = baits;
        this.targetFile   = targets;
        this.baits        = IntervalList.fromFile(baits);
        this.targets      = IntervalList.fromFile(targets);

        this.metrics.BAIT_SET = baits.getName();
        final int tmp = this.metrics.BAIT_SET.lastIndexOf(".");
        if (tmp > 0) {
            this.metrics.BAIT_SET = this.metrics.BAIT_SET.substring(0, tmp);
        }

        final List<Interval> uniqueBaits = this.baits.getUniqueIntervals();
        this.baitDetector.addAll(uniqueBaits, uniqueBaits);
        this.metrics.BAIT_TERRITORY = Interval.countBases(uniqueBaits);

        final List<Interval> uniqueTargets = this.targets.getUniqueIntervals();
        this.targetDetector.addAll(uniqueTargets, uniqueTargets);
        this.metrics.TARGET_TERRITORY = Interval.countBases(uniqueTargets);

        for (final SAMSequenceRecord seq : this.baits.getHeader().getSequenceDictionary().getSequences()) {
            this.metrics.GENOME_SIZE += seq.getSequenceLength();
        }

        // Populate the coverage by target map
        this.coverageByTarget = new HashMap<Interval, Coverage>(uniqueTargets.size() * 2, 0.5f);
        for (final Interval target : uniqueTargets) {
            this.coverageByTarget.put(target, new Coverage(target, 0));
        }
    }

    /** Iterates over all records in the file and collects metrics. */
    public void analyze(final Iterator<SAMRecord> records) {
        int i = 0;
        while (records.hasNext()) {
            analyze(records.next());

            if (++i % 1000000 == 0) {
                log.info("Processed " + i + " records so far.");
            }
        }
    }

    /** Adds information about an individual SAMRecord to the statistics. */
    public void analyze(final SAMRecord rec) {
        // Just plain avoid records that are marked as not-primary
        if (rec.getNotPrimaryAlignmentFlag()) return;
        
        this.metrics.TOTAL_READS += 1;

        // Prefetch the list of target and bait overlaps here as they're needed multiple times.
        final Collection<Interval> targets;
        final Collection<Interval> baits;

        if (!rec.getReadUnmappedFlag()) {
            final Interval read = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
            targets = this.targetDetector.getOverlaps(read);
            baits = this.baitDetector.getOverlaps(read);
        }
        else {
            targets = null;
            baits = null;
        }

        // Check for PF reads
        if (rec.getReadFailsVendorQualityCheckFlag()) {
            return;
        }
        else {
            ++this.metrics.PF_READS;
            this.PF_BASES += rec.getReadLength();

            // And now calculate the values we need for HS_LIBRARY_SIZE
            if (rec.getReadPairedFlag() && rec.getFirstOfPairFlag() && !rec.getReadUnmappedFlag() && !rec.getMateUnmappedFlag()) {
                if (baits != null && !baits.isEmpty()) {
                    ++this.PF_SELECTED_PAIRS;
                    if (!rec.getDuplicateReadFlag()) ++this.PF_SELECTED_UNIQUE_PAIRS;
                }
            }
        }

        // Check for reads that are marked as duplicates
        if (rec.getDuplicateReadFlag()) {
            return;
        }
        else {
            ++this.metrics.PF_UNIQUE_READS;
        }

        // Don't bother with reads that didn't align uniquely
        if (rec.getReadUnmappedFlag() || rec.getMappingQuality() == 0) {
            return;
        }

        this.metrics.PF_UQ_READS_ALIGNED += 1;
        for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
            this.metrics.PF_UQ_BASES_ALIGNED += block.getLength();
        }

        final boolean mappedInPair = rec.getReadPairedFlag() && !rec.getMateUnmappedFlag();

        // Find the target overlaps
        if (targets != null && !targets.isEmpty()) {
            for (final Interval target : targets) {
                final Coverage coverage = this.coverageByTarget.get(target);

                for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                    final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                    for (int pos=block.getReferenceStart(); pos<=end; ++ pos) {
                        if (pos >= target.getStart() && pos <= target.getEnd()) {
                            ++this.metrics.ON_TARGET_BASES;
                            if (mappedInPair) ++this.ON_TARGET_FROM_PAIR_BASES;
                            coverage.addBase(pos - target.getStart());
                        }
                    }
                }
            }
        }

        // Now do the bait overlaps
        int mappedBases = 0;
        for (final AlignmentBlock block : rec.getAlignmentBlocks()) mappedBases += block.getLength();
        int onBaitBases = 0;

        if (baits != null && !baits.isEmpty()) {
            for (final Interval bait : baits) {
                for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                    final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());

                    for (int pos=block.getReferenceStart(); pos<=end; ++pos) {
                        if (pos >= bait.getStart() && pos <= bait.getEnd()) ++onBaitBases;
                    }
                }
            }

            this.metrics.ON_BAIT_BASES   += onBaitBases;
            this.metrics.NEAR_BAIT_BASES += (mappedBases - onBaitBases);            
        }
        else {
            this.metrics.OFF_BAIT_BASES += mappedBases;
        }

    }

    /** Calculates a few last summary metrics and then returns the metrics calculated. */
    public HsMetrics getMetrics() {
        this.metrics.PCT_USABLE_BASES_ON_BAIT   = this.metrics.ON_BAIT_BASES / this.PF_BASES;
        this.metrics.PCT_USABLE_BASES_ON_TARGET = this.metrics.ON_TARGET_BASES / this.PF_BASES;
        this.metrics.HS_LIBRARY_SIZE = DuplicationMetrics.estimateLibrarySize(this.PF_SELECTED_PAIRS, this.PF_SELECTED_UNIQUE_PAIRS);

        this.metrics.calculateDerivedMetrics();
        calculateTargetCoverageMetrics();

        this.metrics.HS_PENALTY_10X = calculateHsPenalty(10);
        this.metrics.HS_PENALTY_20X = calculateHsPenalty(20);
        this.metrics.HS_PENALTY_30X = calculateHsPenalty(30);

        return this.metrics;
    }

    /** Calculates how much additional sequencing is needed to raise 80% of bases to the mean for the lane. */
    private void calculateTargetCoverageMetrics() {
        final short[] depths = new short[(int) this.metrics.TARGET_TERRITORY];  // may not use entire array
        int zeroCoverageTargets = 0;
        int depthIndex = 0;
        double totalCoverage = 0;
        int basesConsidered = 0;

        for (final Coverage c : this.coverageByTarget.values()) {
            if (!c.hasCoverage()) {
                ++zeroCoverageTargets;
                continue;
            }

            final short[] targetDepths = c.getDepths();
            basesConsidered += targetDepths.length;

            for (final short depth : targetDepths) {
                depths[depthIndex++] = depth;
                totalCoverage += depth;
            }
        }

        this.metrics.MEAN_TARGET_COVERAGE = totalCoverage / basesConsidered;

        // Sort the array (ASCENDING) and then find the base the coverage value that lies at the 80%
        // line, which is actually at 20% into the array now
        Arrays.sort(depths);
        final int indexOf80thPercentile = (depths.length - 1 - basesConsidered) + (int) (basesConsidered * 0.2);
        final int coverageAt80thPercentile = depths[indexOf80thPercentile];
        this.metrics.FOLD_80_BASE_PENALTY = this.metrics.MEAN_TARGET_COVERAGE / coverageAt80thPercentile;
        this.metrics.ZERO_CVG_TARGETS_PCT = zeroCoverageTargets / (double) this.targets.getIntervals().size();

        // Now do the "how many bases at X" calculations.
        int totalTargetBases = 0;
        int targetBases2x  = 0;
        int targetBases10x = 0;
        int targetBases20x = 0;
        int targetBases30x = 0;

        for (final Coverage c : this.coverageByTarget.values()) {
            for (short depth : c.getDepths()) {
                ++totalTargetBases;

                if (depth >= 2) {
                    ++targetBases2x;
                    if (depth >=10) {
                        ++targetBases10x;
                        if (depth >= 20) {
                            ++targetBases20x;
                            if (depth >=30) {
                                ++targetBases30x;
                            }
                        }
                    }
                }
            }
        }

        this.metrics.PCT_TARGET_BASES_2X  = (double) targetBases2x  / (double) totalTargetBases;
        this.metrics.PCT_TARGET_BASES_10X = (double) targetBases10x / (double) totalTargetBases;
        this.metrics.PCT_TARGET_BASES_20X = (double) targetBases20x / (double) totalTargetBases;
        this.metrics.PCT_TARGET_BASES_30X = (double) targetBases30x / (double) totalTargetBases;
    }

    /**
     * Attempts to calculate the HS penalty incurred by the library in order to get 80%
     * of target bases (in non-zero-covered targets) to a specific target coverage (e.g. 20X).
     *
     * @param coverageGoal the desired coverage target (e.g. 20X)
     * @return the hs penalty - a multiplier that tells if you want, e.g. 20X coverage, then you will
     *         need to produce this many PF aligned bases per target bases in your design!
     */
    private double calculateHsPenalty(final int coverageGoal) {
        if (this.metrics.HS_LIBRARY_SIZE == null) return 0;

        final double meanCoverage  = this.ON_TARGET_FROM_PAIR_BASES / (double) this.metrics.TARGET_TERRITORY;
        final double fold80        = this.metrics.FOLD_80_BASE_PENALTY;
        final long hsLibrarySize   = this.metrics.HS_LIBRARY_SIZE;
        final long pairs           = this.PF_SELECTED_PAIRS;
        final long uniquePairs     = this.PF_SELECTED_UNIQUE_PAIRS;
        final double onTargetPct   = (double) this.metrics.ON_TARGET_BASES / (double) this.metrics.PF_UQ_BASES_ALIGNED;

        final double uniquePairGoalMultiplier = (coverageGoal / meanCoverage) * fold80;
        double pairMultiplier = uniquePairGoalMultiplier;
        double increment = 1;
        boolean goingUp = uniquePairGoalMultiplier >= 1;
        double finalPairMultiplier = -1;

        // Converge "pairMultiplier" to the number that gives us a uniquePairMultiplier equal
        // to the coverage multiplier we desire.  If we can't get there with 1000X coverage,
        // we're not going to get there!
        for (int i=0; i<10000; ++i) {
            double uniquePairMultiplier = DuplicationMetrics.estimateRoi(hsLibrarySize, pairMultiplier, pairs, uniquePairs);

            if (Math.abs(uniquePairMultiplier - uniquePairGoalMultiplier) <= 0.01) {
                finalPairMultiplier  = pairMultiplier;
                break;
            }
            else if ((uniquePairMultiplier > uniquePairGoalMultiplier && goingUp) ||
                     (uniquePairMultiplier < uniquePairGoalMultiplier && !goingUp)){
                increment /= 2;
                goingUp = !goingUp;
            }

            pairMultiplier += (goingUp ? increment : -increment);
        }

        if (finalPairMultiplier == -1) {
            return -1;
        }
        else {
            double uniqueFraction = (uniquePairs * uniquePairGoalMultiplier) / (pairs * finalPairMultiplier);
            return (1 / uniqueFraction) * fold80 * (1 / onTargetPct);
        }
    }

    /**
         * A simple class that is used to store the coverage information about an interval.
         *
         * @author Tim Fennell
         */
        public static class Coverage {
            private final Interval interval;
            private final short[] depths;

            /** Constructs a new coverage object for the provided mapping with the desired padding either side. */
            public Coverage(final Interval i, final int padding) {
                this.interval = i;
                this.depths = new short[interval.length() + 2*padding];
            }

            /** Adds a single point of depth at the desired offset into the coverage array. */
            public void addBase(final int offset) {
                if (offset >= 0 && offset < this.depths.length) {
                    this.depths[offset] += 1;
                }
            }

            /** Returns true if any base in the range has coverage of > 1 */
            public boolean hasCoverage() {
                for (final short s : depths) {
                    if (s > 1) return true;
                }

                return false;
            }

            /** Gets the coverage depths as an array of shorts. */
            public short[] getDepths() { return this.depths; }
        }
}
