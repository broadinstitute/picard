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

package picard.analysis.directed;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.MetricAccumulationLevel;
import picard.filter.CountingMapQFilter;
import picard.metrics.MultilevelMetrics;
import picard.metrics.PerUnitMetricCollector;
import picard.metrics.SAMRecordMultiLevelCollector;
import picard.util.MathUtil;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * TargetMetrics, are metrics to measure how well we hit specific targets (or baits) when using a targeted sequencing process like hybrid selection
 * or Targeted PCR Techniques (TSCA).  TargetMetrics at the moment are the metrics that are shared by both HybridSelection and TargetedPcrMetrics.
 *
 * TargetMetricsCollector collects for a run these common metrics and can be sub-classed to provide metrics more specific to a targeted sequencing
 * run.
 *
 * Note: Probe is the name I've used to indicate the bait set or amplicon set (e.g. the individual technological units used to target specific
 * sites).
 *
 * @author Jonathan Burke
 */
public abstract class TargetMetricsCollector<METRIC_TYPE extends MultilevelMetrics> extends SAMRecordMultiLevelCollector<METRIC_TYPE, Integer> {

    // What is considered "near" to the bait
    private static final int NEAR_PROBE_DISTANCE = 250;

    //If perTargetCoverage != null then coverage is computed for each specified target and output to this file
    private final File perTargetCoverage;

    //The name of the set of probes used
    private final String probeSetName;

    private static final Log log = Log.getInstance(TargetMetricsCollector.class);

    //The interval list indicating the regions targeted by all probes
    private final IntervalList allProbes;

    //The interval list of the the regions we intend to cover
    private final IntervalList allTargets;

    // Overlap detector for finding overlaps between reads and the experimental targets
    private final OverlapDetector<Interval> targetDetector;

    // Overlap detector for finding overlaps between the reads and the baits (and the near bait space)
    private final OverlapDetector<Interval> probeDetector;

    private Map<Interval,Double> intervalToGc = null;

    //The number of bases within all unique intervals in allProbes
    private final long probeTerritory;

    //The number of bases within all unique intervals found in allTargets
    private final long targetTerritory;

    private final long genomeSize;

    private final int minimumMappingQuality;

    private final int minimumBaseQuality;

    private boolean noSideEffects;

    //A map of coverage by target in which coverage is reset every read, this is done
    //so that we can calculate overlap for a read once and the resulting coverage is
    //than added to the cumulative coverage of every collector that collects
    //information on that read
    private Map<Interval, Coverage> coverageByTargetForRead;
    private Coverage [] cov;

    //Converts a targetMetric into a more specific metric of METRIC_TYPE
    public abstract METRIC_TYPE convertMetric(final TargetMetrics targetMetrics);


    /**
     * In the case of ignoring bases in overlapping reads from the same template,
     * we can choose whether to modify the existing record, or modify a clone.  The
     * latter may significantly slow down the performance, but also has no side effects.
     * @param value the boolean value to set.
     */
    public void setNoSideEffects(final boolean value) {
        this.noSideEffects = value;
    }

    /**
     * Since the targeted metrics (HsMetrics, TargetedPcrMetrics,...) share many of the same values as TargetMetrics, this copy will copy all public attributes in targetMetrics
     * to the outputMetrics' attributes of the same name.  If no matching attribute exists in the outputMetrics or the attribute of the target metrics class also is found
     * in targetKeys then it's value is not copied.  Further more, targetKeys and outputKeys are attribute name arrays synchronized by the index.
     * For each target key, targetMetrics.<targetKeys[i]> is assigned to outputMetrics.<outputKeys[i]>
     *
     * @param targetMetrics A metric with values to be copied
     * @param outputMetrics A metrics intended to receive values from targetMetrics
     * @param targetKeys Specific names of attributes of targetMetrics to copy to outputMetrics, each key has a corresponding one in outputKeys
     * @param outputKeys Specific names of the destination attributes of outputMetrics that will be filled with values of outputMetrics, each key has a corresponding one in targetKeys
     * @param <MT> The type of metric of outputMetrics
     */
    protected static <MT extends MetricBase> void reflectiveCopy(final TargetMetrics targetMetrics, final MT outputMetrics, final String [] targetKeys, final String [] outputKeys) {

        if(targetKeys == null || outputKeys == null) {
            if(outputKeys != null) {
                throw new PicardException("Target keys is null but output keys == " + StringUtil.join(",", outputKeys));
            }

            if(targetKeys != null) {
                throw new PicardException("Output keys is null but target keys == " + StringUtil.join(",", targetKeys));
            }
        } else {
            if(targetKeys.length != outputKeys.length) {
                throw new PicardException("Target keys and output keys do not have the same length: " +
                        "targetKeys == (" + StringUtil.join(",", targetKeys) + ") " +
                        "outputKeys == (" + StringUtil.join(",", outputKeys) + ")");
            }
        }

        final Class mtClass = outputMetrics.getClass();
        final Set<Field> targetSet = CollectionUtil.makeSet(TargetMetrics.class.getFields());

        for(final String targetKey : targetKeys) {
            if(targetSet.contains(targetKey)) {
                targetSet.remove(targetKey);
            }
        }

        final Set<String> outputSet = new HashSet<String>();
        for(final Field field : outputMetrics.getClass().getFields()) {
            outputSet.add(field.getName());
        }

        for(final Field field : targetSet) {
            if(outputSet.contains(field.getName())) {
                try {
                    final Field outputField = mtClass.getField(field.getName());
                    outputField.set(outputMetrics, field.get(targetMetrics));
                } catch (Exception e) {
                    throw new PicardException("Exception while copying targetMetrics to " + outputMetrics.getClass().getName(), e);
                }
            }
        }

        for(int i = 0; i < targetKeys.length; i++) {
            try {
                Field targetMetricField = TargetMetrics.class.getField(targetKeys[i]);
                Field outputMetricField = mtClass.getField(outputKeys[i]);
                outputMetricField.set(outputMetrics, targetMetricField.get(targetMetrics));
            } catch(final Exception exc) {
                throw new PicardException("Exception while copying TargetMetrics." + targetKeys[i] + " to " + mtClass.getName() + "." + outputKeys[i], exc);
            }
        }
    }

    public TargetMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                  final List<SAMReadGroupRecord> samRgRecords,
                                  final ReferenceSequenceFile refFile,
                                  final File perTargetCoverage,
                                  final IntervalList targetIntervals,
                                  final IntervalList probeIntervals,
                                  final String probeSetName,
                                  final int minimumMappingQuality,
                                  final int minimumBaseQuality) {
        this(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName, minimumMappingQuality, minimumBaseQuality, false);
    }

    public TargetMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                  final List<SAMReadGroupRecord> samRgRecords,
                                  final ReferenceSequenceFile refFile,
                                  final File perTargetCoverage,
                                  final IntervalList targetIntervals,
                                  final IntervalList probeIntervals,
                                  final String probeSetName,
                                  final int minimumMappingQuality,
                                  final int minimumBaseQuality,
                                  final boolean noSideEffects) {
        this.perTargetCoverage = perTargetCoverage;
        this.probeSetName = probeSetName;

        this.allProbes  = probeIntervals;
        this.allTargets = targetIntervals;

        final List<Interval> uniqueBaits = this.allProbes.uniqued().getIntervals();
        this.probeDetector = new OverlapDetector<Interval>(-NEAR_PROBE_DISTANCE, 0);
        this.probeDetector.addAll(uniqueBaits, uniqueBaits);
        this.probeTerritory = Interval.countBases(uniqueBaits);

        final List<Interval> uniqueTargets = this.allTargets.uniqued().getIntervals();
        targetDetector = new OverlapDetector<Interval>(0,0);
        this.targetDetector.addAll(uniqueTargets, uniqueTargets);
        this.targetTerritory = Interval.countBases(uniqueTargets);

        // Populate the coverage by target map
        int i = 0;
        cov = new Coverage[uniqueTargets.size()];
        this.coverageByTargetForRead = new LinkedHashMap<Interval, Coverage>(uniqueTargets.size() * 2, 0.5f);
        for (final Interval target : uniqueTargets) {
            final Coverage coverage = new Coverage(target, 0);
            this.coverageByTargetForRead.put(target, coverage);
            cov[i++] = coverage;
        }

        long genomeSizeAccumulator = 0;
        for (final SAMSequenceRecord seq : this.allProbes.getHeader().getSequenceDictionary().getSequences()) {
            genomeSizeAccumulator += seq.getSequenceLength();
        }
        this.genomeSize = genomeSizeAccumulator;


        if (refFile != null) {
            intervalToGc = new HashMap<Interval,Double>();
            for (final Interval target : uniqueTargets) {
                final ReferenceSequence rs = refFile.getSubsequenceAt(target.getSequence(), target.getStart(), target.getEnd());
                intervalToGc.put(target,SequenceUtil.calculateGc(rs.getBases()));
            }
        }


        this.minimumMappingQuality = minimumMappingQuality;

        this.minimumBaseQuality = minimumBaseQuality;

        this.noSideEffects = noSideEffects;

        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
        final PerUnitTargetMetricCollector collector =  new PerUnitTargetMetricCollector(probeSetName, coverageByTargetForRead.keySet(),
                                                                                         sample, library, readGroup, probeTerritory, targetTerritory, genomeSize,
                                                                                         intervalToGc, minimumMappingQuality, minimumBaseQuality);
        if (this.probeSetName != null) {
            collector.setBaitSetName(probeSetName);
        }

        return collector;
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeAllReadCollector() {
        final PerUnitTargetMetricCollector collector = (PerUnitTargetMetricCollector) makeChildCollector(null, null, null);
        if (perTargetCoverage != null) {
            collector.setPerTargetOutput(perTargetCoverage);
        }

        return collector;
    }

    /**
     * Collect the Target Metrics for one unit of "accumulation" (i.e. for one sample, or for one library ...)
     */
    public class PerUnitTargetMetricCollector implements PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> {

        private final Map<Interval,Double> intervalToGc;
        private File perTargetOutput;

        // A Map to accumulate per-bait-region (i.e. merge of overlapping targets) coverage. */
        private final Map<Interval, Coverage> coverageByTarget;

        private final TargetMetrics metrics = new TargetMetrics();


        // Filters that should be applied before accepting a record.  Order of filters in this list
        // matches the order applied.
        final List<SamRecordFilter> filters = new ArrayList<SamRecordFilter>();

        private final int minimumBaseQuality;
        private final CountingMapQFilter mapQFilter;

        /**
         * Constructor that parses the squashed reference to genome reference file and stores the
         * information in a map for later use.
         */
        public PerUnitTargetMetricCollector(final String probeSetName, final Set<Interval> coverageTargets,
                                            final String sample, final String library, final String readGroup,
                                            final long probeTerritory, final long targetTerritory, final long genomeSize,
                                            final Map<Interval, Double> intervalToGc,
                                            final int minimumMappingQuality,
                                            final int minimumBaseQuality) {
            this.metrics.SAMPLE           = sample;
            this.metrics.LIBRARY          = library;
            this.metrics.READ_GROUP       = readGroup;
            this.metrics.PROBE_SET        = probeSetName;

            metrics.PROBE_TERRITORY  = probeTerritory;
            metrics.TARGET_TERRITORY = targetTerritory;
            metrics.GENOME_SIZE      = genomeSize;

            this.coverageByTarget = new LinkedHashMap<Interval, Coverage>(coverageTargets.size() * 2, 0.5f);
            for (Interval target : coverageTargets) {
                this.coverageByTarget.put(target, new Coverage(target,0));
            }

            this.mapQFilter = new CountingMapQFilter(minimumMappingQuality);

            this.filters.add(new SecondaryAlignmentFilter());
            this.filters.add(mapQFilter);

            this.minimumBaseQuality = minimumBaseQuality;

            this.intervalToGc = intervalToGc;
        }

        /** If set, the metrics collector will output per target coverage information to this file. */
        public void setPerTargetOutput(final File perTargetOutput) {
            this.perTargetOutput = perTargetOutput;
        }

        /** Sets the name of the bait set explicitly instead of inferring it from the bait file. */
        public void setBaitSetName(final String name) {
            this.metrics.PROBE_SET = name;
        }

        private int getBasesAtMinimumQualityInBlock(final SAMRecord record, final AlignmentBlock block) {
            int basesInBLockAtMinimumQuality = 0;
            final byte[] baseQualities = record.getBaseQualities();
            for (int idx = block.getReadStart(); idx <= CoordMath.getEnd(block.getLength(), block.getReadStart()); idx++) { // idx is one-based
                if (minimumBaseQuality <= baseQualities[idx-1]) basesInBLockAtMinimumQuality++;
            }
            return basesInBLockAtMinimumQuality;
        }

        private SAMRecord clipOverlappingBases(final SAMRecord rec) {
            // NB: ignores how to handle supplemental records when present for both ends by just using the mate information in the record.

            if (!rec.getReadPairedFlag()) return rec;

            // Only clip records that are left-most in genomic order and overlapping and are the first end
            if (rec.getSecondOfPairFlag()) return rec; // second end, ignore.
            if (rec.getMateAlignmentStart() < rec.getAlignmentStart()) return rec; // right-most, so ignore.  Do this first since getAlignmentEnd is expensive
            final int numBasesToClip = rec.getMateAlignmentStart() - rec.getAlignmentEnd() + 1;
            if (0 < numBasesToClip) return rec; // left-most but not overlapping

            try {
                // watch out for when the second read overlaps all of the first read
                if (rec.getMateAlignmentStart() <= rec.getAlignmentStart()) { // make it unmapped
                    if (noSideEffects) {
                        final SAMRecord recClone = ((SAMRecord)rec.clone());
                        recClone.setReadUnmappedFlag(true);
                        return recClone;
                    }
                    else rec.setReadUnmappedFlag(true);
                    return rec;
                }

                // clip it, clip it good
                final Cigar cigar = CigarUtil.addSoftClippedBasesToEndsOfCigar(
                        rec.getCigar(),
                        false, // ignore, since we want 3' end genomic order
                        numBasesToClip, // want 3' end genomic order
                        0
                );
                if (noSideEffects) {
                    final SAMRecord recClone = ((SAMRecord)rec.clone());
                    recClone.setCigar(cigar);
                    return recClone;
                }
                else {
                    rec.setCigar(cigar);
                    return rec;
                }
            } catch (CloneNotSupportedException e) {
                throw new PicardException(e.getMessage(), e);
            }
        }

        /** Adds information about an individual SAMRecord to the statistics. */
        public void acceptRecord(final SAMRecord record) {
            // We should not "Just plain avoid records that are marked as not-primary" since:
            // 1. We should disregard secondary alignments in for both record and base counting.
            // 2. We should disregard supplementary alignments for ONLY record counting
            // This means that we need to consider supplemental alignments for base counting but not record counting.

            // NB: this could modify the record.  See noSideEffects
            final SAMRecord rec = clipOverlappingBases(record);

            // apply any filters
            for (final SamRecordFilter filter : this.filters) {
                if (filter.filterOut(rec)) return;
            }

            if (!rec.getSupplementaryAlignmentFlag()) this.metrics.TOTAL_READS += 1;

            // Check for PF reads
            if (rec.getReadFailsVendorQualityCheckFlag()) {
                return;
            }

            // Prefetch the list of target and bait overlaps here as they're needed multiple times.
            final Collection<Interval> targets;
            final Collection<Interval> probes;

            if (!rec.getReadUnmappedFlag()) {
                final Interval read = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
                targets = targetDetector.getOverlaps(read);
                probes   = probeDetector.getOverlaps(read);
            }
            else {
                targets = null;
                probes = null;
            }

            // get the number of bases above the minimum base quality
            int basesAtMinimumQuality = 0;
            for (final byte baseQuality : rec.getBaseQualities()) {
                if (minimumBaseQuality <= baseQuality) basesAtMinimumQuality++;
                else metrics.PCT_EXC_BASEQ++;
            }

            if (!rec.getSupplementaryAlignmentFlag()) {
                ++this.metrics.PF_READS;
                this.metrics.PF_BASES += basesAtMinimumQuality;

                // And now calculate the values we need for HS_LIBRARY_SIZE
                if (rec.getReadPairedFlag() && rec.getFirstOfPairFlag() && !rec.getReadUnmappedFlag() && !rec.getMateUnmappedFlag()) {
                    if (probes != null && !probes.isEmpty()) {
                        ++this.metrics.PF_SELECTED_PAIRS;
                        if (!rec.getDuplicateReadFlag()) ++this.metrics.PF_SELECTED_UNIQUE_PAIRS;
                    }
                }
            }

            // Check for reads that are marked as duplicates
            if (rec.getDuplicateReadFlag()) {
                return;
            }
            else if (!rec.getSupplementaryAlignmentFlag()) {
                ++this.metrics.PF_UNIQUE_READS;
            }

            // Don't bother with reads that didn't align uniquely
            if (rec.getReadUnmappedFlag() || rec.getMappingQuality() == 0) {
                return;
            }

            if (!rec.getSupplementaryAlignmentFlag()) this.metrics.PF_UQ_READS_ALIGNED += 1;
            for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                this.metrics.PF_UQ_BASES_ALIGNED += getBasesAtMinimumQualityInBlock(rec, block);
            }

            final boolean mappedInPair = rec.getReadPairedFlag() && !rec.getMateUnmappedFlag() && !rec.getSupplementaryAlignmentFlag();

            final byte[] baseQualities = rec.getBaseQualities();

            // Find the target overlaps
            if (targets != null && !targets.isEmpty()) {
                for (final Interval target : targets) {
                    final Coverage coverage = this.coverageByTarget.get(target);

                    for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                        final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                        int readPos = block.getReadStart(); // one-based
                        for (int pos=block.getReferenceStart(); pos<=end; ++ pos, ++readPos) {
                            if (pos >= target.getStart() && pos <= target.getEnd()) {
                                if (baseQualities[readPos-1] < minimumBaseQuality) continue; // ignore low-quality bases
                                ++this.metrics.ON_TARGET_BASES;
                                if (mappedInPair) ++this.metrics.ON_TARGET_FROM_PAIR_BASES;
                                coverage.addBase(pos - target.getStart());
                            }
                        }
                    }
                }
            }

            // Now do the bait overlaps
            int mappedBases = 0;
            for (final AlignmentBlock block : rec.getAlignmentBlocks()) mappedBases += getBasesAtMinimumQualityInBlock(rec, block);
            int onBaitBases = 0;

            if (probes != null && !probes.isEmpty()) {
                for (final Interval bait : probes) {
                    for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                        final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                        int readPos = block.getReadStart(); // one-based
                        for (int pos=block.getReferenceStart(); pos<=end; ++pos, ++readPos) {
                            if (baseQualities[readPos-1] < minimumBaseQuality) continue; // ignore low-quality bases
                            if (pos >= bait.getStart() && pos <= bait.getEnd()) ++onBaitBases;
                        }
                    }
                }

                this.metrics.ON_PROBE_BASES   += onBaitBases;
                this.metrics.NEAR_PROBE_BASES += (mappedBases - onBaitBases);
            }
            else {
                this.metrics.OFF_PROBE_BASES += mappedBases;
            }

        }

        @Override
        public void finish() {
            metrics.PCT_PF_READS         = metrics.PF_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS      = metrics.PF_UNIQUE_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS_ALIGNED = metrics.PF_UQ_READS_ALIGNED / (double) metrics.PF_UNIQUE_READS;

            final double denominator   = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES + metrics.OFF_PROBE_BASES);

            metrics.PCT_SELECTED_BASES = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES) / denominator;
            metrics.PCT_OFF_PROBE         = metrics.OFF_PROBE_BASES / denominator;
            metrics.ON_PROBE_VS_SELECTED = metrics.ON_PROBE_BASES / (double) (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES);
            metrics.MEAN_PROBE_COVERAGE   = metrics.ON_PROBE_BASES / (double) metrics.PROBE_TERRITORY;
            metrics.FOLD_ENRICHMENT       = (metrics.ON_PROBE_BASES/ denominator) / ((double) metrics.PROBE_TERRITORY / metrics.GENOME_SIZE);

            metrics.PCT_EXC_MAPQ = mapQFilter.getFilteredBases();

            calculateTargetCoverageMetrics();
            calculateGcMetrics();
        }

        /** Calculates how much additional sequencing is needed to raise 80% of bases to the mean for the lane. */
        private void calculateTargetCoverageMetrics() {
            final int[] depths = new int[(int) this.metrics.TARGET_TERRITORY];  // may not use entire array
            int zeroCoverageTargets = 0;
            int depthIndex = 0;
            double totalCoverage = 0;
            int basesConsidered = 0;
            final Histogram<Integer> targetCoverageHistogram = new Histogram<Integer>();

            for (final Coverage c : this.coverageByTarget.values()) {
                if (!c.hasCoverage()) {
                    ++zeroCoverageTargets;
                    continue;
                }

                final int[] targetDepths = c.getDepths();
                basesConsidered += targetDepths.length;

                for (final int depth : targetDepths) {
                    depths[depthIndex++] = depth;
                    totalCoverage += depth;
                    targetCoverageHistogram.increment(depth);
                }
            }
            targetCoverageHistogram.increment(0, zeroCoverageTargets); // do not forget zero coverage targets, and also forces this histogram to have a value for the number of 0x target bases

            this.metrics.MEAN_TARGET_COVERAGE = totalCoverage / basesConsidered;
            this.metrics.MEDIAN_TARGET_COVERAGE = targetCoverageHistogram.getMedian();

            // get the number of bases at zero coverage, and subtract that from the total coverage, then divide that by the total coverage
            // to get the percent of target bases at 1x coverage or greater.
            final double numTargetBasesAtZeroCoverage = targetCoverageHistogram.get(0).getValue();
            this.metrics.PCT_TARGET_BASES_1X = (this.metrics.TARGET_TERRITORY - numTargetBasesAtZeroCoverage) / this.metrics.TARGET_TERRITORY;

            // Sort the array (ASCENDING) and then find the base the coverage value that lies at the 80%
            // line, which is actually at 20% into the array now
            Arrays.sort(depths);
            // Note.  basesConsidered can be between 0 and depths.length inclusive.  indexOf80thPercentile will be -1 in the latter case
            final int indexOf80thPercentile = Math.max((depths.length - 1 - basesConsidered) + (int) (basesConsidered * 0.2), 0);
            final int coverageAt80thPercentile;
            if(depths.length > 0) {
                coverageAt80thPercentile = depths[indexOf80thPercentile];
            }
            else {
                throw new PicardException("Interval list only contains one zero-length interval.");
            }
            this.metrics.FOLD_80_BASE_PENALTY = this.metrics.MEAN_TARGET_COVERAGE / coverageAt80thPercentile;
            this.metrics.ZERO_CVG_TARGETS_PCT = zeroCoverageTargets / (double) allTargets.getIntervals().size();

            // Now do the "how many bases at X" calculations.
            int totalTargetBases = 0;
            int targetBases2x  = 0;
            int targetBases10x = 0;
            int targetBases20x = 0;
	        int targetBases30x = 0;
	        int targetBases40x = 0;
	        int targetBases50x = 0;
	        int targetBases100x = 0;

            for (final Coverage c : this.coverageByTarget.values()) {
                for (final int depth : c.getDepths()) {
                    ++totalTargetBases;

                    if (depth >= 2) {
                        ++targetBases2x;
                        if (depth >=10) {
                            ++targetBases10x;
                            if (depth >= 20) {
                                ++targetBases20x;
                                if (depth >=30) {
                                    ++targetBases30x;
	                                if (depth >=40) {
		                                ++targetBases40x;
		                                if (depth >=50) {
			                                ++targetBases50x;
			                                if (depth >=100) {
				                                ++targetBases100x;
			                                }
		                                }
	                                }
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
	        this.metrics.PCT_TARGET_BASES_40X = (double) targetBases40x / (double) totalTargetBases;
	        this.metrics.PCT_TARGET_BASES_50X = (double) targetBases50x / (double) totalTargetBases;
	        this.metrics.PCT_TARGET_BASES_100X = (double) targetBases100x / (double) totalTargetBases;
        }

        private void calculateGcMetrics() {
            if (this.intervalToGc != null) {
                log.info("Calculating GC metrics");

                // Setup the output file if we're outputting per-target coverage
                FormatUtil fmt = new FormatUtil();
                final PrintWriter out;
                try {
                    if (perTargetOutput != null) {
                        out = new PrintWriter(perTargetOutput);
                        out.println("chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage\tnormalized_coverage\tmin_normalized_coverage");
                    }
                    else {
                        out = null;
                    }
                }
                catch (IOException ioe) { throw new RuntimeIOException(ioe); }

                final int bins = 101;
                final long[] targetBasesByGc  = new long[bins];
                final long[] alignedBasesByGc = new long[bins];

                for (final Map.Entry<Interval,Coverage> entry : this.coverageByTarget.entrySet()) {
                    final Interval interval = entry.getKey();
                    final Coverage cov = entry.getValue();

                    if (interval.length() <= 0) {
                        log.warn("interval of length zero found: " + interval + " skipped.");
                        continue;
                    }

                    final double gcDouble = this.intervalToGc.get(interval);
                    final int gc = (int) Math.round(gcDouble * 100);

                    targetBasesByGc[gc]  += interval.length();
                    alignedBasesByGc[gc] += cov.getTotal();

                    if (out != null) {
                        final double coverage = cov.getTotal() / (double) interval.length();
                        final double min = MathUtil.min(cov.getDepths());

                        out.println(interval.getSequence() + "\t" +
                                    interval.getStart() + "\t" +
                                    interval.getEnd() + "\t" +
                                    interval.length() + "\t" +
                                    interval.getName() + "\t" +
                                    fmt.format(gcDouble) + "\t" +
                                    fmt.format(coverage) + "\t" +
                                    fmt.format(coverage / this.metrics.MEAN_TARGET_COVERAGE) + "\t" +
                                    fmt.format(min / this.metrics.MEAN_TARGET_COVERAGE)
                        );
                    }
                }

                if (out != null) out.close();

                // Total things up
                long totalTarget = 0;
                long totalBases  = 0;
                for (int i=0; i<targetBasesByGc.length; ++i) {
                    totalTarget += targetBasesByGc[i];
                    totalBases  += alignedBasesByGc[i];
                }

                // Re-express things as % of the totals and calculate dropout metrics
                for (int i=0; i<targetBasesByGc.length; ++i) {
                    final double targetPct  = targetBasesByGc[i]  / (double) totalTarget;
                    final double alignedPct = alignedBasesByGc[i] / (double) totalBases;

                    double dropout = (alignedPct - targetPct) * 100d;
                    if (dropout < 0) {
                        dropout = Math.abs(dropout);

                        if (i <=50) this.metrics.AT_DROPOUT += dropout;
                        if (i >=50) this.metrics.GC_DROPOUT += dropout;
                    }
                }
            }
        }


        @Override
        public void addMetricsToFile(MetricsFile<METRIC_TYPE, Integer> hsMetricsComparableMetricsFile) {
            hsMetricsComparableMetricsFile.addMetric(convertMetric(this.metrics));
        }
    }

    /**
     * A simple class that is used to store the coverage information about an interval.
     *
     * @author Tim Fennell
     */
    public static class Coverage {
        private final Interval interval;
        private final int[] depths;

        /** Constructs a new coverage object for the provided mapping with the desired padding either side. */
        public Coverage(final Interval i, final int padding) {
            this.interval = i;
            this.depths = new int[interval.length() + 2*padding];
        }

        /** Adds a single point of depth at the desired offset into the coverage array. */
        public void addBase(final int offset) {
            if (offset >= 0 && offset < this.depths.length) {
                if (this.depths[offset] < Integer.MAX_VALUE) {
                    this.depths[offset] += 1;
                }
            }
        }

        /** Returns true if any base in the range has coverage of > 1 */
        public boolean hasCoverage() {
            for (final int s : depths) {
                if (s > 1) return true;
            }

            return false;
        }

        /** Gets the coverage depths as an array of ints. */
        public int[] getDepths() { return this.depths; }

        public int getTotal() {
            int total = 0;
            for (int i=0; i<depths.length; ++i) total += depths[i];
            return total;
        }

        @Override
        public String toString() {
            return "TargetedMetricCollector(interval=" + interval + ", depths = [" + StringUtil.intValuesToString(this.depths) + "])";
        }
    }
}

/**
 * For a sequencing run targeting specific regions of the genome this metric class holds metrics describing
 * how well those regions were targeted.
 */
class TargetMetrics extends MultilevelMetrics {
    /**  The name of the PROBE_SET (BAIT SET, AMPLICON SET, ...) used in this metrics collection run */
    public String PROBE_SET;

    /** The number of unique bases covered by the intervals of all probes in the probe set */
    public long PROBE_TERRITORY;

    /** The number of unique bases covered by the intervals of all targets that should be covered */
    public long TARGET_TERRITORY;

    /** The number of bases in the reference genome used for alignment. */
    public long GENOME_SIZE;

    /** The total number of reads in the SAM or BAM file examined. */
    public long TOTAL_READS;

    /** The number of reads that pass the vendor's filter. */
    public long PF_READS;

    /** The number of bases in the SAM or BAM file to be examined */
    public long PF_BASES;

    /** The number of PF reads that are not marked as duplicates. */
    public long PF_UNIQUE_READS;

    // Tracks the number of read pairs that we see that are PF (used to calculate library size) */
    public long PF_SELECTED_PAIRS;

    // Tracks the number of unique PF reads pairs we see (used to calc library size)
    public long PF_SELECTED_UNIQUE_PAIRS;

    /** The number of PF unique reads that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** The number of PF unique bases that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF aligned probed that mapped to a baited region of the genome. */
    public long ON_PROBE_BASES;

    /** The number of PF aligned bases that mapped to within a fixed interval of a probed region, but not on a baited region. */
    public long NEAR_PROBE_BASES;

    /** The number of PF aligned bases that mapped to neither on or near a probe. */
    public long OFF_PROBE_BASES;

    /** The number of PF aligned bases that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    /** The number of PF aligned bases that are mapped in pair to a targeted region of the genome. */
    public long ON_TARGET_FROM_PAIR_BASES;

    /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
    public double PCT_EXC_BASEQ;

    //metrics below here are derived after collection

    /** PF reads / total reads.  The percent of reads passing filter. */
    public double PCT_PF_READS;

    /** PF Unique Reads / Total Reads. */
    public double PCT_PF_UQ_READS;

    /** PF Reads Aligned / PF Reads. */
    public double PCT_PF_UQ_READS_ALIGNED;

    /** On+Near Bait Bases / PF Bases Aligned. */
    public double PCT_SELECTED_BASES;

    /** The percentage of aligned PF bases that mapped neither on or near a probe. */
    public double PCT_OFF_PROBE;

    /** The percentage of on+near probe bases that are on as opposed to near. */
    public double ON_PROBE_VS_SELECTED;

    /** The mean coverage of all probes in the experiment. */
    public double MEAN_PROBE_COVERAGE;

    /** The fold by which the probed region has been amplified above genomic background. */
    public double FOLD_ENRICHMENT;

    /** The mean coverage of targets that received at least coverage depth = 2 at one base. */
    public double MEAN_TARGET_COVERAGE;

    /** The mean coverage of targets that received at least coverage depth = 2 at one base. */
    public double MEDIAN_TARGET_COVERAGE;

    /** The number of targets that did not reach coverage=2 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
    public double PCT_EXC_MAPQ;

    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;

    /** The percentage of ALL target bases achieving 1X or greater coverage. */
    public double PCT_TARGET_BASES_1X;
    /** The percentage of ALL target bases achieving 2X or greater coverage. */
    public double PCT_TARGET_BASES_2X;
    /** The percentage of ALL target bases achieving 10X or greater coverage. */
    public double PCT_TARGET_BASES_10X;
    /** The percentage of ALL target bases achieving 20X or greater coverage. */
    public double PCT_TARGET_BASES_20X;
	/** The percentage of ALL target bases achieving 30X or greater coverage. */
	public double PCT_TARGET_BASES_30X;
	/** The percentage of ALL target bases achieving 40X or greater coverage. */
	public double PCT_TARGET_BASES_40X;
	/** The percentage of ALL target bases achieving 50X or greater coverage. */
	public double PCT_TARGET_BASES_50X;
	/** The percentage of ALL target bases achieving 100X or greater coverage. */
	public double PCT_TARGET_BASES_100X;

    /**
     * A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC<=50% regions mapped elsewhere.
     */
    public double AT_DROPOUT;

    /**
     * A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100]
     * we calculate a = % of target territory, and b = % of aligned reads aligned to these targets.
     * GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total
     * reads that should have mapped to GC>=50% regions mapped elsewhere.
     */
    public double GC_DROPOUT;
}
