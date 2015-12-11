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
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FormatUtil;
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
    /** Default distance for a read to be considered "selected". */
    public static final int NEAR_PROBE_DISTANCE_DEFAULT = 250;
    private int nearProbeDistance = NEAR_PROBE_DISTANCE_DEFAULT;

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
    private final boolean clipOverlappingReads;
    private boolean noSideEffects;

    //A map of coverage by target in which coverage is reset every read, this is done
    //so that we can calculate overlap for a read once and the resulting coverage is
    //than added to the cumulative coverage of every collector that collects
    //information on that read
    private final Map<Interval, Coverage> coverageByTargetForRead;
    private final Coverage [] cov;

    /** Gets the distance that is allowed between a read and the nearest probe for it to be considered "near probe" and "selected. */
    public int getNearProbeDistance() { return nearProbeDistance; }

    /** Sets the distance that is allowed between a read and the nearest probe for it to be considered "near probe" and "selected. */
    public void setNearProbeDistance(final int nearProbeDistance) { this.nearProbeDistance = nearProbeDistance; }

    //Converts a targetMetric into a more specific metric of METRIC_TYPE
    public abstract METRIC_TYPE convertMetric(final TargetMetrics targetMetrics);


    /**
     * In the case of ignoring bases in overlapping reads from the same template,
     * we choose to internally modify the SAM record's CIGAR to clip overlapping bases.
     * We can either to modify the passed-in record (a side-effect), or modify a
     * an internally clone of the record (no side-effect).  Due to the overhead of
     * cloning a SAMRecord object, we may see significant slow down of the
     * performance to ensure there are no side effects.  Therefore, users of this
     * collector who do not care if the record passed to
     * {@link #acceptRecord(htsjdk.samtools.SAMRecord, htsjdk.samtools.reference.ReferenceSequence)}
     * is modified can pass in false to this method to realize performance gains.
     * @param value the boolean value to set.
     */
    public void setNoSideEffects(final boolean value) {
        this.noSideEffects = value;
    }

    /** Get the the number of bases in the given alignment block and record that have base quality greater or equal to the minimum */
    public static int getNumBasesPassingMinimumBaseQuality(final SAMRecord record, final AlignmentBlock block, final int minimumBaseQuality) {
        int basesInBlockAtMinimumQuality = 0;
        final byte[] baseQualities = record.getBaseQualities();
        for (int idx = block.getReadStart(); idx <= CoordMath.getEnd(block.getLength(), block.getReadStart()); idx++) { // idx is one-based
            if (minimumBaseQuality <= baseQualities[idx-1]) basesInBlockAtMinimumQuality++;
        }
        return basesInBlockAtMinimumQuality;
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
                } catch (final Exception e) {
                    throw new PicardException("Exception while copying targetMetrics to " + outputMetrics.getClass().getName(), e);
                }
            }
        }

        for(int i = 0; i < targetKeys.length; i++) {
            try {
                final Field targetMetricField = TargetMetrics.class.getField(targetKeys[i]);
                final Field outputMetricField = mtClass.getField(outputKeys[i]);
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
                                  final int nearProbeDistance,
                                  final int minimumMappingQuality,
                                  final int minimumBaseQuality,
                                  final boolean clipOverlappingReads) {
        this(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance, minimumMappingQuality, minimumBaseQuality, clipOverlappingReads, false);
    }

    public TargetMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                  final List<SAMReadGroupRecord> samRgRecords,
                                  final ReferenceSequenceFile refFile,
                                  final File perTargetCoverage,
                                  final IntervalList targetIntervals,
                                  final IntervalList probeIntervals,
                                  final String probeSetName,
                                  final int nearProbeDistance,
                                  final int minimumMappingQuality,
                                  final int minimumBaseQuality,
                                  final boolean clipOverlappingReads,
                                  final boolean noSideEffects) {
        this.perTargetCoverage = perTargetCoverage;
        this.probeSetName = probeSetName;
        this.nearProbeDistance = nearProbeDistance;

        this.allProbes  = probeIntervals;
        this.allTargets = targetIntervals;

        final List<Interval> uniqueBaits = this.allProbes.uniqued().getIntervals();
        this.probeDetector = new OverlapDetector<Interval>(-this.nearProbeDistance, 0);
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
        this.clipOverlappingReads = clipOverlappingReads;
        this.noSideEffects = noSideEffects;

        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
        final PerUnitTargetMetricCollector collector =  new PerUnitTargetMetricCollector(probeSetName, coverageByTargetForRead.keySet(),
                                                                                         sample, library, readGroup, probeTerritory, targetTerritory, genomeSize,
                                                                                         intervalToGc, minimumMappingQuality, minimumBaseQuality, clipOverlappingReads);
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

        private final int minimumBaseQuality;
        private final CountingMapQFilter mapQFilter;
        private final boolean clipOverlappingReads;

        /**
         * Constructor that parses the squashed reference to genome reference file and stores the
         * information in a map for later use.
         */
        public PerUnitTargetMetricCollector(final String probeSetName, final Set<Interval> coverageTargets,
                                            final String sample, final String library, final String readGroup,
                                            final long probeTerritory, final long targetTerritory, final long genomeSize,
                                            final Map<Interval, Double> intervalToGc,
                                            final int minimumMappingQuality,
                                            final int minimumBaseQuality,
                                            final boolean clipOverlappingReads) {
            this.metrics.SAMPLE           = sample;
            this.metrics.LIBRARY          = library;
            this.metrics.READ_GROUP       = readGroup;
            this.metrics.PROBE_SET        = probeSetName;

            metrics.PROBE_TERRITORY  = probeTerritory;
            metrics.TARGET_TERRITORY = targetTerritory;
            metrics.GENOME_SIZE      = genomeSize;

            this.coverageByTarget = new LinkedHashMap<Interval, Coverage>(coverageTargets.size() * 2, 0.5f);
            for (final Interval target : coverageTargets) {
                this.coverageByTarget.put(target, new Coverage(target,0));
            }

            this.mapQFilter = new CountingMapQFilter(minimumMappingQuality);
            this.minimumBaseQuality = minimumBaseQuality;
            this.intervalToGc = intervalToGc;
            this.clipOverlappingReads = clipOverlappingReads;
        }

        /** If set, the metrics collector will output per target coverage information to this file. */
        public void setPerTargetOutput(final File perTargetOutput) {
            this.perTargetOutput = perTargetOutput;
        }

        /** Sets the name of the bait set explicitly instead of inferring it from the bait file. */
        public void setBaitSetName(final String name) {
            this.metrics.PROBE_SET = name;
        }

        /** Adds information about an individual SAMRecord to the statistics. */
        public void acceptRecord(final SAMRecord record) {
            // Just ignore secondary alignments altogether
            if (record.getNotPrimaryAlignmentFlag()) return;

            // Cache some things, and compute the total number of bases aligned in the record.
            final boolean mappedInPair = record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && !record.getSupplementaryAlignmentFlag();
            final byte[] baseQualities = record.getBaseQualities();
            int basesAlignedInRecord = 0;
            if (!record.getReadUnmappedFlag()) {
                for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                    basesAlignedInRecord += block.getLength();
                }
            }

            /*
            Count metrics related to # of base and reads. Consider supplemental alignments for base counting but not record counting.
            Do this counting *prior* to applying filters to ensure we match other metrics' computation for these values, such as AlignmentSummaryMetrics.
            */
            if (!record.getNotPrimaryAlignmentFlag()) { // ignore secondary alignments for both # of read and # of base metrics
                // # of bases
                if (!record.getReadFailsVendorQualityCheckFlag()) { // only reads that pass vendor's filters
                    // Strangely enough we should not count supplementals in PF_BASES, assuming that the
                    // main record also contains these bases! But we *do* count the aligned bases, assuming
                    // that those bases are not *aligned* in the primary record
                    if (!record.getSupplementaryAlignmentFlag()) this.metrics.PF_BASES += record.getReadLength();

                    if (!record.getReadUnmappedFlag()) {
                        this.metrics.PF_BASES_ALIGNED += basesAlignedInRecord;
                        if (!record.getDuplicateReadFlag()) {
                            this.metrics.PF_UQ_BASES_ALIGNED += basesAlignedInRecord;
                        }
                    }
                }

                // # of reads
                if (!record.getSupplementaryAlignmentFlag()) { // only consider the primary
                    this.metrics.TOTAL_READS++;
                    if (!record.getReadFailsVendorQualityCheckFlag()) { // only reads that pass vendor's filters
                        this.metrics.PF_READS++;
                        if (!record.getDuplicateReadFlag()) { // ignore duplicates for unique reads/bases
                            this.metrics.PF_UNIQUE_READS++;
                            if (!record.getReadUnmappedFlag()) { // ignore unmapped reads
                                this.metrics.PF_UQ_READS_ALIGNED++;
                            }
                        }
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////
            // Unmapped reads can be totally ignored beyond this point
            ///////////////////////////////////////////////////////////////////
            if (record.getReadUnmappedFlag()) return;

            // Prefetch the list of target and bait overlaps here as they're needed multiple times.
            final Interval read = new Interval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd());
            final Collection<Interval> targets = targetDetector.getOverlaps(read);
            final Collection<Interval> probes  = probeDetector.getOverlaps(read);

            // Calculate the values we need for HS_LIBRARY_SIZE
            if (!record.getSupplementaryAlignmentFlag()) {
                if (record.getReadPairedFlag() && record.getFirstOfPairFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag()) {
                    if (!probes.isEmpty()) {
                        ++this.metrics.PF_SELECTED_PAIRS;
                        if (!record.getDuplicateReadFlag()) ++this.metrics.PF_SELECTED_UNIQUE_PAIRS;
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////
            // Duplicate reads can be totally ignored beyond this point
            ///////////////////////////////////////////////////////////////////
            if (record.getDuplicateReadFlag()) {
                this.metrics.PCT_EXC_DUPE += basesAlignedInRecord;
                return;
            }

            // Compute the bait-related metrics *before* applying the overlap clipping and
            // the map-q threshold, since those would skew the assay-related metrics
            {
                final int mappedBases = basesAlignedInRecord;
                int onBaitBases = 0;

                if (!probes.isEmpty()) {
                    for (final Interval bait : probes) {
                        for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                            final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());

                            for (int pos=block.getReferenceStart(); pos<=end; ++pos) {
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

            ///////////////////////////////////////////////////////////////////
            // And lastly, ignore reads falling below the mapq threshold
            ///////////////////////////////////////////////////////////////////
            if (this.mapQFilter.filterOut(record)) return;

            // NB: this could modify the record.  See noSideEffects.
            final SAMRecord rec;
            if (clipOverlappingReads) {
                final int numOverlappingBasesToClip = SAMUtils.getNumOverlappingAlignedBasesToClip(record);
                rec = SAMUtils.clipOverlappingAlignedBases(record, numOverlappingBasesToClip, noSideEffects);
                metrics.PCT_EXC_OVERLAP += numOverlappingBasesToClip;
            }
            else rec = record;

            // Find the target overlaps
            for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                final int length = block.getLength(), refStart = block.getReferenceStart(), readStart = block.getReadStart();

                for (int offset = 0; offset < length; ++offset) {
                    final int refPos = refStart + offset;
                    final int readPos = readStart + offset;
                    final int qual = baseQualities[readPos - 1];

                    if (qual < minimumBaseQuality) {
                        this.metrics.PCT_EXC_BASEQ++;
                    }
                    else {
                        boolean isOnTarget = false;
                        for (final Interval target : targets) {
                            if (refPos >= target.getStart() && refPos<= target.getEnd()) {
                                isOnTarget = true;
                                ++this.metrics.ON_TARGET_BASES;
                                if (mappedInPair) ++this.metrics.ON_TARGET_FROM_PAIR_BASES;

                                final int targetOffset = refPos - target.getStart();
                                final Coverage coverage = this.coverageByTarget.get(target);
                                coverage.addBase(targetOffset);
                            }
                        }

                        if (!isOnTarget) this.metrics.PCT_EXC_OFF_TARGET++;
                    }
                }
            }
        }

        @Override
        public void finish() {
            metrics.PCT_PF_READS         = metrics.PF_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS      = metrics.PF_UNIQUE_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS_ALIGNED = metrics.PF_UQ_READS_ALIGNED / (double) metrics.PF_UNIQUE_READS;

            final double denominator   = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES + metrics.OFF_PROBE_BASES);

            metrics.PCT_SELECTED_BASES   = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES) / denominator;
            metrics.PCT_OFF_PROBE        = metrics.OFF_PROBE_BASES / denominator;
            metrics.ON_PROBE_VS_SELECTED = metrics.ON_PROBE_BASES / (double) (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES);
            metrics.MEAN_PROBE_COVERAGE  = metrics.ON_PROBE_BASES / (double) metrics.PROBE_TERRITORY;
            metrics.FOLD_ENRICHMENT      = (metrics.ON_PROBE_BASES/ denominator) / ((double) metrics.PROBE_TERRITORY / metrics.GENOME_SIZE);

            metrics.PCT_EXC_DUPE       /= (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_MAPQ        = mapQFilter.getFilteredBases() / (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_BASEQ      /= (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_OVERLAP    /= (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_OFF_TARGET /= (double) metrics.PF_BASES_ALIGNED;

            calculateTargetCoverageMetrics();
            calculateGcMetrics();
        }

        /** Calculates how much additional sequencing is needed to raise 80% of bases to the mean for the lane. */
        private void calculateTargetCoverageMetrics() {
            final int[] depths = new int[(int) this.metrics.TARGET_TERRITORY];
            int zeroCoverageTargets = 0;
            int depthIndex = 0;
            double totalCoverage = 0;

            // The "how many bases at at-least X" calculations.
            final int targetBasesDepth[] = {0, 1, 2, 10, 20, 30, 40, 50, 100}; // NB: this should be in ascending order
            final int targetBases[] = new int[targetBasesDepth.length]; // counts for how many target bases are at at least X coverage, where X corresponds to the value at the same offset in targetBasesDepth

            // consider all bases in calculating the mean, median etc.
            for (final Coverage c : this.coverageByTarget.values()) {
                final boolean hasCoverage = c.hasCoverage();
                final int[] targetDepths = c.getDepths();

                if (!hasCoverage) zeroCoverageTargets++;

                for (final int depth : targetDepths) {
                    if (0 < depth) totalCoverage += depth;
                    if (hasCoverage) depths[depthIndex++] = depth;

                    // Add to the "how many bases at at-least X" calculations.
                    for (int i = 0; i < targetBasesDepth.length; i++) {
                        if (depth >= targetBasesDepth[i]) targetBases[i]++;
                        else break; // NB: assumes that targetBasesDepth is sorted in ascending order
                    }
                }
            }

            // Sort the array (ASCENDING)
            Arrays.sort(depths);

            this.metrics.MEAN_TARGET_COVERAGE = totalCoverage / this.metrics.TARGET_TERRITORY;
            this.metrics.MEDIAN_TARGET_COVERAGE =
                    (0 == (depths.length % 2)) ? ((depths[(depths.length/2) - 1] + depths[(depths.length/2)]) / 2.0) : depths[(depths.length-1)/2];

            // find the base the coverage value that lies at the 80%
            // line, which is actually at 20% into the array now
            final int indexOf80thPercentile = Math.max((int) (depthIndex * 0.2 ) - 1, 0);
            final int coverageAt80thPercentile;
            if (depths.length > 0) {
                coverageAt80thPercentile = depths[indexOf80thPercentile];
            }
            else {
                throw new PicardException("Interval list only contains one zero-length interval.");
            }
            this.metrics.FOLD_80_BASE_PENALTY = this.metrics.MEAN_TARGET_COVERAGE / coverageAt80thPercentile;
            this.metrics.ZERO_CVG_TARGETS_PCT = zeroCoverageTargets / (double) allTargets.getIntervals().size();

            // Store the "how many bases at at-least X" calculations.
            this.metrics.PCT_TARGET_BASES_1X   = (double) targetBases[1] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_2X   = (double) targetBases[2] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_10X  = (double) targetBases[3] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_20X  = (double) targetBases[4] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_30X  = (double) targetBases[5] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_40X  = (double) targetBases[6] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_50X  = (double) targetBases[7] / (double) targetBases[0];
            this.metrics.PCT_TARGET_BASES_100X = (double) targetBases[8] / (double) targetBases[0];
        }

        private void calculateGcMetrics() {
            if (this.intervalToGc != null) {
                log.info("Calculating GC metrics");

                // Setup the output file if we're outputting per-target coverage
                final FormatUtil fmt = new FormatUtil();
                final PrintWriter out;
                try {
                    if (perTargetOutput != null) {
                        out = new PrintWriter(perTargetOutput);
                        out.println("chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage\tnormalized_coverage\tmin_normalized_coverage\tmax_normalized_coverage\tmin_coverage\tmax_coverage\tpct_0x");
                    }
                    else {
                        out = null;
                    }
                }
                catch (final IOException ioe) { throw new RuntimeIOException(ioe); }

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
                        double min = Integer.MAX_VALUE;
                        double max = Integer.MIN_VALUE;
                        double targetBasesAt0x = 0.0;
                        for (final int d : cov.getDepths()) {
                            if (0 == d) targetBasesAt0x++;
                            if (d < min) min = d;
                            if (max < d) max = d;
                        }

                        out.println(interval.getSequence() + "\t" +
                                    interval.getStart() + "\t" +
                                    interval.getEnd() + "\t" +
                                    interval.length() + "\t" +
                                    interval.getName() + "\t" +
                                    fmt.format(gcDouble) + "\t" +
                                    fmt.format(coverage) + "\t" +
                                    fmt.format(coverage / this.metrics.MEAN_TARGET_COVERAGE) + "\t" +
                                    fmt.format(min / this.metrics.MEAN_TARGET_COVERAGE) + "\t" +
                                    fmt.format(max / this.metrics.MEAN_TARGET_COVERAGE) + "\t" +
                                    fmt.format(min) + "\t" +
                                    fmt.format(max) + "\t" +
                                    fmt.format(targetBasesAt0x / interval.length())
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
        public void addMetricsToFile(final MetricsFile<METRIC_TYPE, Integer> hsMetricsComparableMetricsFile) {
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

        /** Returns true if any base in the range has coverage of > 0 */
        public boolean hasCoverage() {
            // NB: if this is expensive, we could easily pre-compute this as we go along in addBase
            for (final int s : depths) {
                if (s > 0) return true;
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

    /** Tracks the number of read pairs that we see that are PF (used to calculate library size) */
    public long PF_SELECTED_PAIRS;

    /** Tracks the number of unique PF reads pairs we see (used to calc library size) */
    public long PF_SELECTED_UNIQUE_PAIRS;

    /** The number of PF unique reads that are aligned with mapping score > 0 to the reference genome. */
    public long PF_UQ_READS_ALIGNED;

    /** The number of PF unique bases that are aligned with mapping score > 0 to the reference genome. */
    public long PF_BASES_ALIGNED;

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

    /** The mean coverage of targets. */
    public double MEAN_TARGET_COVERAGE;

    /** The median coverage of targets. */
    public double MEDIAN_TARGET_COVERAGE;

    /** The fraction of targets that did not reach coverage=1 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
    public double PCT_EXC_DUPE;

    /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality. */
    public double PCT_EXC_MAPQ;

    /** The fraction of aligned bases that were filtered out because they were of low base quality. */
    public double PCT_EXC_BASEQ;

    /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
    public double PCT_EXC_OVERLAP;

    /** The fraction of aligned bases that were filtered out because they did not align over a target base. */
    public double PCT_EXC_OFF_TARGET;

    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;

    /** The percentage of all target bases achieving 1X or greater coverage. */
    public double PCT_TARGET_BASES_1X;
    /** The percentage of all target bases achieving 2X or greater coverage. */
    public double PCT_TARGET_BASES_2X;
    /** The percentage of all target bases achieving 10X or greater coverage. */
    public double PCT_TARGET_BASES_10X;
    /** The percentage of all target bases achieving 20X or greater coverage. */
    public double PCT_TARGET_BASES_20X;
    /** The percentage of all target bases achieving 30X or greater coverage. */
    public double PCT_TARGET_BASES_30X;
    /** The percentage of all target bases achieving 40X or greater coverage. */
    public double PCT_TARGET_BASES_40X;
    /** The percentage of all target bases achieving 50X or greater coverage. */
    public double PCT_TARGET_BASES_50X;
    /** The percentage of all target bases achieving 100X or greater coverage. */
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
