/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.TheoreticalSensitivity;
import picard.filter.CountingAdapterFilter;
import picard.filter.CountingMapQFilter;
import picard.metrics.MultilevelMetrics;
import picard.metrics.PerUnitMetricCollector;
import picard.metrics.SAMRecordMultiLevelCollector;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.*;
import java.util.stream.LongStream;

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

    private final File perTargetCoverage;  // If not null, per-target coverage summaries are written to this file
    private final File perBaseCoverage;    // If not null, per-base(!) coverage summaries are written to this file

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

    private final int coverageCap;

    private final int sampleSize;

    // histogram of depths. does not include bases with quality less than MINIMUM_BASE_QUALITY (default 20)
    // give it the bin label "coverage_or_base_quality" to make clear that in the metrics file the coverage and base quality histograms share the same bin column on the left
    private final Histogram<Integer> highQualityDepthHistogram = new Histogram<>("coverage_or_base_quality", "high_quality_coverage_count");

    // histogram of base qualities. includes all but quality 2 bases. we use this histogram to calculate theoretical het sensitivity.
    private final Histogram<Integer> unfilteredDepthHistogram = new Histogram<>("coverage_or_base_quality", "unfiltered_coverage_count");

    // histogram of base qualities. includes all but quality 2 bases. we use this histogram to calculate theoretical het sensitivity.
    private final Histogram<Integer> unfilteredBaseQHistogram = new Histogram<>("baseq", "unfiltered_baseq_count");

    private static final double LOG_ODDS_THRESHOLD = 3.0;

    private final File theoreticalSensitivityOutput = null;

    private final int minimumMappingQuality;
    private final int minimumBaseQuality;
    private final boolean clipOverlappingReads;
    private boolean noSideEffects;
    private final boolean includeIndels;

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
                                  final File perBaseCoverage,
                                  final IntervalList targetIntervals,
                                  final IntervalList probeIntervals,
                                  final String probeSetName,
                                  final int nearProbeDistance,
                                  final int minimumMappingQuality,
                                  final int minimumBaseQuality,
                                  final boolean clipOverlappingReads,
                                  final int coverageCap,
                                  final int sampleSize) {
        this(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance, minimumMappingQuality, minimumBaseQuality, clipOverlappingReads, false, false, coverageCap, sampleSize);
    }

    public TargetMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                  final List<SAMReadGroupRecord> samRgRecords,
                                  final ReferenceSequenceFile refFile,
                                  final File perTargetCoverage,
                                  final File perBaseCoverage,
                                  final IntervalList targetIntervals,
                                  final IntervalList probeIntervals,
                                  final String probeSetName,
                                  final int nearProbeDistance,
                                  final int minimumMappingQuality,
                                  final int minimumBaseQuality,
                                  final boolean clipOverlappingReads,
                                  final boolean noSideEffects,
                                  final boolean includeIndels,
                                  final int coverageCap,
                                  final int sampleSize) {
        this.perTargetCoverage = perTargetCoverage;
        this.perBaseCoverage   = perBaseCoverage;
        this.probeSetName = probeSetName;
        this.nearProbeDistance = nearProbeDistance;

        this.allProbes  = probeIntervals;
        this.allTargets = targetIntervals;
        this.coverageCap = coverageCap;
        this.sampleSize = sampleSize;

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
                final ReferenceSequence rs = refFile.getSubsequenceAt(target.getContig(), target.getStart(), target.getEnd());
                intervalToGc.put(target,SequenceUtil.calculateGc(rs.getBases()));
            }
        }

        this.minimumMappingQuality = minimumMappingQuality;
        this.minimumBaseQuality = minimumBaseQuality;
        this.clipOverlappingReads = clipOverlappingReads;
        this.noSideEffects = noSideEffects;
        this.includeIndels = includeIndels;

        setup(accumulationLevels, samRgRecords);
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
        final PerUnitTargetMetricCollector collector =  new PerUnitTargetMetricCollector(probeSetName, coverageByTargetForRead.keySet(),
                sample, library, readGroup, probeTerritory, targetTerritory, genomeSize,
                intervalToGc, minimumMappingQuality, minimumBaseQuality, clipOverlappingReads, includeIndels);
        if (this.probeSetName != null) {
            collector.setBaitSetName(probeSetName);
        }

        return collector;
    }

    @Override
    protected PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> makeAllReadCollector() {
        final PerUnitTargetMetricCollector collector = (PerUnitTargetMetricCollector) makeChildCollector(null, null, null);
        if (perTargetCoverage != null) collector.setPerTargetOutput(perTargetCoverage);
        if (perBaseCoverage   != null) collector.setPerBaseOutput(perBaseCoverage);

        return collector;
    }

    /**
     * Collect the Target Metrics for one unit of "accumulation" (i.e. for one sample, or for one library ...)
     */
    public class PerUnitTargetMetricCollector implements PerUnitMetricCollector<METRIC_TYPE, Integer, SAMRecord> {
        private final Map<Interval, Double> intervalToGc;
        private File perTargetOutput;
        private File perBaseOutput;

        final long[] baseQHistogramArray = new long[Byte.MAX_VALUE];
        // A Map to accumulate per-bait-region (i.e. merge of overlapping targets) coverage
        // excludes bases with qualities lower than minimumBaseQuality (default 20)
        private final Map<Interval, Coverage> highQualityCoverageByTarget;

        // only excludes bases with quality 2. collected for theoretical set sensitivity
        private final Map<Interval, Coverage> unfilteredCoverageByTarget;

        private final TargetMetrics metrics = new TargetMetrics();
        private final int minimumBaseQuality;
        private final CountingAdapterFilter adapterFilter;
        private final CountingMapQFilter mapQFilter;
        private final boolean clipOverlappingReads;
        private final boolean includeIndels;

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
                                            final boolean clipOverlappingReads,
                                            final boolean includeIndels) {
            this.metrics.SAMPLE      = sample;
            this.metrics.LIBRARY     = library;
            this.metrics.READ_GROUP  = readGroup;
            this.metrics.PROBE_SET   = probeSetName;

            metrics.PROBE_TERRITORY  = probeTerritory;
            metrics.TARGET_TERRITORY = targetTerritory;
            metrics.GENOME_SIZE      = genomeSize;

            highQualityCoverageByTarget = new LinkedHashMap<>(coverageTargets.size() * 2, 0.5f);
            unfilteredCoverageByTarget =  new LinkedHashMap<>(coverageTargets.size() * 2, 0.5f);

            for (final Interval target : coverageTargets) {
                highQualityCoverageByTarget.put(target, new Coverage(target, 0));
                unfilteredCoverageByTarget.put(target, new Coverage(target, 0));
            }

            this.mapQFilter = new CountingMapQFilter(minimumMappingQuality);
            this.adapterFilter = new CountingAdapterFilter();
            this.minimumBaseQuality = minimumBaseQuality;
            this.intervalToGc = intervalToGc;
            this.clipOverlappingReads = clipOverlappingReads;
            this.includeIndels = includeIndels;
        }

        /** Sets the (optional) File to write per-target coverage information to. If null (the default), no file is produced. */
        public void setPerTargetOutput(final File perTargetOutput) {
            this.perTargetOutput = perTargetOutput;
        }

        /** Sets the (optional) File to write per-base coverage information to. If null (the default), no file is produced. */
        public void setPerBaseOutput(final File perBaseOutput) {
            this.perBaseOutput = perBaseOutput;
        }

        /** Sets the name of the bait set explicitly instead of inferring it from the bait file. */
        public void setBaitSetName(final String name) {
            this.metrics.PROBE_SET = name;
        }

        /**
         * Returns the accumulated coverage per target.  Note that while the returned Map is
         * immutable, it is possible that the underlying Map will continue to be mutated if
         * the map is retrieved prior to additional calls to {@link #acceptRecord(SAMRecord)}.
         */
        public Map<Interval, Coverage> getCoverageByTarget() {
            return Collections.unmodifiableMap(this.highQualityCoverageByTarget);
        }

        /** Adds information about an individual SAMRecord to the statistics. */
        public void acceptRecord(final SAMRecord record) {
            // Just ignore secondary alignments altogether
            if (record.isSecondaryAlignment()) return;

            // Cache some things, and compute the total number of bases aligned in the record.
            final boolean mappedInPair = record.getReadPairedFlag() && !record.getReadUnmappedFlag() && !record.getMateUnmappedFlag() && !record.getSupplementaryAlignmentFlag();
            final byte[] baseQualities = record.getBaseQualities();
            int basesAlignedInRecord = 0;
            if (!record.getReadUnmappedFlag()) {
                for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                    basesAlignedInRecord += block.getLength();
                }
            }

            /* Count metrics related to # of base and reads. Consider supplemental alignments for base counting but not record counting.
               Do this counting *prior* to applying filters to ensure we match other metrics' computation for these values, such as AlignmentSummaryMetrics. */

            // READ Based Metrics
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

            ///////////////////////////////////////////////////////////////////
            // Non-PF reads can be totally ignored beyond this point
            ///////////////////////////////////////////////////////////////////
            if (record.getReadFailsVendorQualityCheckFlag()) return;

            // BASE Based Metrics
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

            ///////////////////////////////////////////////////////////////////
            // Unmapped reads can be totally ignored beyond this point
            ///////////////////////////////////////////////////////////////////
            if (record.getReadUnmappedFlag()) return;

            // Prefetch the list of target and bait overlaps here as they're needed multiple times.
            final Interval read = new Interval(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd());
            final Collection<Interval> targets = targetDetector.getOverlaps(read);
            final Collection<Interval> probes  = probeDetector.getOverlaps(read);

            // Calculate the values we need for HS_LIBRARY_SIZE
            if (!record.getSupplementaryAlignmentFlag() &&
                    record.getReadPairedFlag() &&
                    record.getFirstOfPairFlag() &&
                    !record.getReadUnmappedFlag() &&
                    !record.getMateUnmappedFlag() &&
                    !probes.isEmpty()) {
                ++this.metrics.PF_SELECTED_PAIRS;
                if (!record.getDuplicateReadFlag()) ++this.metrics.PF_SELECTED_UNIQUE_PAIRS;
            }

            // Compute the bait-related metrics *before* applying the duplicate read
            // filtering, overlap clipping and the map-q threshold, since those would
            // skew the assay-related metrics
            {
                int onBaitBases = 0;

                if (!probes.isEmpty()) {
                    for (final Interval bait : probes) {
                        for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                            final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());

                            for (int pos = block.getReferenceStart(); pos <= end; ++pos) {
                                if (pos >= bait.getStart() && pos <= bait.getEnd()) ++onBaitBases;
                            }
                        }
                    }

                    this.metrics.ON_PROBE_BASES += onBaitBases;
                    this.metrics.NEAR_PROBE_BASES += (basesAlignedInRecord - onBaitBases);
                } else {
                    this.metrics.OFF_PROBE_BASES += basesAlignedInRecord;
                }
            }

            ///////////////////////////////////////////////////////////////////
            // Duplicate reads can be totally ignored beyond this point
            ///////////////////////////////////////////////////////////////////
            if (record.getDuplicateReadFlag()) {
                this.metrics.PCT_EXC_DUPE += basesAlignedInRecord;
                return;
            }

            ///////////////////////////////////////////////////////////////////
            // MapQ 0 adapter reads can be ignored beyond this point
            // but first, make sure we count the (aligned) adapters.
            ///////////////////////////////////////////////////////////////////
            if (this.adapterFilter.filterOut(record)) {
                return;
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

                // If clipping resulted in the read becoming unmapped (because all bases were clipped), return here
                if (rec.getReadUnmappedFlag()) {
                    return;
                }
            } else {
                rec = record;
            }

            // Calculate all the things that require examining individual bases in the read. This includes:
            //   1. Per-base coverage
            //   2. The number of reads contributing to per-base coverage per target
            //   3. Unfiltered coverage information for het sensitivity
            //   4. The count of bases rejected for being low baseq or off-target
            //   5. The count of overall on-target bases, and on-target bases from paired reads
            final Set<Interval> coveredTargets = new HashSet<>(); // Each target is added to this the first time the read covers it
            int readOffset = 0;
            int refOffset  = rec.getAlignmentStart() - 1;

            for (final CigarElement cig : rec.getCigar()) {
                final CigarOperator op = cig.getOperator();
                final int len = cig.getLength();

                for (int i=0; i<len; ++i) {
                    if (op.isAlignment() || (this.includeIndels && op.isIndel())) {
                        final int refPos       = refOffset + 1;
                        final int qual         = baseQualities[readOffset];
                        final boolean highQual = qual >= this.minimumBaseQuality;
                        final boolean onTarget = overlapsAny(refPos, targets);
                        final boolean incrementPerTargetCoverage = op != CigarOperator.INSERTION;  // Inserted bases don't have a target position

                        // Firstly handle all the summary metrics
                        if (!highQual) {
                            metrics.PCT_EXC_BASEQ++;
                        } else if (!onTarget) {
                            metrics.PCT_EXC_OFF_TARGET++;
                        } else {
                            metrics.ON_TARGET_BASES++;
                            if (mappedInPair) metrics.ON_TARGET_FROM_PAIR_BASES++;
                        }

                        // Then go through the per-target/per-base hq and unfiltered coverage
                        // The cutoff of > 2 is because even the unfilteredCoverage doesn't want those bases
                        if (qual > 2 && incrementPerTargetCoverage && onTarget) {
                            for (final Interval target : targets) {
                                if (overlapsInterval(refPos, target)) {
                                    final int targetOffset = refPos - target.getStart();

                                    // Unfiltered first (for theoretical het sensitivity)
                                    final Coverage ufCoverage = unfilteredCoverageByTarget.get(target);
                                    ufCoverage.addBase(targetOffset);
                                    if (ufCoverage.getDepths()[targetOffset] <= coverageCap) baseQHistogramArray[qual]++;

                                    // Then filtered
                                    if (highQual) {
                                        final Coverage hqCoverage = highQualityCoverageByTarget.get(target);
                                        hqCoverage.addBase(targetOffset);

                                        if (coveredTargets.add(target)) {
                                            hqCoverage.incrementReadCount();
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Finally update the offsets!
                    if (op.consumesReadBases()) readOffset += 1;
                    if (op.consumesReferenceBases()) refOffset += 1;
                }
            }
        }

        /* Returns true if the `pos` is between the `start` and `end` of at least one interval. */
        private boolean overlapsAny(final int pos, final Collection<Interval> intervals) {
            for (final Interval interval : intervals) {
                if (overlapsInterval(pos, interval)) return true;
            }
            return false;
        }

        /** Returns true if the position is within the start-end range inclusive of the given interval. */
        private boolean overlapsInterval(final int pos, final Interval interval) {
            return pos >= interval.getStart() && pos <= interval.getEnd();
        }

        @Override
        public void finish() {
            metrics.PCT_PF_READS            = metrics.PF_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS         = metrics.PF_UNIQUE_READS / (double) metrics.TOTAL_READS;
            metrics.PCT_PF_UQ_READS_ALIGNED = metrics.PF_UQ_READS_ALIGNED / (double) metrics.PF_UNIQUE_READS;

            final double denominator        = metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES + metrics.OFF_PROBE_BASES;

            metrics.PCT_SELECTED_BASES      = (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES) / denominator;
            metrics.PCT_OFF_PROBE           = metrics.OFF_PROBE_BASES / denominator;
            metrics.ON_PROBE_VS_SELECTED    = metrics.ON_PROBE_BASES / (double) (metrics.ON_PROBE_BASES + metrics.NEAR_PROBE_BASES);
            metrics.MEAN_PROBE_COVERAGE     = metrics.ON_PROBE_BASES / (double) metrics.PROBE_TERRITORY;
            metrics.FOLD_ENRICHMENT         = (metrics.ON_PROBE_BASES/ denominator) / ((double) metrics.PROBE_TERRITORY / metrics.GENOME_SIZE);

            metrics.PCT_EXC_DUPE           /= (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_ADAPTER         = adapterFilter.getFilteredBases() / (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_MAPQ            = mapQFilter.getFilteredBases() / (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_BASEQ          /= (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_OVERLAP        /= (double) metrics.PF_BASES_ALIGNED;
            metrics.PCT_EXC_OFF_TARGET     /= (double) metrics.PF_BASES_ALIGNED;

            calculateTargetCoverageMetrics();
            calculateTheoreticalHetSensitivity();
            calculateGcMetrics();
            emitPerBaseCoverageIfRequested();
        }

        /** Calculates how much additional sequencing is needed to raise 80% of bases to the mean for the lane. */
        private void calculateTargetCoverageMetrics() {
            final long[] highQualityCoverageHistogramArray = new long[coverageCap+1];
            int zeroCoverageTargets = 0;

            // the number of bases we counted towards the depth histogram plus those that got thrown out by the coverage cap
            long totalCoverage = 0;

            // the maximum depth at any target base
            long maxDepth = 0;

            // the minimum depth at any target base
            long minDepth = Long.MAX_VALUE;

            // The "how many target bases at at-least X" calculations.
            // downstream code relies on this array being sorted in ascending order
            final int[] targetBasesDepth = {0, 1, 2, 10, 20, 30, 40, 50, 100};

            // counts for how many target bases are at at least X coverage,
            // where X corresponds to the value at the same offset in targetBasesDepth
            final int[] targetBases = new int[targetBasesDepth.length];

            // for each target, count up the depth for each base and increment the depth histogram array
            for (final Coverage c : this.highQualityCoverageByTarget.values()) {
                if (!c.hasCoverage()) {
                    zeroCoverageTargets++;
                    highQualityCoverageHistogramArray[0] += c.interval.length();
                    targetBases[0] += c.interval.length();
                    minDepth = 0;
                    continue;
                }

                for (final int depth : c.getDepths()) {
                    totalCoverage += depth;
                    highQualityCoverageHistogramArray[Math.min(depth, coverageCap)]++;
                    maxDepth = Math.max(maxDepth, depth);
                    minDepth = Math.min(minDepth, depth);

                    // Add to the "how many target bases at at-least X" calculations.
                    for (int i = 0; i < targetBasesDepth.length; i++) {
                        if (depth >= targetBasesDepth[i]) targetBases[i]++;
                        else break; // NB: assumes that targetBasesDepth is sorted in ascending order
                    }
                }
            }

            if (targetBases[0] !=  highQualityCoverageByTarget.keySet().stream().mapToInt(Interval::length).sum()) {
                throw new PicardException("the number of target bases with at least 0x coverage does not equal the number of target bases");
            }

            for (int i = 0; i < highQualityCoverageHistogramArray.length; ++i) {
                highQualityDepthHistogram.increment(i, highQualityCoverageHistogramArray[i]);
            }

            // we do this instead of highQualityDepthHistogram.getMean() because the histogram imposes a coverage cap
            metrics.MEAN_TARGET_COVERAGE = (double) totalCoverage / metrics.TARGET_TERRITORY;
            metrics.MEDIAN_TARGET_COVERAGE = highQualityDepthHistogram.getMedian();
            metrics.MAX_TARGET_COVERAGE = maxDepth;
            // Use Math.min() to account for edge case where highQualityCoverageByTarget is empty (minDepth=Long.MAX_VALUE)
            metrics.MIN_TARGET_COVERAGE = Math.min(minDepth, maxDepth);

            // compute the coverage value such that 80% of target bases have better coverage than it i.e. 20th percentile
            // this roughly measures how much we must sequence extra such that 80% of target bases have coverage at least as deep as the current mean coverage
            metrics.FOLD_80_BASE_PENALTY = metrics.MEAN_TARGET_COVERAGE / highQualityDepthHistogram.getPercentile(0.2);
            metrics.ZERO_CVG_TARGETS_PCT = zeroCoverageTargets / (double) allTargets.getIntervals().size();

            // Store the "how many bases at at-least X" calculations.
            metrics.PCT_TARGET_BASES_1X   = (double) targetBases[1] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_2X   = (double) targetBases[2] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_10X  = (double) targetBases[3] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_20X  = (double) targetBases[4] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_30X  = (double) targetBases[5] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_40X  = (double) targetBases[6] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_50X  = (double) targetBases[7] / (double) targetBases[0];
            metrics.PCT_TARGET_BASES_100X = (double) targetBases[8] / (double) targetBases[0];
        }

        private void calculateTheoreticalHetSensitivity(){
            final long[] unfilteredDepthHistogramArray = new long[coverageCap + 1];

            // collect the unfiltered coverages (i.e. only quality 2 bases excluded) for all targets into a histogram array
            for (final Coverage c : this.unfilteredCoverageByTarget.values()) {
                if (!c.hasCoverage()) {
                    unfilteredDepthHistogramArray[0] += c.interval.length();
                    continue;
                }

                for (final int depth : c.getDepths()) {
                    unfilteredDepthHistogramArray[Math.min(depth, coverageCap)]++;
                }
            }

            if (LongStream.of(baseQHistogramArray).sum() != LongStream.rangeClosed(0, coverageCap).map(i -> i * unfilteredDepthHistogramArray[(int)i]).sum()) {
                throw new PicardException("numbers of bases in the base quality histogram and the coverage histogram are not equal");
            }

            // TODO: normalize the arrays directly. then we don't have to convert to Histograms
            for (int i=0; i<baseQHistogramArray.length; ++i) {
                unfilteredBaseQHistogram.increment(i, baseQHistogramArray[i]);
            }

            for (int i = 0; i < unfilteredDepthHistogramArray.length; i++){
                unfilteredDepthHistogram.increment(i, unfilteredDepthHistogramArray[i]);
            }

            final double [] depthDoubleArray = TheoreticalSensitivity.normalizeHistogram(unfilteredDepthHistogram);
            final double [] baseQDoubleArray = TheoreticalSensitivity.normalizeHistogram(unfilteredBaseQHistogram);
            metrics.HET_SNP_SENSITIVITY = TheoreticalSensitivity.hetSNPSensitivity(depthDoubleArray, baseQDoubleArray, sampleSize, LOG_ODDS_THRESHOLD);
            metrics.HET_SNP_Q = QualityUtil.getPhredScoreFromErrorProbability((1 - metrics.HET_SNP_SENSITIVITY));

        }

        /** Emits a file with per base coverage if an output file has been set. */
        private void emitPerBaseCoverageIfRequested() {
            if (this.perBaseOutput == null) return;

            final PrintWriter out = new PrintWriter(IOUtil.openFileForBufferedWriting(this.perBaseOutput));
            out.println("chrom\tpos\ttarget\tcoverage");
            for (final Map.Entry<Interval,Coverage> entry : this.highQualityCoverageByTarget.entrySet()) {
                final Interval interval = entry.getKey();
                final String chrom = interval.getContig();
                final int firstBase = interval.getStart();

                final int[] cov = entry.getValue().getDepths();
                for (int i = 0; i < cov.length; ++i) {
                    out.print(chrom);
                    out.print('\t');
                    out.print(firstBase + i);
                    out.print('\t');
                    out.print(interval.getName());
                    out.print('\t');
                    out.print(cov[i]);
                    out.println();
                }
            }

            out.close();
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
                        out.println("chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage\tnormalized_coverage\tmin_normalized_coverage\tmax_normalized_coverage\tmin_coverage\tmax_coverage\tpct_0x\tread_count");
                    }
                    else {
                        out = null;
                    }
                }
                catch (final IOException ioe) { throw new RuntimeIOException(ioe); }

                final int bins = 101;
                final long[] targetBasesByGc  = new long[bins];
                final long[] alignedBasesByGc = new long[bins];

                for (final Map.Entry<Interval,Coverage> entry : this.highQualityCoverageByTarget.entrySet()) {
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

                        out.println(interval.getContig() + "\t" +
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
                                fmt.format(targetBasesAt0x / interval.length()) + "\t" +
                                fmt.format(cov.readCount)
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
            hsMetricsComparableMetricsFile.addHistogram(highQualityDepthHistogram);
            hsMetricsComparableMetricsFile.addHistogram(unfilteredBaseQHistogram);
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
        public long readCount = 0;

        /** Constructs a new coverage object for the provided mapping with the desired padding either side. */
        public Coverage(final Interval i, final int padding) {
            this.interval = i;
            this.depths = new int[interval.length() + 2*padding];
        }

        /** Adds a single point of depth at the desired offset into the coverage array. */
        public void addBase(final int offset) {
            addBase(offset, 1);
        }

        /** Adds some depth at the desired offset into the coverage array. */
        public void addBase(final int offset, final int depth) {
            if (offset >= 0 && offset < this.depths.length && this.depths[offset] < Integer.MAX_VALUE - depth) {
                this.depths[offset] += depth;
            }
        }

        /** Increments the # of reads mapping to this target. */
        public void incrementReadCount() {
            this.readCount++;
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

        public long getTotal() {
            long total = 0;
            for (int i=0; i<depths.length; ++i) {
                total += (total < Long.MAX_VALUE - depths[i]) ? depths[i] : Long.MAX_VALUE - total;
            }
            return total;
        }

        @Override
        public String toString() {
            return "TargetedMetricCollector(interval=" + interval + ", depths = [" + StringUtil.intValuesToString(this.depths) + "])";
        }
    }

    public Histogram<Integer> getBaseQualityHistogram() {
        return unfilteredBaseQHistogram;
    }

    public Histogram<Integer> getDepthHistogram() {
        return unfilteredDepthHistogram;
    }
}
