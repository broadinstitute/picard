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
package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.filter.CountingDuplicateFilter;
import picard.filter.CountingFilter;
import picard.filter.CountingMapQFilter;
import picard.filter.CountingPairedFilter;
import picard.util.MathUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 *
 * @author tfennell
 */
@CommandLineProgramProperties(
        usage = CollectWgsMetrics.USAGE_SUMMARY + CollectWgsMetrics.USAGE_DETAILS,
        usageShort = CollectWgsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectWgsMetrics extends CommandLineProgram {
static final String USAGE_SUMMARY = "Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.";
static final String USAGE_DETAILS = "<p>This tool collects metrics about the fractions of reads that pass base- and mapping-quality "+
"filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum "+
"read depths (coverage cap) are user defined.</p>" +

"<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>"+
"<h4>Usage Example:</h4>"+
"<pre>"  +
"java -jar picard.jar CollectWgsMetrics \\<br /> " +
"      I=input.bam \\<br /> "+
"      O=collect_wgs_metrics.txt \\<br /> " +
"      R=reference_sequence.fasta " +
"</pre>" +
"Please see "+
"<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics'>CollectWgsMetrics</a> "+
"for detailed explanations of the output metrics." +
"<hr />"
;

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "The reference sequence fasta aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage. N bases will be treated as having a base quality " +
            "of negative infinity and will therefore be excluded from coverage regardless of the value of this parameter.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(shortName = "CAP", doc = "Treat positions with coverage exceeding this value as if they had coverage at this value (but calculate the difference for PCT_EXC_CAPPED).", overridable = true)
    public int COVERAGE_CAP = 250;

    @Option(doc="At positions with coverage exceeding this value, completely ignore reads that accumulate beyond this value (so that they will not be considered for PCT_EXC_CAPPED).  Used to keep memory consumption in check, but could create bias if set too low", overridable = true)
    public int LOCUS_ACCUMULATION_CAP = 100000;

    @Option(doc = "For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;

    @Option(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    @Option(doc="If true, count unpaired reads, and paired reads with one end unmapped")
    public boolean COUNT_UNPAIRED = false;

    @Option(doc="Sample Size used for Theoretical Het Sensitivity sampling. Default is 10000.", optional = true)
    public int SAMPLE_SIZE=10000;

    @Option(doc = "An interval list file that contains the positions to restrict the assessment. Please note that " +
            "all bases of reads that overlap these intervals will be considered, even if some of those bases extend beyond the boundaries of " +
            "the interval. The ideal use case for this argument is to use it to restrict the calculation to a subset of (whole) contigs. To " +
            "restrict the calculation just to individual positions without overlap, please see CollectWgsMetricsFromSampledSites.",
            optional = true, overridable = true)
    public File INTERVALS = null;

    private SAMFileHeader header = null;

    private final Log log = Log.getInstance(CollectWgsMetrics.class);
    private static final double LOG_ODDS_THRESHOLD = 3.0;

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetrics extends MergeableMetricBase {

        /** The intervals over which this metric was computed. */
        @MergingIsManual
        protected IntervalList intervals;

        /** The count of sites with a given observed depth. */
        @MergingIsManual
        protected final Histogram<Integer> depthHistogram;

        /** The count of bases observed with a given base quality. */
        @MergingIsManual
        protected final Histogram<Integer> baseQHistogram;

        /** The maximum depth/coverage to consider. */
        @MergeByAssertEquals
        protected final int coverageCap;

        /** The sample size used for theoretical het sensitivity. */
        @NoMergingKeepsValue
        protected final int theoreticalHetSensitivitySampleSize;

        /**
         * Create an instance of this metric that is not mergeable.
         */
        public WgsMetrics() {
            intervals                           = null;
            depthHistogram                      = null;
            baseQHistogram                      = null;
            theoreticalHetSensitivitySampleSize = -1;
            coverageCap                         = -1;
        }

        /**
         * Create an instance of this metric that is mergeable.
         *
         * @param depthHistogram the count of genomic positions observed for each observed depth.
         * @param pctExcludedByMapq the fraction of aligned bases that were filtered out because they were in reads with low mapping quality.
         * @param pctExcludedByDupes the fraction of aligned bases that were filtered out because they were in reads marked as duplicates.
         * @param pctExcludedByPairing the fraction of bases that were filtered out because they were in reads without a mapped mate pair.
         * @param pctExcludedByBaseq the fraction of aligned bases that were filtered out because they were of low base quality.
         * @param pctExcludedByOverlap the fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
         * @param pctExcludedByCapping the fraction of aligned bases that were filtered out because they would have raised coverage above the capped value.
         * @param pctExcludeTotal the fraction of bases excluded across all filters.
         * @param coverageCap Treat positions with coverage exceeding this value as if they had coverage at this value.
         * @param baseQHistogram the count of bases observed with a given quality.
         * @param theoreticalHetSensitivitySampleSize the sample size used for theoretical het sensitivity sampling.
         */
        public WgsMetrics(final IntervalList intervals,
                          final Histogram<Integer> depthHistogram,
                          final double pctExcludedByMapq,
                          final double pctExcludedByDupes,
                          final double pctExcludedByPairing,
                          final double pctExcludedByBaseq,
                          final double pctExcludedByOverlap,
                          final double pctExcludedByCapping,
                          final double pctExcludeTotal,
                          final int coverageCap,
                          final Histogram<Integer> baseQHistogram,
                          final int theoreticalHetSensitivitySampleSize) {
            this.intervals      = intervals.uniqued();
            this.depthHistogram = depthHistogram;
            this.baseQHistogram = baseQHistogram;
            this.coverageCap    = coverageCap;
            this.theoreticalHetSensitivitySampleSize = theoreticalHetSensitivitySampleSize;

            PCT_EXC_MAPQ     = pctExcludedByMapq;
            PCT_EXC_DUPE     = pctExcludedByDupes;
            PCT_EXC_UNPAIRED = pctExcludedByPairing;
            PCT_EXC_BASEQ    = pctExcludedByBaseq;
            PCT_EXC_OVERLAP  = pctExcludedByOverlap;
            PCT_EXC_CAPPED   = pctExcludedByCapping;
            PCT_EXC_TOTAL    = pctExcludeTotal;

            calculateDerivedFields();
        }

        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        @NoMergingIsDerived
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        @NoMergingIsDerived
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        @NoMergingIsDerived
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        @NoMergingIsDerived
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        @NoMergingIsDerived
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        @NoMergingIsDerived
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        @NoMergingIsDerived
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        @NoMergingIsDerived
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        @NoMergingIsDerived
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        @NoMergingIsDerived
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        @NoMergingIsDerived
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        @NoMergingIsDerived
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 1X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_1X;
        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
         @NoMergingIsDerived
         public double PCT_100X;

        /** The theoretical HET SNP sensitivity. */
        @NoMergingIsDerived
        public double HET_SNP_SENSITIVITY;

        /** The Phred Scaled Q Score of the theoretical HET SNP sensitivity. */
        @NoMergingIsDerived
        public double HET_SNP_Q;

        /**
         * Merges the various PCT_EXC_* metrics.
         * @param other metric to merge into this one.
         */
        @Override
        public void merge(final MergeableMetricBase other) {
            final WgsMetrics otherMetric = (WgsMetrics) other;

            if (depthHistogram == null || otherMetric.depthHistogram == null) {
                throw new PicardException("Depth histogram is required when deriving metrics.");
            }

            // Union the intervals over which bases are called.  They should have no overlaps!
            // NB: interval lists are already uniqued.
            final long genomeTerritory = this.intervals.getBaseCount() + otherMetric.intervals.getBaseCount();
            this.intervals.addall(otherMetric.intervals.getIntervals());
            this.intervals = this.intervals.uniqued();
            if (this.intervals.getBaseCount() != genomeTerritory) {
                throw new PicardException("Trying to merge WgsMetrics calculated on intervals that overlap.");
            }

            // NB:
            // Since: PCT_EXC_TOTAL     = (totalWithExcludes - thisMetricTotal) / totalWithExcludes;
            // Thus:  totalWithExcludes = total / (1 - PCT_EXC_TOTAL)
            // Proof: Exercise is left to the reader.
            final long thisMetricTotal        = (long) depthHistogram.getSum();
            final long otherMetricTotal       = (long) otherMetric.depthHistogram.getSum();
            final long total                  = thisMetricTotal + otherMetricTotal;
            final long thisTotalWithExcludes  = (long) (thisMetricTotal / (1.0 - PCT_EXC_TOTAL));
            final long otherTotalWithExcludes = (long) (otherMetricTotal / (1.0 - otherMetric.PCT_EXC_TOTAL));
            final double totalWithExcludes    = thisTotalWithExcludes + otherTotalWithExcludes;

            if (0 < totalWithExcludes) {
                PCT_EXC_DUPE     = (PCT_EXC_DUPE * thisTotalWithExcludes + otherMetric.PCT_EXC_DUPE * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_MAPQ     = (PCT_EXC_MAPQ * thisTotalWithExcludes + otherMetric.PCT_EXC_MAPQ * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_UNPAIRED = (PCT_EXC_UNPAIRED * thisTotalWithExcludes + otherMetric.PCT_EXC_UNPAIRED * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_BASEQ    = (PCT_EXC_BASEQ * thisTotalWithExcludes + otherMetric.PCT_EXC_BASEQ * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_OVERLAP  = (PCT_EXC_OVERLAP * thisTotalWithExcludes + otherMetric.PCT_EXC_OVERLAP * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_CAPPED   = (PCT_EXC_CAPPED * thisTotalWithExcludes + otherMetric.PCT_EXC_CAPPED * otherTotalWithExcludes) / totalWithExcludes;
                PCT_EXC_TOTAL    = (totalWithExcludes - total) / totalWithExcludes;
            }

            // do any merging that are dictated by the annotations.
            super.merge(other);

            // merge the histograms
            this.depthHistogram.addHistogram(otherMetric.depthHistogram);
            if (baseQHistogram != null && otherMetric.baseQHistogram != null) this.baseQHistogram.addHistogram(otherMetric.baseQHistogram);
        }

        @Override
        public void calculateDerivedFields() {
            if (depthHistogram == null) throw new PicardException("Depth histogram is required when deriving metrics.");
            if (baseQHistogram != null && theoreticalHetSensitivitySampleSize <= 0) {
                throw new PicardException("Sample size is required when a baseQ histogram is given when deriving metrics.");
            }

            final long[] depthHistogramArray = new long[coverageCap+1];
            for (final Histogram.Bin<Integer> bin : depthHistogram.values()) {
                final int depth = Math.min((int) bin.getIdValue(), coverageCap);
                depthHistogramArray[depth] += bin.getValue();
            }

            GENOME_TERRITORY = (long) depthHistogram.getSumOfValues();
            MEAN_COVERAGE    = depthHistogram.getMean();
            SD_COVERAGE      = depthHistogram.getStandardDeviation();
            MEDIAN_COVERAGE  = depthHistogram.getMedian();
            MAD_COVERAGE     = depthHistogram.getMedianAbsoluteDeviation();

            PCT_1X    = MathUtil.sum(depthHistogramArray, 1, depthHistogramArray.length)   / (double) GENOME_TERRITORY;
            PCT_5X    = MathUtil.sum(depthHistogramArray, 5, depthHistogramArray.length)   / (double) GENOME_TERRITORY;
            PCT_10X   = MathUtil.sum(depthHistogramArray, 10, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_15X   = MathUtil.sum(depthHistogramArray, 15, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_20X   = MathUtil.sum(depthHistogramArray, 20, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_25X   = MathUtil.sum(depthHistogramArray, 25, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_30X   = MathUtil.sum(depthHistogramArray, 30, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_40X   = MathUtil.sum(depthHistogramArray, 40, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_50X   = MathUtil.sum(depthHistogramArray, 50, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_60X   = MathUtil.sum(depthHistogramArray, 60, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_70X   = MathUtil.sum(depthHistogramArray, 70, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_80X   = MathUtil.sum(depthHistogramArray, 80, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_90X   = MathUtil.sum(depthHistogramArray, 90, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
            PCT_100X  = MathUtil.sum(depthHistogramArray, 100, depthHistogramArray.length) / (double) GENOME_TERRITORY;

            // Get Theoretical Het SNP Sensitivity
            if (baseQHistogram != null) {
                final double[] depthDoubleArray = TheoreticalSensitivity.normalizeHistogram(depthHistogram);
                final double[] baseQDoubleArray = TheoreticalSensitivity.normalizeHistogram(baseQHistogram);
                HET_SNP_SENSITIVITY = TheoreticalSensitivity.hetSNPSensitivity(depthDoubleArray, baseQDoubleArray, theoreticalHetSensitivitySampleSize, LOG_ODDS_THRESHOLD);
                HET_SNP_Q = QualityUtil.getPhredScoreFromErrorProbability((1 - HET_SNP_SENSITIVITY));
            }
        }
    }

    public static void main(final String[] args) {
        new CollectWgsMetrics().instanceMainWithExit(args);
    }

    /** Gets the SamReader from which records will be examined.  This will also set the header so that it is available in
     *  */
    protected SamReader getSamReader() {
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        this.header        = in.getFileHeader();
        return in;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
        }

        // it doesn't make sense for the locus accumulation cap to be lower than the coverage cap
        if (LOCUS_ACCUMULATION_CAP < COVERAGE_CAP) {
            log.warn("Setting the LOCUS_ACCUMULATION_CAP to be equal to the COVERAGE_CAP (" + COVERAGE_CAP + ") because it should not be lower");
            LOCUS_ACCUMULATION_CAP = COVERAGE_CAP;
        }

        // Setup all the inputs
        final ProgressLogger progress = new ProgressLogger(log, 10000000, "Processed", "loci");
        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = getSamReader();
        final SamLocusIterator iterator = getLocusIterator(in);

        final List<SamRecordFilter> filters = new ArrayList<>();
        final CountingFilter mapqFilter = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingPairedFilter pairFilter = new CountingPairedFilter();
        // The order in which filters are added matters!
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        if (!COUNT_UNPAIRED) {
            filters.add(pairFilter);
        }
        iterator.setSamFilters(filters);
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setQualityScoreCutoff(0);        // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);
        iterator.setMaxReadsToAccumulatePerLocus(LOCUS_ACCUMULATION_CAP);

        final WgsMetricsCollector collector = getCollector(COVERAGE_CAP, getIntervalsToExamine());

        final boolean usingStopAfter = STOP_AFTER > 0;
        final long stopAfter = STOP_AFTER - 1;
        long counter = 0;

        // Loop through all the loci
        while (iterator.hasNext()) {
            final SamLocusIterator.LocusInfo info = iterator.next();
            final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());

            // Check that the reference is not N
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (SequenceUtil.isNoCall(base)) continue;

            // add to the collector
            collector.addInfo(info, ref);

            // Record progress and perhaps stop
            progress.record(info.getSequenceName(), info.getPosition());
            if (usingStopAfter && ++counter > stopAfter) break;
        }


        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        collector.addToMetricsFile(out, INCLUDE_BQ_HISTOGRAM, dupeFilter, mapqFilter, pairFilter);
        out.write(OUTPUT);

        return 0;
    }

    /** Gets the intervals over which we will calculate metrics. */
    protected IntervalList getIntervalsToExamine() {
        final IntervalList intervals;
        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            intervals = IntervalList.fromFile(INTERVALS);
        } else {
            intervals = new IntervalList(this.header);
            for (final SAMSequenceRecord rec : this.header.getSequenceDictionary().getSequences()) {
                final Interval interval = new Interval(rec.getSequenceName(), 1, rec.getSequenceLength());
                intervals.add(interval);
            }
        }
        return intervals;
    }

    /** This method should only be called after {@link this.getSamReader()} is called. */
    protected SAMFileHeader getSamFileHeader() {
        if (this.header == null) throw new IllegalStateException("getSamFileHeader() was called but this.header is null");
        return this.header;
    }

    protected WgsMetrics generateWgsMetrics(final IntervalList intervals,
                                          final Histogram<Integer> depthHistogram,
                                          final double pctExcludedByMapq,
                                          final double pctExcludedByDupes,
                                          final double pctExcludedByPairing,
                                          final double pctExcludedByBaseq,
                                          final double pctExcludedByOverlap,
                                          final double pctExcludedByCapping,
                                          final double pctTotal,
                                          final int coverageCap,
                                          final Histogram<Integer> baseQHistogram,
                                          final int theoreticalHetSensitivitySampleSize) {
        return new WgsMetrics(
                intervals,
                depthHistogram,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                baseQHistogram,
                theoreticalHetSensitivitySampleSize
        );
    }
    
    private WgsMetrics generateWgsMetrics(final IntervalList intervals,
                                            final Histogram<Integer> depthHistogram,
                                            final long basesExcludedByMapq,
                                            final long basesExcludedByDupes,
                                            final long basesExcludedByPairing,
                                            final long basesExcludedByBaseq,
                                            final long basesExcludedByOverlap,
                                            final long basesExcludedByCapping,
                                            final int coverageCap,
                                            final Histogram<Integer> baseQHistogram,
                                            final int theoreticalHetSensitivitySampleSize) {
        final double total              = depthHistogram.getSum();
        final double totalWithExcludes  = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;

        final double pctExcludedByMapq    = (0 == totalWithExcludes) ? 0 : (basesExcludedByMapq         / totalWithExcludes);
        final double pctExcludedByDupes   = (0 == totalWithExcludes) ? 0 : (basesExcludedByDupes        / totalWithExcludes);
        final double pctExcludedByPairing = (0 == totalWithExcludes) ? 0 : (basesExcludedByPairing      / totalWithExcludes);
        final double pctExcludedByBaseq   = (0 == totalWithExcludes) ? 0 : (basesExcludedByBaseq        / totalWithExcludes);
        final double pctExcludedByOverlap = (0 == totalWithExcludes) ? 0 : (basesExcludedByOverlap      / totalWithExcludes);
        final double pctExcludedByCapping = (0 == totalWithExcludes) ? 0 : (basesExcludedByCapping      / totalWithExcludes);
        final double pctTotal             = (0 == totalWithExcludes) ? 0 : ((totalWithExcludes - total) / totalWithExcludes);

        return generateWgsMetrics(
                intervals,
                depthHistogram,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                baseQHistogram,
                theoreticalHetSensitivitySampleSize
        );
    }

    /**
     * If INTERVALS is specified, this will count bases beyond the interval list when the read overlaps the intervals and extends beyond the
     * edge. Ideally INTERVALS should only include regions that have hard edges without reads that could extend beyond the boundary (such as a whole contig).
     */
    protected long getBasesExcludedBy(final CountingFilter filter) {
        return filter.getFilteredBases();
    }

    protected SamLocusIterator getLocusIterator(final SamReader in) {
        return (INTERVALS != null) ? new SamLocusIterator(in, IntervalList.fromFile(INTERVALS)) : new SamLocusIterator(in);
    }

    /**
     * @param coverageCap the maximum depth/coverage to consider.
     * @param intervals the intervals over which metrics are collected.
     * @return
     */
    protected WgsMetricsCollector getCollector(final int coverageCap, final IntervalList intervals) {
        return new WgsMetricsCollector(coverageCap, intervals);
    }

    protected class WgsMetricsCollector {

        protected final long[] depthHistogramArray;
        private   final long[] baseQHistogramArray;

        private long basesExcludedByBaseq = 0;
        private long basesExcludedByOverlap = 0;
        private long basesExcludedByCapping = 0;
        protected final IntervalList intervals;
        protected final int coverageCap;

        public WgsMetricsCollector(final int coverageCap, final IntervalList intervals) {
            depthHistogramArray = new long[coverageCap + 1];
            baseQHistogramArray = new long[Byte.MAX_VALUE];
            this.coverageCap    = coverageCap;
            this.intervals      = intervals;
        }

        public void addInfo(final SamLocusIterator.LocusInfo info, final ReferenceSequence ref) {

            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<>(info.getRecordAndPositions().size());
            int pileupSize = 0;
            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {

                if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY ||
                        SequenceUtil.isNoCall(recs.getReadBase()))                  { ++basesExcludedByBaseq;   continue; }
                if (!readNames.add(recs.getRecord().getReadName()))                 { ++basesExcludedByOverlap; continue; }

                pileupSize++;
                if (pileupSize <= coverageCap) {
                    baseQHistogramArray[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                }
            }

            final int depth = Math.min(pileupSize, coverageCap);
            if (depth < pileupSize) basesExcludedByCapping += pileupSize - coverageCap;
            depthHistogramArray[depth]++;
        }

        public void addToMetricsFile(final MetricsFile<WgsMetrics, Integer> file,
                                     final boolean includeBQHistogram,
                                     final CountingFilter dupeFilter,
                                     final CountingFilter mapqFilter,
                                     final CountingPairedFilter pairFilter) {
            addMetricsToFile(file, dupeFilter, mapqFilter, pairFilter);

            if (includeBQHistogram) {
                addBaseQHistogram(file);
            }
        }

        protected void addBaseQHistogram(final MetricsFile<WgsMetrics, Integer> file) {
            file.addHistogram(getBaseQHistogram());
        }

        protected void addMetricsToFile(final MetricsFile<WgsMetrics, Integer> file,
                                        final CountingFilter dupeFilter,
                                        final CountingFilter mapqFilter,
                                        final CountingPairedFilter pairFilter) {
            // get the depth histogram and metrics
            final Histogram<Integer> depthHistogram = getDepthHistogram();
            final WgsMetrics metrics = getMetrics(depthHistogram, dupeFilter, mapqFilter, pairFilter);

            // add them to the file
            file.addMetric(metrics);
            file.addHistogram(depthHistogram);
        }

        protected Histogram<Integer> getDepthHistogram() {
            return getHistogram(depthHistogramArray,"coverage", "count");
        }

        protected Histogram<Integer> getBaseQHistogram() {
            return getHistogram(baseQHistogramArray, "value", "baseq_count");
        }

        protected Histogram<Integer> getHistogram(final long[] array, final String binLabel, final String valueLabel) {
            final Histogram<Integer> histogram = new Histogram<>(binLabel, valueLabel);
            for (int i = 0; i < array.length; ++i) {
                histogram.increment(i, array[i]);
            }
            return histogram;
        }

        protected WgsMetrics getMetrics(final Histogram<Integer> depthHistogram,
                                        final CountingFilter dupeFilter,
                                        final CountingFilter mapqFilter,
                                        final CountingPairedFilter pairFilter) {
            return generateWgsMetrics(
                    this.intervals,
                    depthHistogram,
                    getBasesExcludedBy(mapqFilter),
                    getBasesExcludedBy(dupeFilter),
                    getBasesExcludedBy(pairFilter),
                    basesExcludedByBaseq,
                    basesExcludedByOverlap,
                    basesExcludedByCapping,
                    coverageCap,
                    getBaseQHistogram(),
                    SAMPLE_SIZE
            );
        }
    }
}

