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
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.QualityUtil;
import htsjdk.samtools.util.SamLocusIterator;
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
import java.util.stream.IntStream;

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

    @Option(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.", overridable = true)
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
    public static class WgsMetrics extends MetricBase {

        /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
        public long GENOME_TERRITORY;
        /** The mean coverage in bases of the genome territory, after all filters are applied. */
        public double MEAN_COVERAGE;
        /** The standard deviation of coverage of the genome after all filters are applied. */
        public double SD_COVERAGE;
        /** The median coverage in bases of the genome territory, after all filters are applied. */
        public double MEDIAN_COVERAGE;
        /** The median absolute deviation of coverage of the genome after all filters are applied. */
        public double MAD_COVERAGE;

        /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
        public double PCT_EXC_MAPQ;
        /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
        public double PCT_EXC_DUPE;
        /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
        public double PCT_EXC_UNPAIRED;
        /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
        public double PCT_EXC_BASEQ;
        /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
        public double PCT_EXC_OVERLAP;
        /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
        public double PCT_EXC_CAPPED;
        /** The total fraction of aligned bases excluded due to all filters. */
        public double PCT_EXC_TOTAL;

        /** The fraction of bases that attained at least 1X sequence coverage in post-filtering bases. */
        public double PCT_1X;
        /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
        public double PCT_5X;
        /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
        public double PCT_10X;
        /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
        public double PCT_15X;
        /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
        public double PCT_20X;
        /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
        public double PCT_25X;
        /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
        public double PCT_30X;
        /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
        public double PCT_40X;
        /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
        public double PCT_50X;
        /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
        public double PCT_60X;
        /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
        public double PCT_70X;
        /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
        public double PCT_80X;
        /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
        public double PCT_90X;
        /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
        public double PCT_100X;

        /** The theoretical HET SNP sensitivity. */
        public double HET_SNP_SENSITIVITY;

        /** The Phred Scaled Q Score of the theoretical HET SNP sensitivity. */
        public double HET_SNP_Q;
    }

    public static void main(final String[] args) {
        new CollectWgsMetrics().instanceMainWithExit(args);
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
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SamLocusIterator iterator = getLocusIterator(in);
        this.header = in.getFileHeader();

        final List<SamRecordFilter> filters = new ArrayList<>();
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingFilter mapqFilter = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
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

        final WgsMetricsCollector collector = getCollector(COVERAGE_CAP);

        final boolean usingStopAfter = STOP_AFTER > 0;
        final long stopAfter = STOP_AFTER - 1;
        long counter = 0;

        // Loop through all the loci
        while (iterator.hasNext()) {
            final SamLocusIterator.LocusInfo info = iterator.next();
            final ReferenceSequence ref = refWalker.get(info.getSequenceIndex());

            // Check that the reference is not N
            final byte base = ref.getBases()[info.getPosition() - 1];
            if (base == 'N') continue;

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

    protected SAMFileHeader getSamFileHeader() {
        return this.header;
    }

    protected WgsMetrics generateWgsMetrics() {
        return new WgsMetrics();
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

    protected WgsMetricsCollector getCollector(final int coverageCap) {
        return new WgsMetricsCollector(coverageCap);
    }

    protected class WgsMetricsCollector {

        protected final long[] histogramArray;
        private final long[] baseQHistogramArray;
        private final long[] baseQHetSensHistogram;

        private long basesExcludedByBaseq = 0;
        private long basesExcludedByOverlap = 0;
        private long basesExcludedByCapping = 0;
        protected final int coverageCap;

        public WgsMetricsCollector(final int coverageCap) {
            histogramArray = new long[coverageCap + 1];
            baseQHistogramArray = new long[Byte.MAX_VALUE];
            baseQHetSensHistogram = new long[Byte.MAX_VALUE];
            this.coverageCap = coverageCap;
        }

        public void addInfo(final SamLocusIterator.LocusInfo info, final ReferenceSequence ref) {

            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<>(info.getRecordAndPositions().size());
            int pileupSize = 0;
            int pileupSizeForBaseQHetSens = 0;
            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndPositions()) {
                pileupSizeForBaseQHetSens++;
                if(pileupSizeForBaseQHetSens <= coverageCap) {
                    baseQHetSensHistogram[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                }

                if (recs.getBaseQuality() < MINIMUM_BASE_QUALITY)                   { ++basesExcludedByBaseq;   continue; }
                if (!readNames.add(recs.getRecord().getReadName()))                 { ++basesExcludedByOverlap; continue; }

                pileupSize++;
                if (pileupSize <= coverageCap) {
                    baseQHistogramArray[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                }
            }

            final int depth = Math.min(readNames.size(), coverageCap);
            if (depth < readNames.size()) basesExcludedByCapping += readNames.size() - coverageCap;
            histogramArray[depth]++;
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
            // Construct and write the outputs
            final Histogram<Integer> baseQHistogram = new Histogram<>("value", "baseq_count");

            for (int i = 0; i < baseQHistogramArray.length; ++i) {
                baseQHistogram.increment(i, baseQHistogramArray[i]);
            }

            if (INCLUDE_BQ_HISTOGRAM) {
                file.addHistogram(baseQHistogram);
            }
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
            final Histogram<Integer> depthHistogram = new Histogram<>("coverage", "count");
            for (int i=0; i<histogramArray.length; ++i) {
                depthHistogram.increment(i, histogramArray[i]);
            }
            return depthHistogram;
        }

        protected WgsMetrics getMetrics(final Histogram<Integer> depthHistogram,
                                      final CountingFilter dupeFilter,
                                      final CountingFilter mapqFilter,
                                      final CountingPairedFilter pairFilter) {

            // the base q het histogram
            final Histogram<Integer> baseQHetHistogram = new Histogram<>("value", "baseq_count");
            final int BASEQ_MAX = 50;
            final Integer[] x = new Integer[BASEQ_MAX];
            IntStream.range(0, BASEQ_MAX).forEach(i -> x[i] = i);
            baseQHetHistogram.prefillBins(x);

            // Haplotype caller uses 17 as a baseQ cut off, so we are too. Everything below 17 is squashed into the '0' bin.
            final int BASEQ_MIN_CUTOFF = 17;
            for (int i=0; i<baseQHetSensHistogram.length; ++i) {
                baseQHetHistogram.increment( i < BASEQ_MIN_CUTOFF ? 0 : i, baseQHetSensHistogram[i]);
            }

            final WgsMetrics metrics = generateWgsMetrics();
            metrics.GENOME_TERRITORY = (long) depthHistogram.getSumOfValues();
            metrics.MEAN_COVERAGE    = depthHistogram.getMean();
            metrics.SD_COVERAGE      = depthHistogram.getStandardDeviation();
            metrics.MEDIAN_COVERAGE  = depthHistogram.getMedian();
            metrics.MAD_COVERAGE     = depthHistogram.getMedianAbsoluteDeviation();

            final long basesExcludedByDupes   = getBasesExcludedBy(dupeFilter);
            final long basesExcludedByMapq    = getBasesExcludedBy(mapqFilter);
            final long basesExcludedByPairing = getBasesExcludedBy(pairFilter);
            final double total                = depthHistogram.getSum();
            final double totalWithExcludes    = total + basesExcludedByDupes + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;

            metrics.PCT_EXC_DUPE     = basesExcludedByDupes / totalWithExcludes;
            metrics.PCT_EXC_MAPQ     = basesExcludedByMapq / totalWithExcludes;
            metrics.PCT_EXC_UNPAIRED = basesExcludedByPairing / totalWithExcludes;
            metrics.PCT_EXC_BASEQ    = basesExcludedByBaseq   / totalWithExcludes;
            metrics.PCT_EXC_OVERLAP  = basesExcludedByOverlap / totalWithExcludes;
            metrics.PCT_EXC_CAPPED   = basesExcludedByCapping / totalWithExcludes;
            metrics.PCT_EXC_TOTAL    = (totalWithExcludes - total) / totalWithExcludes;


            metrics.PCT_1X    = MathUtil.sum(histogramArray, 1, histogramArray.length)   / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_5X    = MathUtil.sum(histogramArray, 5, histogramArray.length)   / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_10X   = MathUtil.sum(histogramArray, 10, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_15X   = MathUtil.sum(histogramArray, 15, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_20X   = MathUtil.sum(histogramArray, 20, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_25X   = MathUtil.sum(histogramArray, 25, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_30X   = MathUtil.sum(histogramArray, 30, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_40X   = MathUtil.sum(histogramArray, 40, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_50X   = MathUtil.sum(histogramArray, 50, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_60X   = MathUtil.sum(histogramArray, 60, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_70X   = MathUtil.sum(histogramArray, 70, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_80X   = MathUtil.sum(histogramArray, 80, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_90X   = MathUtil.sum(histogramArray, 90, histogramArray.length)  / (double) metrics.GENOME_TERRITORY;
            metrics.PCT_100X  = MathUtil.sum(histogramArray, 100, histogramArray.length) / (double) metrics.GENOME_TERRITORY;

            // Get Theoretical Het SNP Sensitivity
            final double [] depthDoubleArray = TheoreticalSensitivity.normalizeHistogram(depthHistogram);
            final double [] baseQDoubleArray = TheoreticalSensitivity.normalizeHistogram(baseQHetHistogram);
            metrics.HET_SNP_SENSITIVITY = TheoreticalSensitivity.hetSNPSensitivity(depthDoubleArray, baseQDoubleArray, SAMPLE_SIZE, LOG_ODDS_THRESHOLD);
            metrics.HET_SNP_Q = QualityUtil.getPhredScoreFromErrorProbability((1-metrics.HET_SNP_SENSITIVITY));

            return metrics;
        }
    }

}

