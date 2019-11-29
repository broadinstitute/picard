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
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.util.AbstractLocusIterator;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.EdgeReadIterator;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.IntervalArgumentCollection;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.filter.CountingAdapterFilter;
import picard.filter.CountingDuplicateFilter;
import picard.filter.CountingFilter;
import picard.filter.CountingMapQFilter;
import picard.filter.CountingPairedFilter;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments.
 * Two algorithms are available for this metrics: default and fast. The fast algorithm is enabled by USE_FAST_ALGORITHM option.
 * The fast algorithm works better for regions of BAM file with coverage at least 10 reads per locus,
 * for lower coverage the algorithms perform the same.
 * @author tfennell
 */
@CommandLineProgramProperties(
        summary = CollectWgsMetrics.USAGE_SUMMARY + CollectWgsMetrics.USAGE_DETAILS,
        oneLineSummary = CollectWgsMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
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

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Argument(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "Minimum mapping quality for a read to contribute coverage.")
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Argument(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage. N bases will be treated as having a base quality " +
            "of negative infinity and will therefore be excluded from coverage regardless of the value of this parameter.")
    public int MINIMUM_BASE_QUALITY = 20;

    @Argument(shortName = "CAP", doc = "Treat positions with coverage exceeding this value as if they had coverage at this value (but calculate the difference for PCT_EXC_CAPPED).")
    public int COVERAGE_CAP = 250;

    @Argument(doc="At positions with coverage exceeding this value, completely ignore reads that accumulate beyond this value (so that they will not be considered for PCT_EXC_CAPPED).  Used to keep memory consumption in check, but could create bias if set too low")
    public int LOCUS_ACCUMULATION_CAP = 100000;

    @Argument(doc = "For debugging purposes, stop after processing this many genomic bases.")
    public long STOP_AFTER = -1;

    @Argument(doc = "Determines whether to include the base quality histogram in the metrics file.")
    public boolean INCLUDE_BQ_HISTOGRAM = false;

    @Argument(doc="If true, count unpaired reads, and paired reads with one end unmapped")
    public boolean COUNT_UNPAIRED = false;

    @Argument(doc="Sample Size used for Theoretical Het Sensitivity sampling. Default is 10000.", optional = true)
    public int SAMPLE_SIZE=10000;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArugmentCollection = makeIntervalArgumentCollection();

    @Argument(doc="Output for Theoretical Sensitivity metrics.", optional = true)
    public File THEORETICAL_SENSITIVITY_OUTPUT;

    @Argument(doc="Allele fraction for which to calculate theoretical sensitivity.", optional = true)
    public List<Double> ALLELE_FRACTION = new ArrayList<>(Arrays.asList(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5));

    @Argument(doc = "If true, fast algorithm is used.")
    public boolean USE_FAST_ALGORITHM = false;

    @Argument(doc = "Average read length in the file. Default is 150.", optional = true)
    public int READ_LENGTH = 150;

    protected File INTERVALS = null;

    private SAMFileHeader header = null;

    private final Log log = Log.getInstance(CollectWgsMetrics.class);

    @Override
    protected boolean requiresReference() {
        return true;
    }

    /**
     * @return An interval argument collection to be used for this tool. Subclasses can override this
     * to provide an argument collection with alternative arguments or argument annotations.
     */
    protected IntervalArgumentCollection makeIntervalArgumentCollection() {
        return new CollectWgsMetricsIntervalArgumentCollection();
    }

    public static class CollectWgsMetricsIntervalArgumentCollection implements IntervalArgumentCollection {
        @Argument(doc = "An interval list file that contains the positions to restrict the assessment. Please note that " +
                "all bases of reads that overlap these intervals will be considered, even if some of those bases extend beyond the boundaries of " +
                "the interval. The ideal use case for this argument is to use it to restrict the calculation to a subset of (whole) contigs.",
                optional = true)
        public File INTERVALS;

        public File getIntervalFile() { return INTERVALS; }
    };



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
        INTERVALS = intervalArugmentCollection.getIntervalFile();
        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
        }
        if (THEORETICAL_SENSITIVITY_OUTPUT != null) {
            IOUtil.assertFileIsWritable(THEORETICAL_SENSITIVITY_OUTPUT);
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
        final AbstractLocusIterator iterator = getLocusIterator(in);

        final List<SamRecordFilter> filters = new ArrayList<>();
        final CountingFilter adapterFilter = new CountingAdapterFilter();
        final CountingFilter mapqFilter = new CountingMapQFilter(MINIMUM_MAPPING_QUALITY);
        final CountingFilter dupeFilter = new CountingDuplicateFilter();
        final CountingPairedFilter pairFilter = new CountingPairedFilter();
        // The order in which filters are added matters!
        filters.add(new SecondaryAlignmentFilter()); // Not a counting filter because we never want to count reads twice
        filters.add(adapterFilter);
        filters.add(mapqFilter);
        filters.add(dupeFilter);
        if (!COUNT_UNPAIRED) {
            filters.add(pairFilter);
        }
        iterator.setSamFilters(filters);
        iterator.setMappingQualityScoreCutoff(0); // Handled separately because we want to count bases
        iterator.setIncludeNonPfReads(false);

        final AbstractWgsMetricsCollector<?> collector = getCollector(COVERAGE_CAP, getIntervalsToExamine());
        final WgsMetricsProcessor processor = getWgsMetricsProcessor(progress, refWalker, iterator, collector);
        processor.processFile();

        final MetricsFile<WgsMetrics, Integer> out = getMetricsFile();
        processor.addToMetricsFile(out, INCLUDE_BQ_HISTOGRAM, dupeFilter, adapterFilter, mapqFilter, pairFilter);
        out.write(OUTPUT);

        if (THEORETICAL_SENSITIVITY_OUTPUT != null) {
            // Write out theoretical sensitivity results.
            final MetricsFile<TheoreticalSensitivityMetrics, ?> theoreticalSensitivityMetrics = getMetricsFile();
            log.info("Calculating theoretical sentitivity at " + ALLELE_FRACTION.size() + " allele fractions.");
            List<TheoreticalSensitivityMetrics> tsm = TheoreticalSensitivity.calculateSensitivities(SAMPLE_SIZE, collector.getUnfilteredDepthHistogram(), collector.getUnfilteredBaseQHistogram(), ALLELE_FRACTION);
            theoreticalSensitivityMetrics.addAllMetrics(tsm);
            theoreticalSensitivityMetrics.write(THEORETICAL_SENSITIVITY_OUTPUT);
        }

        return 0;
    }

    private <T extends AbstractRecordAndOffset> WgsMetricsProcessorImpl<T> getWgsMetricsProcessor(
            ProgressLogger progress, ReferenceSequenceFileWalker refWalker,
            AbstractLocusIterator<T, AbstractLocusInfo<T>> iterator, AbstractWgsMetricsCollector<T> collector) {
        return new WgsMetricsProcessorImpl<>(iterator, refWalker, collector, progress);
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
                                            final Histogram<Integer> highQualityDepthHistogram,
                                            final Histogram<Integer> unfilteredDepthHistogram,
                                            final double pctExcludedByAdapter,
                                            final double pctExcludedByMapq,
                                            final double pctExcludedByDupes,
                                            final double pctExcludedByPairing,
                                            final double pctExcludedByBaseq,
                                            final double pctExcludedByOverlap,
                                            final double pctExcludedByCapping,
                                            final double pctTotal,
                                            final int coverageCap,
                                            final Histogram<Integer> unfilteredBaseQHistogram,
                                            final int theoreticalHetSensitivitySampleSize) {
        return new WgsMetrics(
                intervals,
                highQualityDepthHistogram,
                unfilteredDepthHistogram,
                pctExcludedByAdapter,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                unfilteredBaseQHistogram,
                theoreticalHetSensitivitySampleSize
        );
    }

    WgsMetrics generateWgsMetrics(final IntervalList intervals,
                                          final Histogram<Integer> highQualityDepthHistogram,
                                          final Histogram<Integer> unfilteredDepthHistogram,
                                          final long basesExcludedByAdapter,
                                          final long basesExcludedByMapq,
                                          final long basesExcludedByDupes,
                                          final long basesExcludedByPairing,
                                          final long basesExcludedByBaseq,
                                          final long basesExcludedByOverlap,
                                          final long basesExcludedByCapping,
                                          final int coverageCap,
                                          final Histogram<Integer> unfilteredBaseQHistogram,
                                          final int theoreticalHetSensitivitySampleSize) {
        final double total = highQualityDepthHistogram.getSum();
        final double totalWithExcludes = total + basesExcludedByDupes + basesExcludedByAdapter + basesExcludedByMapq + basesExcludedByPairing + basesExcludedByBaseq + basesExcludedByOverlap + basesExcludedByCapping;

        final double pctExcludedByAdapter = totalWithExcludes == 0 ? 0 : basesExcludedByAdapter      / totalWithExcludes;
        final double pctExcludedByMapq    = totalWithExcludes == 0 ? 0 : basesExcludedByMapq         / totalWithExcludes;
        final double pctExcludedByDupes   = totalWithExcludes == 0 ? 0 : basesExcludedByDupes        / totalWithExcludes;
        final double pctExcludedByPairing = totalWithExcludes == 0 ? 0 : basesExcludedByPairing      / totalWithExcludes;
        final double pctExcludedByBaseq   = totalWithExcludes == 0 ? 0 : basesExcludedByBaseq        / totalWithExcludes;
        final double pctExcludedByOverlap = totalWithExcludes == 0 ? 0 : basesExcludedByOverlap      / totalWithExcludes;
        final double pctExcludedByCapping = totalWithExcludes == 0 ? 0 : basesExcludedByCapping      / totalWithExcludes;
        final double pctTotal             = totalWithExcludes == 0 ? 0 : (totalWithExcludes - total) / totalWithExcludes;

        return generateWgsMetrics(
                intervals,
                highQualityDepthHistogram,
                unfilteredDepthHistogram,
                pctExcludedByAdapter,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                unfilteredBaseQHistogram,
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

    /**
     * Creates {@link htsjdk.samtools.util.AbstractLocusIterator} implementation according to {@link this#USE_FAST_ALGORITHM} value.
     *
     * @param in inner {@link htsjdk.samtools.SamReader}
     * @return if {@link this#USE_FAST_ALGORITHM} is enabled, returns {@link htsjdk.samtools.util.EdgeReadIterator} implementation,
     * otherwise default algorithm is used and {@link htsjdk.samtools.util.SamLocusIterator} is returned.
     */
    protected AbstractLocusIterator getLocusIterator(final SamReader in) {
        if (USE_FAST_ALGORITHM) {
            return (INTERVALS != null) ? new EdgeReadIterator(in, IntervalList.fromFile(INTERVALS)) : new EdgeReadIterator(in);
        }
        SamLocusIterator iterator = (INTERVALS != null) ? new SamLocusIterator(in, IntervalList.fromFile(INTERVALS)) : new SamLocusIterator(in);
        iterator.setMaxReadsToAccumulatePerLocus(LOCUS_ACCUMULATION_CAP);
        iterator.setEmitUncoveredLoci(true);
        iterator.setQualityScoreCutoff(0);
        return iterator;
    }

    /**
     * Creates {@link picard.analysis.AbstractWgsMetricsCollector} implementation according to {@link this#USE_FAST_ALGORITHM} value.
     *
     * @param coverageCap the maximum depth/coverage to consider.
     * @param intervals the intervals over which metrics are collected.
     * @return if {@link this#USE_FAST_ALGORITHM} is enabled, returns {@link picard.analysis.FastWgsMetricsCollector} implementation,
     * otherwise default algorithm is used and {@link picard.analysis.CollectWgsMetrics.WgsMetricsCollector} is returned.
     */
    protected AbstractWgsMetricsCollector getCollector(final int coverageCap, final IntervalList intervals) {
        return USE_FAST_ALGORITHM ? new FastWgsMetricsCollector(this, coverageCap, intervals) :
                new WgsMetricsCollector(this, coverageCap, intervals);
    }

    protected static class WgsMetricsCollector extends AbstractWgsMetricsCollector<SamLocusIterator.RecordAndOffset> {

        public WgsMetricsCollector(final CollectWgsMetrics metrics, final int coverageCap, final IntervalList intervals) {
            super(metrics, coverageCap, intervals);
        }

        @Override
        public void addInfo(final AbstractLocusInfo<SamLocusIterator.RecordAndOffset> info, final ReferenceSequence ref, boolean referenceBaseN) {

            if (referenceBaseN) {
                return;
            }
            // Figure out the coverage while not counting overlapping reads twice, and excluding various things
            final HashSet<String> readNames = new HashSet<>(info.getRecordAndOffsets().size());
            int pileupSize = 0;
            int unfilteredDepth = 0;

            for (final SamLocusIterator.RecordAndOffset recs : info.getRecordAndOffsets()) {
                if (recs.getBaseQuality() <= 2) { ++basesExcludedByBaseq;   continue; }

                // we add to the base quality histogram any bases that have quality > 2
                // the raw depth may exceed the coverageCap before the high-quality depth does. So stop counting once we reach the coverage cap.
                if (unfilteredDepth < coverageCap) {
                    unfilteredBaseQHistogramArray[recs.getRecord().getBaseQualities()[recs.getOffset()]]++;
                    unfilteredDepth++;
                }

                if (recs.getBaseQuality() < collectWgsMetrics.MINIMUM_BASE_QUALITY ||
                        SequenceUtil.isNoCall(recs.getReadBase())) {
                    ++basesExcludedByBaseq;
                    continue;
                }
                if (!readNames.add(recs.getRecord().getReadName())) {
                    ++basesExcludedByOverlap;
                    continue;
                }
                pileupSize++;
            }
            final int highQualityDepth = Math.min(pileupSize, coverageCap);
            if (highQualityDepth < pileupSize) basesExcludedByCapping += pileupSize - coverageCap;
            highQualityDepthHistogramArray[highQualityDepth]++;
            unfilteredDepthHistogramArray[unfilteredDepth]++;
        }
    }
}
