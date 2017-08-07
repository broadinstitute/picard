/*
 * The MIT License
 *
 * Copyright (c) 2016 Nils Homer
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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Alpha;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;
import picard.util.RExecutor;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        summary = CollectWgsMetricsWithNonZeroCoverage.USAGE_SUMMARY + CollectWgsMetricsWithNonZeroCoverage.USAGE_DETAILS,
        oneLineSummary = CollectWgsMetricsWithNonZeroCoverage.USAGE_SUMMARY,
        programGroup = Alpha.class
)
public class CollectWgsMetricsWithNonZeroCoverage extends CollectWgsMetrics {

    static final String USAGE_SUMMARY = "Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.  ";
    static final String USAGE_DETAILS = "This tool collects metrics about the percentages of reads that pass base- and mapping- quality " +
            "filters as well as coverage (read-depth) levels. Both minimum base- and mapping-quality values as well as the maximum " +
            "read depths (coverage cap) are user defined.  This extends CollectWgsMetrics by including metrics related only to sites" +
            "with non-zero (>0) coverage." +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage Example:</h4>" +
            "<pre>"  +
            "java -jar picard.jar CollectWgsMetricsWithNonZeroCoverage \\<br /> " +
            "      I=input.bam \\<br /> "+
            "      O=collect_wgs_metrics.txt \\<br /> " +
            "      CHART=collect_wgs_metrics.pdf  \\<br /> " +
            "      R=reference_sequence.fasta " +
            "</pre>" +
            "Please see the " +
            "<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetricsWithNonZeroCoverage.WgsMetricsWithNonZeroCoverage'>" +
            "WgsMetricsWithNonZeroCoverage</a> documentation for detailed explanations of the output metrics." +
            "<hr />";

    @Argument(shortName = "CHART", doc = "A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    private final Log log = Log.getInstance(CollectWgsMetricsWithNonZeroCoverage.class);

    // Store this here since we need access to it in the doWork method
    private WgsMetricsWithNonZeroCoverageCollector collector = null;

    private SamReader samReader = null;

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class WgsMetricsWithNonZeroCoverage extends WgsMetrics {
        public enum Category { WHOLE_GENOME, NON_ZERO_REGIONS }

        /** One of either WHOLE_GENOME or NON_ZERO_REGIONS */
        @MergeByAssertEquals
        public Category CATEGORY;

        public WgsMetricsWithNonZeroCoverage() {
            super();
        }

        public WgsMetricsWithNonZeroCoverage(final IntervalList intervals,
                                             final Histogram<Integer> highQualityDepthHistogram,
                                             final Histogram<Integer> unfilteredDepthHistogram,
                                             final double pctExcludedByMapq,
                                             final double pctExcludedByDupes,
                                             final double pctExcludedByPairing,
                                             final double pctExcludedByBaseq,
                                             final double pctExcludedByOverlap,
                                             final double pctExcludedByCapping,
                                             final double pctTotal,
                                             final int coverageCap,
                                             final Histogram<Integer> unfilteredBaseQHistogram,
                                             final int sampleSize) {
            super(intervals, highQualityDepthHistogram, unfilteredDepthHistogram, pctExcludedByMapq, pctExcludedByDupes, pctExcludedByPairing, pctExcludedByBaseq,
                    pctExcludedByOverlap, pctExcludedByCapping, pctTotal, coverageCap, unfilteredBaseQHistogram, sampleSize);
        }
    }

    public static void main(final String[] args) {
        new CollectWgsMetrics().instanceMainWithExit(args);
    }

    @Override
    protected SamReader getSamReader() {
        if (this.samReader == null) {
            this.samReader = super.getSamReader();
        }
        return this.samReader;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        IOUtil.assertFileIsReadable(INPUT);

        // Initialize the SamReader, so the header is available prior to super.doWork, for getIntervalsToExamine call. */
        getSamReader();

        this.collector = new WgsMetricsWithNonZeroCoverageCollector(this, COVERAGE_CAP, getIntervalsToExamine());
        super.doWork();

        final List<SAMReadGroupRecord> readGroups = getSamFileHeader().getReadGroups();
        final String plotSubtitle = (readGroups.size() == 1) ? StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary()) : "";

        if (collector.areHistogramsEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        } else {
            final int rResult = RExecutor.executeFromClasspath("picard/analysis/wgsHistogram.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName(),
                    plotSubtitle);
            if (rResult != 0) {
                throw new PicardException("R script wgsHistogram.R failed with return code " + rResult);
            }
        }

        return 0;
    }

    @Override
    protected WgsMetrics generateWgsMetrics(final IntervalList intervals,
                                            final Histogram<Integer> highQualityDepthHistogram,
                                            final Histogram<Integer> unfilteredDepthHistogram,
                                            final double pctExcludedByMapq,
                                            final double pctExcludedByDupes,
                                            final double pctExcludedByPairing,
                                            final double pctExcludedByBaseq,
                                            final double pctExcludedByOverlap,
                                            final double pctExcludedByCapping,
                                            final double pctTotal,
                                            final int coverageCap,
                                            final Histogram<Integer> unfilteredBaseQHistogram,
                                            final int sampleSize) {
        return new WgsMetricsWithNonZeroCoverage(
                intervals,
                highQualityDepthHistogram,
                unfilteredDepthHistogram,
                pctExcludedByMapq,
                pctExcludedByDupes,
                pctExcludedByPairing,
                pctExcludedByBaseq,
                pctExcludedByOverlap,
                pctExcludedByCapping,
                pctTotal,
                coverageCap,
                unfilteredBaseQHistogram,
                sampleSize);
    }

    @Override
    protected WgsMetricsCollector getCollector(final int coverageCap, final IntervalList intervals) {
        assert(coverageCap == this.collector.coverageCap);
        return this.collector;
    }

    protected class WgsMetricsWithNonZeroCoverageCollector extends WgsMetricsCollector {
        Histogram<Integer> highQualityDepthHistogram;
        Histogram<Integer> highQualityDepthHistogramNonZero;

        public WgsMetricsWithNonZeroCoverageCollector(final CollectWgsMetricsWithNonZeroCoverage metrics,
                final int coverageCap, final IntervalList intervals) {
            super(metrics, coverageCap, intervals);
        }

        @Override
        public void addToMetricsFile(final MetricsFile<WgsMetrics, Integer> file,
                                     final boolean includeBQHistogram,
                                     final CountingFilter dupeFilter,
                                     final CountingFilter mapqFilter,
                                     final CountingPairedFilter pairFilter) {
            highQualityDepthHistogram = getDepthHistogram();
            highQualityDepthHistogramNonZero = getDepthHistogramNonZero();

            // calculate metrics the same way as in CollectWgsMetrics
            final WgsMetricsWithNonZeroCoverage metrics = (WgsMetricsWithNonZeroCoverage) getMetrics(dupeFilter, mapqFilter, pairFilter);
            metrics.CATEGORY = WgsMetricsWithNonZeroCoverage.Category.WHOLE_GENOME;

            // set count of the coverage-zero bin to 0 and re-calculate metrics
            // note we don't need to update the base quality histogram; there are no bases in the depth = 0 bin
            highQualityDepthHistogramArray[0] = 0;
            unfilteredDepthHistogramArray[0] = 0;

            final WgsMetricsWithNonZeroCoverage metricsNonZero = (WgsMetricsWithNonZeroCoverage) getMetrics(dupeFilter, mapqFilter, pairFilter);
            metricsNonZero.CATEGORY = WgsMetricsWithNonZeroCoverage.Category.NON_ZERO_REGIONS;

            file.addMetric(metrics);
            file.addMetric(metricsNonZero);
            file.addHistogram(highQualityDepthHistogram);
            file.addHistogram(highQualityDepthHistogramNonZero);

            if (includeBQHistogram) {
                addBaseQHistogram(file);
            }
        }

        protected Histogram<Integer> getDepthHistogram() {
            return getHistogram(highQualityDepthHistogramArray, "coverage", "count_WHOLE_GENOME");
        }

        private Histogram<Integer> getDepthHistogramNonZero() {
            final Histogram<Integer> depthHistogram = new Histogram<>("coverage", "count_NON_ZERO_REGIONS");
            // do not include the zero-coverage bin
            for (int i = 1; i < highQualityDepthHistogramArray.length; ++i) {
                depthHistogram.increment(i, highQualityDepthHistogramArray[i]);
            }
            return depthHistogram;
        }

        public boolean areHistogramsEmpty() {
            return (highQualityDepthHistogram.isEmpty() || highQualityDepthHistogramNonZero.isEmpty());
        }
    }
}
