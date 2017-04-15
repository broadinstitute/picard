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

import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing
 * experiments, same implementation as CollectWgsMetrics, with different defaults: lacks baseQ and mappingQ filters
 * and has much higher coverage cap.
 *
 * @author farjoun
 */
@CommandLineProgramProperties(
        summary = CollectRawWgsMetrics.USAGE_SUMMARY + CollectRawWgsMetrics.USAGE_DETAILS,
        oneLineSummary = CollectRawWgsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectRawWgsMetrics extends CollectWgsMetrics{
    static final String USAGE_SUMMARY = "Collect whole genome sequencing-related metrics.  ";
    static final String USAGE_DETAILS = "This tool computes metrics that are useful for evaluating coverage and performance " +
            "of whole genome sequencing experiments. These metrics include the percentages of reads that pass" +
            " minimal base- and mapping- quality filters as well as coverage (read-depth) levels. " +
            "<br /><br />  " +
            "The histogram output is optional and for a given run, displays two separate outputs on the y-axis while using a single set" +
            " of values for the x-axis.  Specifically, the first column in the histogram table (x-axis) is labeled 'coverage' and " +
            "represents different possible coverage depths.  However, it also represents the range of values for the base quality scores " +
            "and thus should probably be labeled 'sequence depth and base quality scores'. The second and third columns (y-axes) " +
            "correspond to the numbers of bases at a specific sequence depth 'count' and the numbers of bases at a particular base " +
            "quality score 'baseq_count' respectively." +
            "<br /><br />" +
            "Although similar to the CollectWgsMetrics tool, the default thresholds for CollectRawWgsMetrics are less stringent.  " +
            "For example, the CollectRawWgsMetrics have base and mapping quality score thresholds set to '3' and '0' respectively, " +
            "while the CollectWgsMetrics tool has the default threshold values set to '20' (at time of writing).  Nevertheless, both " +
            "tools enable the user to input specific threshold values." +
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectRawWgsMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=raw_wgs_metrics.txt \\<br />" +
            "      R=reference_sequence.fasta \\<br />" +
            "      INCLUDE_BQ_HISTOGRAM=true" +
            "</pre>" +
            "<hr />" +
            "Please see " +
            "<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectWgsMetrics.WgsMetrics'>" +
            "the WgsMetrics documentation</a> for detailed explanations of the output metrics." +
            "<hr />";

    public CollectRawWgsMetrics() {
        //Override default values inherited from CollectWgsMetrics
        MINIMUM_MAPPING_QUALITY = 0;
        MINIMUM_BASE_QUALITY = 3;
        COVERAGE_CAP = 100000;
        LOCUS_ACCUMULATION_CAP = 200000;
    }

    // rename the class so that in the metric file it is annotated differently.
    public static class RawWgsMetrics extends WgsMetrics {
        public RawWgsMetrics() {
            super();
        }

        public RawWgsMetrics(final IntervalList intervals,
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
        return new RawWgsMetrics(
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

}
