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

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Collects a set of HS metrics from a sam or bam file.  See HsMetricsCollector and CollectTargetedMetrics for more details.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = CollectHsMetrics.USAGE_SUMMARY + CollectHsMetrics.USAGE_DETAILS,
        usageShort = CollectHsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectHsMetrics extends CollectTargetedMetrics<HsMetrics, HsMetricCollector> {
        static final String USAGE_SUMMARY = "<b>Collects hybrid-selection (HS) metrics for a SAM or BAM file.</b>  ";
        static final String USAGE_DETAILS = "This tool takes a SAM/BAM file input and collects metrics that are specific for sequence " +
                "datasets generated through hybrid-selection. Hybrid-selection (HS) is the most commonly used technique to capture " +
                "exon-specific sequences for targeted sequencing experiments such as exome sequencing; for more information, please " +
                "see the corresponding <a href='http://www.broadinstitute.org/gatk/guide/article?id=6331'>GATK Dictionary entry</a>. </p> "+

                "<p>This tool requires an aligned SAM or BAM file as well as bait and target interval files in Picard interval_list format. " +
                "You should use the bait and interval files that correspond to the capture kit that was used to generate the capture " +
                "libraries for sequencing, which can generally be obtained from the manufacturer of the kit. If the baits and target " +
                "intervals are provided in BED format, you can convert them to the Picard interval_list format using Picard's " +
                "<a href='http://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList'>BedToInterval</a> tool. </p>" +

                "<p>If a reference sequence is provided, this program will calculate both AT_DROPOUT and GC_DROPOUT metrics. Dropout " +
                "metrics are an attempt to measure the reduced representation of reads, in regions that deviate from 50% G/C content. " +
                "This reduction in the number of aligned reads is due to the increased numbers of errors associated with sequencing " +
                "regions with excessive or deficient numbers of G/C bases, ultimately leading to poor mapping efficiencies and low" +
                "coverage in the affected regions. </p>" +

                "<p>If you are interested in getting G/C content and mean sequence depth information for every target interval, use the " +
                "PER_TARGET_COVERAGE option. </p>" +

                "<p>Note that coverage measurements are capped at ~32K to constrain memory usage.</p>"+

                "<h4>Usage Example:</h4>"+
                "<pre>" +
                "java -jar picard.jar CollectHsMetrics \\<br />" +
                "      I=input.bam \\<br />" +
                "      O=hs_metrics.txt \\<br />" +
                "      R=reference_sequence.fasta \\<br />" +
                "      BAIT_INTERVALS=bait.interval_list \\<br />" +
                "      TARGET_INTERVALS=target.interval_list" +
                "</pre> "   +

                "<p>Please see the CollectHsMetrics " +
                "<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics'>definitions</a> for a " +
                "complete description of the metrics produced by this tool.</p> "+
                "<hr />";

    @Option(shortName = "BI", doc = "An interval list file that contains the locations of the baits used.", minElements=1)
    public List<File> BAIT_INTERVALS;

    @Option(shortName = "N", doc = "Bait set name. If not provided it is inferred from the filename of the bait intervals.", optional = true)
    public String BAIT_SET_NAME;

    @Option(shortName = "MQ", doc = "Minimum mapping quality for a read to contribute coverage.", overridable = true)
    public int MINIMUM_MAPPING_QUALITY = 20;

    @Option(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.", overridable = true)
    public int MINIMUM_BASE_QUALITY = 20;

    @Option(doc = "True if we are to clip overlapping reads, false otherwise.", optional=true, overridable = true)
    public boolean CLIP_OVERLAPPING_READS = true;

    @Override
    protected IntervalList getProbeIntervals() {
        for (final File file : BAIT_INTERVALS) IOUtil.assertFileIsReadable(file);
        return IntervalList.fromFiles(BAIT_INTERVALS);
    }

    @Override
    protected String getProbeSetName() {
        if (BAIT_SET_NAME != null) {
            return BAIT_SET_NAME;
        } else {
            final SortedSet<String> baitSetNames = new TreeSet<String>();
            for (final File file : BAIT_INTERVALS) {
                baitSetNames.add(CollectTargetedMetrics.renderProbeNameFromFile(file));
            }
            return StringUtil.join(".", baitSetNames);
        }
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    @Override
    protected HsMetricCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                              final List<SAMReadGroupRecord> samRgRecords,
                                              final ReferenceSequenceFile refFile,
                                              final File perTargetCoverage,
                                              final File perBaseCoverage,
                                              final IntervalList targetIntervals,
                                              final IntervalList probeIntervals,
                                              final String probeSetName,
                                              final int nearProbeDistance) {
        return new HsMetricCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance,
                MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, CLIP_OVERLAPPING_READS, true, COVERAGE_CAP, SAMPLE_SIZE);
    }
}