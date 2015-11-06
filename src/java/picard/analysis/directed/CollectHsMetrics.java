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
        usage = "Collects a set of Hybrid Selection specific metrics from an aligned SAM" +
                "or BAM file. If a reference sequence is provided, AT/GC dropout metrics will " +
                "be calculated, and the PER_TARGET_COVERAGE option can be used to output GC and " +
                "mean coverage information for every target.",
        usageShort = "Collects Hybrid Selection-specific metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectHsMetrics extends CollectTargetedMetrics<HsMetrics, HsMetricCollector> {

    @Option(shortName = "BI", doc = "An interval list file that contains the locations of the baits used.", minElements=1)
    public List<File> BAIT_INTERVALS;

    @Option(shortName = "N", doc = "Bait set name. If not provided it is inferred from the filename of the bait intervals.", optional = true)
    public String BAIT_SET_NAME;

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
                                              final IntervalList targetIntervals,
                                              final IntervalList probeIntervals,
                                              final String probeSetName) {
        return new HsMetricCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName,
                MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, true);
    }
}