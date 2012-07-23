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

package net.sf.picard.analysis.directed;

import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;

import java.io.File;
import java.util.*;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMReadGroupRecord;

/**
 * Calculates a set of HS metrics from a sam or bam file.  See HsMetricsCollector and CollectTargetedMetrics for more details.
 *
 * @author Tim Fennell
 */
public class CalculateHsMetrics extends CollectTargetedMetrics {

    @Usage public final String USAGE =
            "Calculates a set of Hybrid Selection specific metrics from an aligned SAM" +
            "or BAM file. If a reference sequence is provided, AT/GC dropout metrics will " +
            "be calculated, and the PER_TARGET_COVERAGE option can be used to output GC and " +
            "mean coverage information for every target.";
    @Option(shortName="BI", doc="An interval list file that contains the locations of the baits used.")
    public File BAIT_INTERVALS;

    @Option(shortName="N",  doc="Bait set name. If not provided it is inferred from the filename of the bait intervals.", optional=true)
    public String BAIT_SET_NAME;

    /**
     * @return BAIT_INTERVALS file
     */
    @Override
    protected File getProbeIntervals() {
        return BAIT_INTERVALS;
    }

    /**
     * @return BAIT_SET_NAME
     */
    @Override
    protected String getProbeSetName() {
        return BAIT_SET_NAME;
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    @Override
    protected TargetMetricsCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                                    final List<SAMReadGroupRecord> samRgRecords,
                                                    final ReferenceSequenceFile refFile,
                                                    final File perTargetCoverage,
                                                    final File targetIntervals,
                                                    final File probeIntervals,
                                                    final String probeSetName) {
        return new HsMetricCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName);
    }
}
