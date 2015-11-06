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
import htsjdk.samtools.util.IntervalList;
import picard.analysis.MetricAccumulationLevel;

import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * Calculates HS metrics for a given SAM or BAM file. Requires the input of a list of
 * target intervals and a list of bait intervals. Can be invoked either on an entire
 * iterator of SAMRecords or be passed SAMRecords one at a time.
 *
 * @author Jonathan Burke
 */
public class TargetedPcrMetricsCollector extends TargetMetricsCollector<TargetedPcrMetrics> {
    //maybe instead just inject this into the TargetedMetricCollector ->

    public TargetedPcrMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                       final List<SAMReadGroupRecord> samRgRecords,
                                       final ReferenceSequenceFile refFile,
                                       final File perTargetCoverage,
                                       final IntervalList targetIntervals,
                                       final IntervalList probeIntervals,
                                       final String probeSetName,
                                       final int minimumMappingQuality,
                                       final int minimumBaseQuality) {
        super(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName, minimumMappingQuality, minimumBaseQuality);
    }

    public TargetedPcrMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                       final List<SAMReadGroupRecord> samRgRecords,
                                       final ReferenceSequenceFile refFile,
                                       final File perTargetCoverage,
                                       final IntervalList targetIntervals,
                                       final IntervalList probeIntervals,
                                       final String probeSetName,
                                       final int minimumMappingQuality,
                                       final int minimumBaseQuality,
                                       final boolean noSideEffects) {
        super(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName, minimumMappingQuality, minimumBaseQuality, noSideEffects);
    }

    @Override
    public TargetedPcrMetrics convertMetric(TargetMetrics targetMetrics) {
        final TargetedPcrMetrics pcrMetrics = new TargetedPcrMetrics();
        TargetMetricsCollector.reflectiveCopy(targetMetrics, pcrMetrics,
                new String[]{"PROBE_SET",           "PROBE_TERRITORY",    "ON_PROBE_BASES",    "NEAR_PROBE_BASES",    "OFF_PROBE_BASES",    "PCT_SELECTED_BASES",  "PCT_OFF_PROBE",    "ON_PROBE_VS_SELECTED",    "MEAN_PROBE_COVERAGE"},
                new String[]{"CUSTOM_AMPLICON_SET", "AMPLICON_TERRITORY", "ON_AMPLICON_BASES", "NEAR_AMPLICON_BASES", "OFF_AMPLICON_BASES", "PCT_AMPLIFIED_BASES", "PCT_OFF_AMPLICON", "ON_AMPLICON_VS_SELECTED", "MEAN_AMPLICON_COVERAGE"}
        );

        return pcrMetrics;
    }
}
