/*
 * The MIT License
 *
 * Copyright (c) 2021 The Broad Institute
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
package picard.util;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.CollectQualityYieldMetrics;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.List;

/**
 * Combines multiple Picard QualityYieldMetrics files into a single file.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary = "Combines multiple QualityYieldMetrics files into a single file. This tool is used in cases where the metrics are calculated" +
                " separately on shards of the same read-group.",
        oneLineSummary = "Combines multiple QualityYieldMetrics files into a single file.",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class AccumulateQualityYieldMetrics extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input QualityYieldMetrics files to merge.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output QualityYieldMetric file to write.")
    public File OUTPUT;

    @Override
    protected int doWork() {

        IOUtil.assertFileIsWritable(OUTPUT);

        // set up the output metric
        // note that useOriginalQualities does not matter here
        CollectQualityYieldMetrics.QualityYieldMetrics finalMetric = new CollectQualityYieldMetrics.QualityYieldMetrics(false);

        INPUT.forEach(file -> finalMetric.merge(MetricsFile.<CollectQualityYieldMetrics.QualityYieldMetrics>readBeans(file).get(0)));

        finalMetric.calculateDerivedFields();

        MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> outputMetricsFile = new MetricsFile<>();
        outputMetricsFile.addMetric(finalMetric);
        outputMetricsFile.write(OUTPUT);
        return 0;
    }
}
