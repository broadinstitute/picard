/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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
package picard.vcf;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Combines multiple Variant Calling Metrics files into a single file.
 * @author Eric Banks
 */
@CommandLineProgramProperties(
        summary = "Combines multiple Variant Calling Metrics files into a single file.  This tool is used in cases where the metrics are calculated" +
                " separately for different (genomic) shards of the same callset and we want to combine them into a single result over the entire callset." +
                " The shards are expected to contain the same samples (although it will not fail if they do not) and to not have been run over overlapping genomic positions.",
        oneLineSummary = "Combines multiple Variant Calling Metrics files into a single file",
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class AccumulateVariantCallingMetrics extends CommandLineProgram {

    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Paths (except for the file extensions) of Variant Calling Metrics files to read and merge.", minElements=1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Path (except for the file extension) of output metrics files to write.")
    public File OUTPUT;

    @Override
    protected int doWork() {

        final String outputPrefix = OUTPUT.getAbsolutePath() + ".";
        final File detailOutputFile = new File(outputPrefix + CollectVariantCallingMetrics.VariantCallingDetailMetrics.getFileExtension());
        final File summaryOutputFile = new File(outputPrefix + CollectVariantCallingMetrics.VariantCallingSummaryMetrics.getFileExtension());
        IOUtil.assertFileIsWritable(detailOutputFile);
        IOUtil.assertFileIsWritable(summaryOutputFile);

        // set up the collectors
        final Map<String, CollectVariantCallingMetrics.VariantCallingDetailMetrics> collapsedSampleDetailsMap = new HashMap<>();
        final CollectVariantCallingMetrics.VariantCallingSummaryMetrics collapsedSummary = new CollectVariantCallingMetrics.VariantCallingSummaryMetrics();

        for (final File file : INPUT) {
            final String inputPrefix = file.getAbsolutePath() + ".";

            try {
                // read in the detailed metrics file
                final File detail = new File(inputPrefix + CollectVariantCallingMetrics.VariantCallingDetailMetrics.getFileExtension());
                IOUtil.assertFileIsReadable(detail);
                MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, ?> detailedMetricsFile = getMetricsFile();
                detailedMetricsFile.read(new FileReader(detail));

                // for each sample in the detailed metrics...
                long totalHetDepth = 0L;
                for (final CollectVariantCallingMetrics.VariantCallingDetailMetrics detailedMetrics : detailedMetricsFile.getMetrics()) {
                    // re-calculate internal fields from derived fields
                    detailedMetrics.calculateFromDerivedFields();
                    totalHetDepth += detailedMetrics.TOTAL_HET_DEPTH;

                    // collapse it into to the main metrics for that sample
                    final CollectVariantCallingMetrics.VariantCallingDetailMetrics sampleDetails = collapsedSampleDetailsMap.computeIfAbsent(detailedMetrics.SAMPLE_ALIAS, f -> new CollectVariantCallingMetrics.VariantCallingDetailMetrics());
                    sampleDetails.merge(detailedMetrics);
                }

                // next, read in the new summary metrics
                final File summary = new File(inputPrefix + CollectVariantCallingMetrics.VariantCallingSummaryMetrics.getFileExtension());
                IOUtil.assertFileIsReadable(summary);
                MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
                summaryMetricsFile.read(new FileReader(summary));
                if (summaryMetricsFile.getMetrics().size() != 1) {
                    throw new PicardException(String.format("Expected 1 row in the summary metrics file but saw %d", summaryMetricsFile.getMetrics().size()));
                }

                // re-calculate internal fields from derived fields and merge it into the main summary metrics
                final CollectVariantCallingMetrics.VariantCallingSummaryMetrics summaryMetrics = summaryMetricsFile.getMetrics().get(0);
                summaryMetrics.calculateFromDerivedFields(totalHetDepth);
                collapsedSummary.merge(summaryMetrics);
            } catch (IOException e) {
                throw new PicardException(String.format("Cannot read from metrics files with prefix %s", inputPrefix));
            }
        }

        // prepare and write the finalized merged metrics
        final MetricsFile<CollectVariantCallingMetrics.VariantCallingDetailMetrics, Integer> detail = getMetricsFile();
        collapsedSampleDetailsMap.values().forEach((v)-> {
            v.calculateDerivedFields();
            detail.addMetric(v);
        });
        detail.write(detailOutputFile);

        final MetricsFile<CollectVariantCallingMetrics.VariantCallingSummaryMetrics, Integer> summary = getMetricsFile();
        collapsedSummary.calculateDerivedFields();
        summary.addMetric(collapsedSummary);
        summary.write(summaryOutputFile);

        return 0;
    }
}
