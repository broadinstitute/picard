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

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.metrics.GcBiasMetrics;
import picard.util.RExecutor;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.List;

/**
 * Tool to collect information about GC bias in the reads in a given BAM file. Computes
 * the number of windows (of size specified by SCAN_WINDOW_SIZE) in the genome at each GC%
 * and counts the number of read starts in each GC bin.  What is output and plotted is
 * the "normalized coverage" in each bin - i.e. the number of reads per window normalized
 * to the average number of reads per window across the whole genome.
 *
 * @author Tim Fennell
 * edited by Kylee Bergin
 */
@CommandLineProgramProperties(
        summary = picard.analysis.CollectGcBiasMetrics.USAGE_SUMMARY + picard.analysis.CollectGcBiasMetrics.USAGE_DETAILS,
        oneLineSummary = picard.analysis.CollectGcBiasMetrics.USAGE_SUMMARY,
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class EstimateFlowCoverage extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Collect metrics regarding Flow. ";
    static final String USAGE_DETAILS = "";
    /**
     * The location of the R script to do the plotting.
     */
    private static final String R_SCRIPT = "picard/analysis/flowBias.R";

    // Usage and parameters

    @Argument(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Argument(shortName = "S", doc = "The text file to write summary metrics to.")
    public File SUMMARY_OUTPUT;

    @Argument(shortName = "FLOW", doc = "The flow order used in sequencing (will be used repeatedly)")
    public String FLOW_ORDER;

    @Argument(shortName = "SIZE", doc = "The number of flows that were used")
    public Integer FLOW_LENGTH;

    /////////////////////////////////////////////////////////////////////////////
    // Setup calculates windowsByGc for the entire reference. Must be done at
    // startup to avoid missing reference contigs in the case of small files
    // that may not have reads aligning to every reference contig.
    /////////////////////////////////////////////////////////////////////////////

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        return super.customCommandLineValidation();
    }

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {
        try (ReferenceSequenceFileWalker referenceSequenceWalker = new ReferenceSequenceFileWalker(OUTPUT)) {
            for (SAMSequenceRecord samSequenceRecord : referenceSequenceWalker.getSequenceDictionary().getSequences()) {
                final ReferenceSequence referenceSequence = referenceSequenceWalker.get(samSequenceRecord.getSequenceIndex());

                updateMetrics(referenceSequence);

            }

        } catch (IOException e) {
            throw new PicardException("Error while reading reference " + REFERENCE_SEQUENCE, e);
        }
        return 0;
    }

    private void updateMetrics(final ReferenceSequence referenceSequence) {
        .....
    }

    /////////////////////////////////////////////////////////////////////////////
    // Write out all levels of normalized coverage metrics to a file
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void finish() {
        multiCollector.finish();
        writeResultsToFiles();
    }

    private void writeResultsToFiles() {
        final MetricsFile<GcBiasMetrics, Integer> file = getMetricsFile();
        final MetricsFile<GcBiasDetailMetrics, ?> detailMetricsFile = getMetricsFile();
        final MetricsFile<GcBiasSummaryMetrics, ?> summaryMetricsFile = getMetricsFile();
        multiCollector.addAllLevelsToFile(file);
        final List<GcBiasMetrics> gcBiasMetricsList = file.getMetrics();
        for (final GcBiasMetrics gcbm : gcBiasMetricsList) {
            final List<GcBiasDetailMetrics> gcDetailList = gcbm.DETAILS.getMetrics();
            for (final GcBiasDetailMetrics d : gcDetailList) {
                detailMetricsFile.addMetric(d);
            }
            summaryMetricsFile.addMetric(gcbm.SUMMARY);
        }
        detailMetricsFile.write(OUTPUT);
        summaryMetricsFile.write(SUMMARY_OUTPUT);

        final NumberFormat fmt = NumberFormat.getIntegerInstance();
        fmt.setGroupingUsed(true);
        RExecutor.executeFromClasspath(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                SUMMARY_OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                String.valueOf(SCAN_WINDOW_SIZE));
    }
}


