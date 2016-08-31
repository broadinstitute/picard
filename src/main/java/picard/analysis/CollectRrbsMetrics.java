/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Calculates and reports QC metrics for RRBS data based on the methylation status at individual C/G bases as well
 * as CpG sites across all reads in the input BAM/SAM file.
 *
 * @author jgentry@broadinstitute.org
 */

@CommandLineProgramProperties(
        summary = CollectRrbsMetrics.USAGE_SUMMARY + CollectRrbsMetrics.USAGE_DETAILS,
        oneLineSummary = CollectRrbsMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
@DocumentedFeature
public class CollectRrbsMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "<b>Collects metrics from reduced representation bisulfite sequencing (Rrbs) data.</b>  ";
    static final String USAGE_DETAILS = "<p>This tool uses reduced representation bisulfite sequencing (Rrbs) data to determine cytosine " +
            "methylation status across all reads of a genomic DNA sequence.  For a primer on bisulfite sequencing and cytosine methylation, " +
            "please see the corresponding <a href='https://www.broadinstitute.org/gatk/guide/article?id=6330'>GATK Dictionary entry</a>. </p>" +

            "<p>Briefly, bisulfite reduction converts un-methylated cytosine (C) to uracil (U) bases.  Methylated sites are not converted " +
            "because they are resistant to bisulfite reduction.  Subsequent to PCR amplification of the reaction products, bisulfite " +
            "reduction manifests as [C -> T (+ strand) or G -> A (- strand)] base conversions.  Thus, conversion rates" +
            " can be calculated from the reads as follows: [CR = converted/(converted + unconverted)]. Since methylated cytosines are " +
            "protected against Rrbs-mediated conversion, the methylation rate (MR) is as follows:" +
            "[MR = unconverted/(converted + unconverted) = (1 - CR)].</p>" +

            "<p>The CpG CollectRrbsMetrics tool outputs three files including summary and detail metrics tables as well as a PDF file containing " +
            "four graphs. These graphs are derived from the summary table and include a comparison of conversion rates for both CpG and non-CpG sites, " +
            "the distribution of total numbers of CpG sites as a function of the CpG conversion rates, the distribution of CpG sites by the level of " +
            "read coverage (depth), and the numbers of reads discarded resulting from either exceeding the mismatch rate or size (too short).  " +
            "The detailed metrics table includes the coordinates of all of the CpG sites for the experiment as well as the conversion rates " +
            "observed for each site.</p>" +

            "<p>It is possible to launch CollectRrbsMetrics tool by setting all output files or by indicating a base name for metrics files only. See examples:</p>" +

            "<h4>Usage example 1:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectRrbsMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=detail_output.bam \\<br />" +
            "      CHART=chart_output.pdf \\<br />" +
            "      S=summary_output.bam \\<br />" +
            "      R=reference_sequence.fasta" +
            "</pre>" +
            "<h4>Usage example 2:</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectRrbsMetrics \\<br />" +
            "      R=reference_sequence.fasta \\<br />" +
            "      I=input.bam \\<br />" +
            "      M=basename_for_metrics_files" +
            "</pre>" +

            "<p>Please see the CollectRrbsMetrics " +
            "<a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#RrbsCpgDetailMetrics'>definitions</a>" +
            " for a complete description of both the detail and summary metrics produced by this tool.</p>" +
            "<hr />";

    // Path to R file for plotting purposes

    private static final String R_SCRIPT = "picard/analysis/rrbsQc.R";

    @Argument(doc = "Base name for output files", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, mutex = {"OUTPUT", "CHART_OUTPUT", "SUMMARY_OUTPUT"})
    public String METRICS_FILE_PREFIX;
    @Argument(doc = "Minimum read length")
    public int MINIMUM_READ_LENGTH = 5;
    @Argument(doc = "Threshold for base quality of a C base before it is considered")
    public int C_QUALITY_THRESHOLD = 20;
    @Argument(doc = "Threshold for quality of a base next to a C before the C base is considered")
    public int NEXT_BASE_QUALITY_THRESHOLD = 10;
    @Argument(doc = "Maximum percentage of mismatches in a read for it to be considered, with a range of 0-1")
    public double MAX_MISMATCH_RATE = 0.1;
    @Argument(doc = "Set of sequence names to consider, if not specified all sequences will be used", optional = true)
    public Set<String> SEQUENCE_NAMES = new HashSet<>();
    @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);
    @Argument(shortName = "CHART", doc = "The PDF file to render the chart to.", mutex = {"METRICS_FILE_PREFIX"})
    public File CHART_OUTPUT;
    @Argument(shortName = "S", doc = "The text file to write summary metrics to.", mutex = {"METRICS_FILE_PREFIX"})
    public File SUMMARY_OUTPUT;

    public static final String DETAIL_FILE_EXTENSION = "rrbs_detail_metrics";
    public static final String SUMMARY_FILE_EXTENSION = "rrbs_summary_metrics";
    public static final String PDF_FILE_EXTENSION = "rrbs_qc.pdf";

    private RrbsMetricsCollector metricsCollector;

    // return a custom argument collection since this tool uses a (required) argument name
    // of "REFERENCE", not "REFERENCE_SEQUENCE"
    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new CollectRrbsMetricsReferenceArgumentCollection();
    }

    public static class CollectRrbsMetricsReferenceArgumentCollection implements ReferenceArgumentCollection {
        @Argument(doc = "The reference sequence fasta file", shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME)
        public File REFERENCE;

        @Override
        public File getReferenceFile() {
            return REFERENCE;
        };
    }

    public static void main(final String[] args) {
        new CollectRrbsMetrics().instanceMainWithExit(args);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        if (METRICS_FILE_PREFIX != null) {
            if (!METRICS_FILE_PREFIX.endsWith(".")) {
                METRICS_FILE_PREFIX = METRICS_FILE_PREFIX + ".";
            }
            OUTPUT = new File(METRICS_FILE_PREFIX + DETAIL_FILE_EXTENSION);
            SUMMARY_OUTPUT = new File(METRICS_FILE_PREFIX + SUMMARY_FILE_EXTENSION);
            CHART_OUTPUT = new File(METRICS_FILE_PREFIX + PDF_FILE_EXTENSION);
        }
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        metricsCollector = new RrbsMetricsCollector(
                METRIC_ACCUMULATION_LEVEL,
                header.getReadGroups(),
                C_QUALITY_THRESHOLD,
                NEXT_BASE_QUALITY_THRESHOLD,
                MINIMUM_READ_LENGTH,
                MAX_MISMATCH_RATE
        );
    }

    @Override
    protected void acceptRead(final SAMRecord samRecord, final ReferenceSequence ref) {
        if (!samRecord.getReadUnmappedFlag() && !isSequenceFiltered(samRecord.getReferenceName())) {
            metricsCollector.acceptRecord(samRecord, ref);
        }
    }

    private boolean isSequenceFiltered(final String sequenceName) {
        return SEQUENCE_NAMES != null
                && !SEQUENCE_NAMES.isEmpty()
                && !SEQUENCE_NAMES.contains(sequenceName);
    }

    @Override
    protected void finish() {
        metricsCollector.finish();
        final MetricsFile<RrbsMetrics, Comparable<?>> rrbsMetrics = getMetricsFile();
        metricsCollector.addAllLevelsToFile(rrbsMetrics);

        // Using RrbsMetrics as a way to get both of the metrics objects through the MultiLevelCollector. Once
        // we get it out split it apart to the two separate MetricsFiles and write them to file
        final MetricsFile<RrbsSummaryMetrics, ?> summaryFile = getMetricsFile();
        final MetricsFile<RrbsCpgDetailMetrics, ?> detailsFile = getMetricsFile();
        rrbsMetrics.getMetrics().forEach(rrbsMetric -> {
            summaryFile.addMetric(rrbsMetric.getSummaryMetrics());
            rrbsMetric.getDetailMetrics().forEach(detailsFile::addMetric);
        });

        summaryFile.write(SUMMARY_OUTPUT);
        detailsFile.write(OUTPUT);

        RExecutor.executeFromClasspath(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                SUMMARY_OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath());
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errorMsgs = new ArrayList<>();
        if (MAX_MISMATCH_RATE < 0 || MAX_MISMATCH_RATE > 1) {
            errorMsgs.add("MAX_MISMATCH_RATE must be in the range of 0-1");
        }

        if (C_QUALITY_THRESHOLD < 0) {
            errorMsgs.add("C_QUALITY_THRESHOLD must be >= 0");
        }

        if (NEXT_BASE_QUALITY_THRESHOLD < 0) {
            errorMsgs.add("NEXT_BASE_QUALITY_THRESHOLD must be >= 0");
        }

        if (MINIMUM_READ_LENGTH <= 0) {
            errorMsgs.add("MINIMUM_READ_LENGTH must be > 0");
        }

        return errorMsgs.isEmpty() ? null : errorMsgs.toArray(new String[errorMsgs.size()]);
    }
}
