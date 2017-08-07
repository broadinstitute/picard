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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.metrics.GcBiasMetrics;
import picard.util.RExecutor;

import java.io.File;
import java.text.NumberFormat;
import java.util.List;
import java.util.Set;

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
        summary = CollectGcBiasMetrics.USAGE_SUMMARY + CollectGcBiasMetrics.USAGE_DETAILS,
        oneLineSummary = CollectGcBiasMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectGcBiasMetrics extends SinglePassSamProgram {

    static final String USAGE_SUMMARY = "Collect metrics regarding GC bias. ";
    static final String USAGE_DETAILS = "This tool collects information about the relative proportions of guanine (G) and cytosine (C)" +
            " nucleotides in a sample.  Regions of high and low G + C content have been shown to interfere with mapping/aligning," +
            " ultimately leading to fragmented genome assemblies and poor coverage in a phenomenon known as 'GC bias'.  " +
            "Detailed information on the effects of GC bias on the collection and analysis of sequencing data can be found at " +
            "DOI: 10.1371/journal.pone.0062856/.<br /><br />" +
            "" +
            "<p>The GC bias statistics are always output in a detailed long-form version, but a summary can also be produced. Both the " +
            "detailed metrics and the summary metrics are output as tables '.txt' files) and an accompanying chart that plots the " +
            "data ('.pdf' file). </p> " +
            "" +
            "<h4>Detailed metrics</h4>" +
            "The table of detailed metrics includes GC percentages for each bin (GC), the percentage of WINDOWS corresponding to each " +
            "GC bin of the reference sequence, the numbers of reads that start within a particular %GC content bin (READ_STARTS), and " +
            "the mean base quality of the reads that correspond to a specific GC content distribution window (MEAN_BASE_QUALITY).  " +
            "NORMALIZED_COVERAGE is a relative measure of sequence coverage by the reads at a particular GC content." +
            "" +
            "For each run, the corresponding reference sequence is divided into bins or windows based on the percentage of G + C" +
            " content ranging from 0 - 100%.  The percentages of G + C are determined from a defined length of sequence; the default " +
            "value is set at 100 bases. The mean of the distribution will vary among organisms; human DNA has a mean GC content " +
            "of 40%, suggesting a slight preponderance of AT-rich regions.  <br /><br />" +
            "" +
            "<h4>Summary metrics</h4>" +
            "The table of summary metrics captures run-specific bias information including WINDOW_SIZE, ALIGNED_READS, TOTAL_CLUSTERS, " +
            "AT_DROPOUT, and GC_DROPOUT.  While WINDOW_SIZE refers to the numbers of bases used for the distribution (see above), the " +
            "ALIGNED_READS and TOTAL_CLUSTERS are the total number of aligned reads and the total number of reads (after filtering) " +
            "produced in a run. In addition, the tool produces both AT_DROPOUT and GC_DROPOUT metrics, which indicate the percentage of " +
            "misaligned reads that correlate with low (%-GC is &lt; 50%) or high (%-GC is &gt; 50%) GC content respectively.  <br /><br />" +
            "" +
            "The percentage of 'coverage' or depth in a GC bin is calculated by dividing the number of reads of a particular GC content " +
            "by the mean number of reads of all GC bins.  A number of 1 represents mean coverage, a number less than 1 represents lower " +
            "than mean coverage (e.g. 0.5 means half as much coverage as average) while a number greater than 1 represents higher than " +
            "mean coverage (e.g. 3.1 means this GC bin has 3.1 times more reads per window than average).  " +
            "" +
            "This tool also tracks mean base-quality scores of the reads within each GC content bin, enabling the user to determine " +
            "how base quality scores vary with GC content.  <br /> <br />"+
            "" +
            "The chart output associated with this data table plots the NORMALIZED_COVERAGE, the distribution of WINDOWs corresponding " +
            "to GC percentages, and base qualities corresponding to each %GC bin."+
            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "<h4>Usage Example:</h4>"+
            "<pre>" +
            "java -jar picard.jar CollectGcBiasMetrics \\<br />"+
            "      I=input.bam \\<br />"+
            "      O=gc_bias_metrics.txt \\<br />"+
            "      CHART=gc_bias_metrics.pdf \\<br />"+
            "      S=summary_metrics.txt \\<br />"+
            "      R=reference_sequence.fasta"+
            "</pre>"+
            "Please see <a href='https://broadinstitute.github.io/picard/picard-metric-definitions.html#GcBiasMetrics'>" +
            "the GcBiasMetrics documentation</a> for further explanations of each metric." +
            "<hr />";
    /** The location of the R script to do the plotting. */
    private static final String R_SCRIPT = "picard/analysis/gcBias.R";

    // Usage and parameters

    @Argument(shortName = "CHART", doc = "The PDF file to render the chart to.")
    public File CHART_OUTPUT;

    @Argument(shortName = "S", doc = "The text file to write summary metrics to.")
    public File SUMMARY_OUTPUT;

    @Argument(shortName = "WINDOW_SIZE", doc = "The size of the scanning windows on the reference genome that are used to bin reads.")
    public int SCAN_WINDOW_SIZE = 100;

    @Argument(shortName = "MGF", doc = "For summary metrics, exclude GC windows that include less than this fraction of the genome.")
    public double MINIMUM_GENOME_FRACTION = 0.00001;

    @Argument(shortName = "BS", doc = "Whether the SAM or BAM file consists of bisulfite sequenced reads.")
    public boolean IS_BISULFITE_SEQUENCED = false;

    @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(shortName = "ALSO_IGNORE_DUPLICATES", doc = "Use to get additional results without duplicates. This option " +
            "allows to gain two plots per level at the same time: one is the usual one and the other excludes duplicates.")
    public boolean ALSO_IGNORE_DUPLICATES = false;

    // Calculates GcBiasMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private GcBiasMetricsCollector multiCollector;

    // Bins for the histograms to track the number of windows at each GC, and the number of read starts
    // at bins of each GC %. Need 101 to get from 0-100.
    private static final int BINS = 101;

    ////////////////////////////////////////////////////////////////////////////
    // Stock main method
    ////////////////////////////////////////////////////////////////////////////
    public static void main(final String[] args) {
        System.exit(new CollectGcBiasMetrics().instanceMain(args));
    }

    /////////////////////////////////////////////////////////////////////////////
    // Setup calculates windowsByGc for the entire reference. Must be done at
    // startup to avoid missing reference contigs in the case of small files
    // that may not have reads aligning to every reference contig.
    /////////////////////////////////////////////////////////////////////////////
    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        IOUtil.assertFileIsWritable(SUMMARY_OUTPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        //Calculate windowsByGc for the reference sequence
        final int[] windowsByGc = GcBiasUtils.calculateRefWindowsByGc(BINS, REFERENCE_SEQUENCE, SCAN_WINDOW_SIZE);

        //Delegate actual collection to GcBiasMetricCollector
        multiCollector = new GcBiasMetricsCollector(METRIC_ACCUMULATION_LEVEL, windowsByGc, header.getReadGroups(), SCAN_WINDOW_SIZE, IS_BISULFITE_SEQUENCED, ALSO_IGNORE_DUPLICATES);
    }

    ////////////////////////////////////////////////////////////////////////////
    // MultiCollector acceptRead
    ////////////////////////////////////////////////////////////////////////////
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        multiCollector.acceptRecord(rec, ref);
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


