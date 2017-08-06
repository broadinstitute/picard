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
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;



@CommandLineProgramProperties(
        summary = CollectBaseDistributionByCycle.USAGE_SUMMARY + CollectBaseDistributionByCycle.USAGE_DETAILS,
        oneLineSummary = CollectBaseDistributionByCycle.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectBaseDistributionByCycle extends SinglePassSamProgram {
        static final String USAGE_SUMMARY = "Chart the nucleotide distribution per cycle in a SAM or BAM file";
        static final String USAGE_DETAILS = "This tool produces a chart of the nucleotide distribution per cycle in a SAM or BAM file " +
                "in order to enable assessment of systematic errors at specific positions in the reads.<br /><br />" +
                "" +
                "<h4>Interpretation notes</h4>" +
                "Increased numbers of miscalled bases will be reflected in base distribution changes and increases in the number of Ns.  "+
                "In general, we expect that for any given cycle, or position within reads, the relative proportions of A, T, C and G "+
                "should reflect the AT:GC content of the organism's genome.  Thus, for all four nucleotides, flattish lines would be "+
                "expected.  Deviations from this expectation, for example a spike of A at a particular cycle (position within reads), "+
                "would suggest a systematic sequencing error."+
                "" +
                "<h4>Note on quality trimming</h4>" +
                "In the past, many sequencing data processing workflows included discarding the low-quality tails of reads by applying "+
                "hard-clipping at some arbitrary base quality threshold value. This is no longer useful because most sophisticated "+
                "analysis tools (such as the GATK variant discovery tools) are quality-aware, meaning that they are able to take base "+
                "quality into account when weighing evidence provided by sequencing reads. Unnecessary clipping may interfere with other "+
                "quality control evaluations and may lower the quality of analysis results. For example, trimming reduces the effectiveness "+
                "of the Base Recalibration (BQSR) pre-processing step of the "+
                "<a href='https://www.broadinstitute.org/gatk/guide/best-practices'>GATK Best Practices for Variant Discovery</a>, "+
                "which aims to correct some types of systematic biases that affect the accuracy of base quality scores."+
                
                "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>"+
                
                "<h4>Usage example:</h4>" +
                "<pre>" +
                "java -jar picard.jar CollectBaseDistributionByCycle \\<br />" +
                "      CHART=collect_base_dist_by_cycle.pdf \\<br />" +
                "      I=input.bam \\<br />" +
                "      O=output.txt" +
                "</pre>" +
                "<hr />"
                ;
    @Argument(shortName = "CHART", doc = "A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Argument(doc = "If set to true, calculate the base distribution over aligned reads only.")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(doc = "If set to true, calculate the base distribution over PF reads only (Illumina specific). PF reads are reads that passed the internal quality filters applied by Illumina sequencers.")
    public boolean PF_READS_ONLY = false;

    private HistogramGenerator hist;
    private String plotSubtitle = "";
    private final Log log = Log.getInstance(CollectBaseDistributionByCycle.class);

    public static void main(String[] args) {
        System.exit(new CollectBaseDistributionByCycle().instanceMain(args));
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(CHART_OUTPUT);
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            plotSubtitle = StringUtil.asEmptyIfNull(readGroups.get(0).getLibrary());
        }
        hist = new HistogramGenerator();
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        if ((PF_READS_ONLY) && (rec.getReadFailsVendorQualityCheckFlag())) {
            return;
        }
        if ((ALIGNED_READS_ONLY) && (rec.getReadUnmappedFlag())) {
            return;
        }
        if (rec.isSecondaryOrSupplementary()) {
            return;
        }
        hist.addRecord(rec);
    }

    @Override
    protected void finish() {
        final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics = getMetricsFile();
        hist.addToMetricsFile(metrics);
        metrics.write(OUTPUT);
        if (hist.isEmpty()) {
            log.warn("No valid bases found in input file. No plot will be produced.");
        } else {
            final int rResult = RExecutor.executeFromClasspath("picard/analysis/baseDistributionByCycle.R",
                    OUTPUT.getAbsolutePath(),
                    CHART_OUTPUT.getAbsolutePath(),
                    INPUT.getName(),
                    plotSubtitle);
            if (rResult != 0) {
                throw new PicardException("R script nucleotideDistributionByCycle.R failed with return code " + rResult);
            }
        }
    }

    private class HistogramGenerator {
        private int maxLengthSoFar = 0;
        final private long[][] firstReadTotalsByCycle = new long[5][maxLengthSoFar];
        private long[] firstReadCountsByCycle = new long[maxLengthSoFar];
        final private long[][] secondReadTotalsByCycle = new long[5][maxLengthSoFar];
        private long[] secondReadCountsByCycle = new long[maxLengthSoFar];
        private boolean seenSecondEnd = false;

        private int baseToInt(final byte base) {
            switch (base) {
                case 'A':
                case 'a':
                    return 0;
                case 'C':
                case 'c':
                    return 1;
                case 'G':
                case 'g':
                    return 2;
                case 'T':
                case 't':
                    return 3;
            }
            return 4;
        }

        void addRecord(final SAMRecord rec) {
            final byte[] bases = rec.getReadBases();
            if (bases == null) {
                return;
            }
            final int length = bases.length;
            final boolean rc = rec.getReadNegativeStrandFlag();
            ensureArraysBigEnough(length + 1);
            if ((rec.getReadPairedFlag()) && (rec.getSecondOfPairFlag())) {
                seenSecondEnd = true;
                for (int i = 0; i < length; i++) {
                    final int cycle = rc ? length - i : i + 1;
                    secondReadTotalsByCycle[baseToInt(bases[i])][cycle] += 1;
                    secondReadCountsByCycle[cycle] += 1;
                }
            } else {
                for (int i = 0; i < length; i++) {
                    final int cycle = rc ? length - i : i + 1;
                    firstReadTotalsByCycle[baseToInt(bases[i])][cycle] += 1;
                    firstReadCountsByCycle[cycle] += 1;
                }
            }
        }

        private void ensureArraysBigEnough(final int length) {
            if (length > maxLengthSoFar) {
                for (int i = 0; i < 5; i++) {
                    firstReadTotalsByCycle[i] = Arrays.copyOf(firstReadTotalsByCycle[i], length);
                    secondReadTotalsByCycle[i] = Arrays.copyOf(secondReadTotalsByCycle[i], length);
                }
                firstReadCountsByCycle = Arrays.copyOf(firstReadCountsByCycle, length);
                secondReadCountsByCycle = Arrays.copyOf(secondReadCountsByCycle, length);
                maxLengthSoFar = length;
            }
        }

        boolean isEmpty() {
            return maxLengthSoFar == 0;
        }

        public void addToMetricsFile(final MetricsFile<BaseDistributionByCycleMetrics, ?> metrics) {
            int firstReadLength = 0;
            for (int i = 0; i < maxLengthSoFar; i++) {
                if (0 != firstReadCountsByCycle[i]) {
                    final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                    metric.READ_END = 1;
                    metric.CYCLE = i;
                    metric.PCT_A = (100.0 * firstReadTotalsByCycle[0][i] / firstReadCountsByCycle[i]);
                    metric.PCT_C = (100.0 * firstReadTotalsByCycle[1][i] / firstReadCountsByCycle[i]);
                    metric.PCT_G = (100.0 * firstReadTotalsByCycle[2][i] / firstReadCountsByCycle[i]);
                    metric.PCT_T = (100.0 * firstReadTotalsByCycle[3][i] / firstReadCountsByCycle[i]);
                    metric.PCT_N = (100.0 * firstReadTotalsByCycle[4][i] / firstReadCountsByCycle[i]);
                    metrics.addMetric(metric);
                    firstReadLength = i;
                }
            }
            if (seenSecondEnd) {
                for (int i = 0; i < maxLengthSoFar; i++) {
                    if (0 != secondReadCountsByCycle[i]) {
                        final BaseDistributionByCycleMetrics metric = new BaseDistributionByCycleMetrics();
                        metric.READ_END = 2;
                        metric.CYCLE = (i + firstReadLength);
                        metric.PCT_A = (100.0 * secondReadTotalsByCycle[0][i] / secondReadCountsByCycle[i]);
                        metric.PCT_C = (100.0 * secondReadTotalsByCycle[1][i] / secondReadCountsByCycle[i]);
                        metric.PCT_G = (100.0 * secondReadTotalsByCycle[2][i] / secondReadCountsByCycle[i]);
                        metric.PCT_T = (100.0 * secondReadTotalsByCycle[3][i] / secondReadCountsByCycle[i]);
                        metric.PCT_N = (100.0 * secondReadTotalsByCycle[4][i] / secondReadCountsByCycle[i]);
                        metrics.addMetric(metric);
                    }
                }
            }
        }
    }
}
