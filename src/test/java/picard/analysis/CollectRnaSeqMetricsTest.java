/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;
import picard.annotation.RefFlatReader.RefFlatColumns;
import picard.util.RExecutor;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

public class CollectRnaSeqMetricsTest extends CommandLineProgramTest {
    public String getCommandLineProgramName() {
        return CollectRnaSeqMetrics.class.getSimpleName();
    }

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/RnaSeqMetrics");

    @Test
    public void testBasic() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addPair("pair2", sequenceIndex, 90, 225);
        builder.addPair("pair3", sequenceIndex, 120, 600);
        builder.addFrag("frag1", sequenceIndex, 150, true);
        builder.addFrag("frag2", sequenceIndex, 450, true);
        builder.addFrag("frag3", sequenceIndex, 225, false);
        builder.addPair("rrnaPair", sequenceIndex, 400, 500);

        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT=" +               samFile.getAbsolutePath(),
                "OUTPUT=" +              metricsFile.getAbsolutePath(),
                "REF_FLAT=" +            getRefFlatFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(metricsFile));

        final RnaSeqMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
        Assert.assertEquals(metrics.PF_BASES, 432);
        Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
        Assert.assertEquals(metrics.CODING_BASES, 136);
        Assert.assertEquals(metrics.UTR_BASES, 51);
        Assert.assertEquals(metrics.INTRONIC_BASES, 50);
        Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
        Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
        Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
        Assert.assertEquals(metrics.IGNORED_READS, 1);
        Assert.assertEquals(metrics.NUM_R1_TRANSCRIPT_STRAND_READS, 1);
        Assert.assertEquals(metrics.NUM_R2_TRANSCRIPT_STRAND_READS, 2);
        Assert.assertEquals(metrics.NUM_UNEXPLAINED_READS, 2);
        Assert.assertEquals(metrics.PCT_R1_TRANSCRIPT_STRAND_READS, 0.333333);
        Assert.assertEquals(metrics.PCT_R2_TRANSCRIPT_STRAND_READS, 0.666667);
    }

    @Test
    public void testChartFailureGATKLite() throws Exception {

        final PrintStream stderr = System.err;
        final String gatkLiteDockerProperty = System.getProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR);

        try {
            final ByteArrayOutputStream stderrCapture = new ByteArrayOutputStream();
            System.setErr(new PrintStream(stderrCapture));

            System.setProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR, "true");

            final String sequence = "chr1";
            final String ignoredSequence = "chrM";

            // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
            final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
            // Set seed so that strandedness is consistent among runs.
            builder.setRandomSeed(0);
            final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
            builder.addPair("pair1", sequenceIndex, 45, 475);
            builder.addPair("pair2", sequenceIndex, 90, 225);
            builder.addPair("pair3", sequenceIndex, 120, 600);
            builder.addFrag("frag1", sequenceIndex, 150, true);
            builder.addFrag("frag2", sequenceIndex, 450, true);
            builder.addFrag("frag3", sequenceIndex, 225, false);
            builder.addPair("rrnaPair", sequenceIndex, 400, 500);

            builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

            final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
            samFile.deleteOnExit();

            final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
            for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
            samWriter.close();

            // Create an interval list with one ribosomal interval.
            final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
            final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
            rRnaIntervalList.add(rRnaInterval);
            final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
            rRnaIntervalsFile.deleteOnExit();
            rRnaIntervalList.write(rRnaIntervalsFile);

            // Generate the metrics.
            final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
            metricsFile.deleteOnExit();

            final File chartFile = File.createTempFile("testChartFailureGATKLite", ".pdf");
            chartFile.deleteOnExit();

            final String[] args = new String[] {
                    "INPUT=" +               samFile.getAbsolutePath(),
                    "OUTPUT=" +              metricsFile.getAbsolutePath(),
                    "REF_FLAT=" +            getRefFlatFile(sequence).getAbsolutePath(),
                    "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                    "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                    "CHART_OUTPUT=" + chartFile.getAbsolutePath(),
                    "IGNORE_SEQUENCE=" + ignoredSequence
            };
            Assert.assertEquals(runPicardCommandLine(args), 1);

            Assert.assertTrue(stderrCapture.toString().contains("The histogram file cannot be written because it requires R, which is not available in the GATK Lite Docker image.")); 
        }
        finally {
            System.setErr(stderr);
            if(gatkLiteDockerProperty != null) {
                System.setProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR, gatkLiteDockerProperty);
            }
            else{
                System.clearProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR);
            } 
        }
    }

    @DataProvider(name = "rRnaIntervalsFiles")
    public static Object[][] rRnaIntervalsFiles() throws IOException {
        return new Object[][] {
                {null},
                {File.createTempFile("tmp.rRna.", ".interval_list")}
        };
    }

    @Test(dataProvider = "rRnaIntervalsFiles")
    public void testNoIntevalsNoFragPercentage(final File rRnaIntervalsFile) throws Exception {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);

        // Add a header but no intervals
        if ( rRnaIntervalsFile != null  ) {
            final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
            rRnaIntervalList.write(rRnaIntervalsFile);
            rRnaIntervalsFile.deleteOnExit();
        }

        // Create some alignments
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String rRnaIntervalsPath = rRnaIntervalsFile != null ? rRnaIntervalsFile.getAbsolutePath() : null;
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "REF_FLAT=" + getRefFlatFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsPath,
                "RRNA_FRAGMENT_PERCENTAGE=" + 0.0,
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence
        };
        try {
            Assert.assertEquals(runPicardCommandLine(args), 0);
        } catch(final Exception e) {
            if (rRnaIntervalsFile != null) {
                Assert.assertEquals(e.getCause().getClass(), PicardException.class);
            }
        }
    }

    @Test
    public void testMultiLevel() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, false);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);
        final SAMReadGroupRecord rg1 = new SAMReadGroupRecord("2");
        rg1.setSample("Sample");
        rg1.setLibrary("foo");
        builder.setReadGroup(rg1);
        builder.addPair("pair1", sequenceIndex, 45, 475);
        builder.addPair("pair2", sequenceIndex, 90, 225);
        builder.addFrag("frag1", sequenceIndex, 150, true);
        builder.addFrag("frag2", sequenceIndex, 450, true);

        final SAMReadGroupRecord rg2 = new SAMReadGroupRecord("3");
        rg2.setSample("Sample");
        rg2.setLibrary("bar");
        builder.setReadGroup(rg2);
        builder.addPair("pair3", sequenceIndex, 120, 600);
        builder.addFrag("frag3", sequenceIndex, 225, false);
        builder.addPair("rrnaPair", sequenceIndex, 400, 500);

        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT=" +               samFile.getAbsolutePath(),
                "OUTPUT=" +              metricsFile.getAbsolutePath(),
                "REF_FLAT=" +            getRefFlatFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence,
                "LEVEL=null",
                "LEVEL=SAMPLE",
                "LEVEL=LIBRARY"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
        output.read(new FileReader(metricsFile));

        for (final RnaSeqMetrics metrics : output.getMetrics()) {
            if (metrics.LIBRARY == null) {
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
                Assert.assertEquals(metrics.PF_BASES, 432);
                Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
                Assert.assertEquals(metrics.CODING_BASES, 136);
                Assert.assertEquals(metrics.UTR_BASES, 51);
                Assert.assertEquals(metrics.INTRONIC_BASES, 50);
                Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
                Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
                Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
                Assert.assertEquals(metrics.IGNORED_READS, 1);
                Assert.assertEquals(metrics.NUM_R1_TRANSCRIPT_STRAND_READS, 1);
                Assert.assertEquals(metrics.NUM_R2_TRANSCRIPT_STRAND_READS, 2);
                Assert.assertEquals(metrics.NUM_UNEXPLAINED_READS, 2);
                Assert.assertEquals(metrics.PCT_R1_TRANSCRIPT_STRAND_READS, 0.333333);
                Assert.assertEquals(metrics.PCT_R2_TRANSCRIPT_STRAND_READS, 0.666667);
            }
            else if (metrics.LIBRARY.equals("foo")) {
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 216);
                Assert.assertEquals(metrics.PF_BASES, 216);
                Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 36L);
                Assert.assertEquals(metrics.CODING_BASES, 89);
                Assert.assertEquals(metrics.UTR_BASES, 51);
                Assert.assertEquals(metrics.INTRONIC_BASES, 25);
                Assert.assertEquals(metrics.INTERGENIC_BASES, 15);
                Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
                Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 2);
                Assert.assertEquals(metrics.IGNORED_READS, 0);
                Assert.assertEquals(metrics.NUM_R1_TRANSCRIPT_STRAND_READS, 0);
                Assert.assertEquals(metrics.NUM_R2_TRANSCRIPT_STRAND_READS, 2);
                Assert.assertEquals(metrics.NUM_UNEXPLAINED_READS, 1);
                Assert.assertEquals(metrics.PCT_R1_TRANSCRIPT_STRAND_READS, 0.0);
                Assert.assertEquals(metrics.PCT_R2_TRANSCRIPT_STRAND_READS, 1.0);
            }
            else if (metrics.LIBRARY.equals("bar")) {
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 180);
                Assert.assertEquals(metrics.PF_BASES, 216);
                Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 72L);
                Assert.assertEquals(metrics.CODING_BASES, 47);
                Assert.assertEquals(metrics.UTR_BASES, 0);
                Assert.assertEquals(metrics.INTRONIC_BASES, 25);
                Assert.assertEquals(metrics.INTERGENIC_BASES, 36);
                Assert.assertEquals(metrics.CORRECT_STRAND_READS, 0);
                Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 2);
                Assert.assertEquals(metrics.IGNORED_READS, 1);
                Assert.assertEquals(metrics.NUM_R1_TRANSCRIPT_STRAND_READS, 1);
                Assert.assertEquals(metrics.NUM_R2_TRANSCRIPT_STRAND_READS, 0);
                Assert.assertEquals(metrics.NUM_UNEXPLAINED_READS, 1);
                Assert.assertEquals(metrics.PCT_R1_TRANSCRIPT_STRAND_READS, 1.0);
                Assert.assertEquals(metrics.PCT_R2_TRANSCRIPT_STRAND_READS, 0.0);
            }
        }
    }

    @Test
    public void testTranscriptionStrandMetrics() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        builder.setReadLength(50);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);

        // Reads that start *before* the gene: count as unexamined
        builder.addPair("pair_prior_1", sequenceIndex, 45, 50, false, false, "50M", "50M", true, false, -1);
        builder.addPair("pair_prior_2", sequenceIndex, 49, 50, false, false, "50M", "50M", true, false, -1);

        // Read that is enclosed in the gene, but one end does not map: count as unexamined
        builder.addPair("read_one_end_unmapped", sequenceIndex, 50, 51, false, true, "50M", null, false, false, -1);

        // Reads that are enclosed in the gene, paired and frag: one count per template
        builder.addFrag("frag_enclosed_forward", sequenceIndex, 150, false);
        builder.addFrag("frag_enclosed_reverse", sequenceIndex, 150, true);
        builder.addPair("pair_enclosed_forward", sequenceIndex, 150, 200, false, false, "50M", "50M", false, true, -1);
        builder.addPair("pair_enclosed_reverse", sequenceIndex, 200, 150, false, false, "50M", "50M", true, false, -1);

        // Read that starts *after* the gene: not counted
        builder.addPair("pair_after_1", sequenceIndex, 545, 550, false, false, "50M", "50M", true, false, -1);
        builder.addPair("pair_after_2", sequenceIndex, 549, 550, false, false, "50M", "50M", true, false, -1);

        // A read that uses the read length instead of the mate cigar
        builder.addFrag("frag_enclosed_forward_no_mate_cigar", sequenceIndex, 150, false).setAttribute(SAMTag.MC.name(), null);

        // Duplicate reads are counted
        builder.addFrag("frag_duplicate", sequenceIndex, 150, false).setDuplicateReadFlag(true);

        // These reads should be ignored.
        builder.addFrag("frag_secondary", sequenceIndex, 150, false).setNotPrimaryAlignmentFlag(true);
        builder.addFrag("frag_supplementary", sequenceIndex, 150, false).setSupplementaryAlignmentFlag(true);
        builder.addFrag("frag_qc_failure", sequenceIndex, 150, false).setReadFailsVendorQualityCheckFlag(true);

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT=" +               samFile.getAbsolutePath(),
                "OUTPUT=" +              metricsFile.getAbsolutePath(),
                "REF_FLAT=" +            getRefFlatFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
        output.read(new FileReader(metricsFile));

        final RnaSeqMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 900);
        Assert.assertEquals(metrics.PF_BASES, 950);
        Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 0L);
        Assert.assertEquals(metrics.CODING_BASES, 471);
        Assert.assertEquals(metrics.UTR_BASES, 125);
        Assert.assertEquals(metrics.INTRONIC_BASES, 98);
        Assert.assertEquals(metrics.INTERGENIC_BASES, 206);
        Assert.assertEquals(metrics.CORRECT_STRAND_READS, 7);
        Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 6);
        Assert.assertEquals(metrics.IGNORED_READS, 0);
        Assert.assertEquals(metrics.NUM_R1_TRANSCRIPT_STRAND_READS, 4);
        Assert.assertEquals(metrics.NUM_R2_TRANSCRIPT_STRAND_READS, 2);
        Assert.assertEquals(metrics.NUM_UNEXPLAINED_READS, 3);
        Assert.assertEquals(metrics.PCT_R1_TRANSCRIPT_STRAND_READS, 0.666667);
        Assert.assertEquals(metrics.PCT_R2_TRANSCRIPT_STRAND_READS, 0.333333);
    }

    @DataProvider(name = "endBiasValues")
    public static Object[][] endBiasValues() {
        return new Object[][] {
                // Case 1 - 25 nt in, ends should both be 0 coverage.
                {
                    25,
                    0.0,
                    0.0,
                    0.0
                },
                // Case 2 - 50 nt in, 5' should be 0 coverage and 3' should be partially covered.
                {
                    50,
                    1.543624,
                    0.0,
                    0.0
                },
                // Case 3 - 100 nt in, 3' and 5' should both have coverage.
                {
                    100,
                    2.449664,
                    1.677852,
                    0.684932
                }
        };
    }

    @Test(dataProvider = "endBiasValues")
    public void testBiasEndBiasAdjust(int end_bias_size, double expected_3prime, double expected_5prime, double expected_ratio) throws Exception {
        final String sequence = "chr1";

        // Create alignments that are shifted by a fixed distance from the 3' and 5' end of a transcript.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);

        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);

        final int read_length = 150;
        final int n_fragments = 30;
        final String read_cigar = read_length + "M";

        final int txStart = 0;
        final int txEnd = 1000;
        final int five_prime_start = txStart + 51;
        final int three_prime_start = txEnd - 25 - read_length;
        final String refFlatStart = Integer.toString(txStart);
        final String refFlatEnd = Integer.toString(txEnd);

        for (int j=0; j < n_fragments; j++) {
            builder.addPair("transcript_pair_" + j, sequenceIndex, five_prime_start, three_prime_start, false, false, read_cigar, read_cigar, false, true, -1);
        }

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        // Generate the metrics.
        final RnaSeqMetrics metrics;
        final MetricsFile<RnaSeqMetrics, Comparable<?>> output;

        final String refFlatFile = getRefFlatFile(sequence, refFlatStart, refFlatEnd, refFlatStart, refFlatEnd, "1", refFlatStart, refFlatEnd).getAbsolutePath();

        final String[] args = new String[] {
                "INPUT=" +               samFile.getAbsolutePath(),
                "OUTPUT=" +              metricsFile.getAbsolutePath(),
                "REF_FLAT=" +            refFlatFile,
                "STRAND_SPECIFICITY=" +     "NONE",
                "END_BIAS_BASES=" +     end_bias_size
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        output = new MetricsFile<>();
        output.read(new FileReader(metricsFile));
        metrics = output.getMetrics().get(0);

        Assert.assertEquals(metrics.MEDIAN_5PRIME_BIAS, expected_5prime);
        Assert.assertEquals(metrics.MEDIAN_3PRIME_BIAS, expected_3prime);
        Assert.assertEquals(metrics.MEDIAN_5PRIME_TO_3PRIME_BIAS, expected_ratio);
    }


    @Test
    public void testWithChimera() throws Exception {
        final File input = new File(TEST_DATA_DIR,"test_chimeras.sam");
        final File outfile   = File.createTempFile("test_with_chimeras", ".rna_seq_metrics");
        outfile.deleteOnExit();
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";
        
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "REF_FLAT=" + getRefFlatFile(sequence).getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "IGNORE_SEQUENCE=" + ignoredSequence,
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        Assert.assertTrue(outfile.exists());
    }
    

    @Test
    public void testPctRibosomalBases() throws Exception {
        final String sequence = "chr1";
        final String ignoredSequence = "chrM";

        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        // Set seed so that strandedness is consistent among runs.
        builder.setRandomSeed(0);
        final int sequenceIndex = builder.getHeader().getSequenceIndex(sequence);

        builder.addPair("rrnaPair", sequenceIndex, 400, 500, false, false, "35I1M", "35I1M", true, true, -1);
        builder.addFrag("frag1", sequenceIndex, 150, true, false, "34I2M", "*", -1);
        builder.addFrag("frag2", sequenceIndex, 190, true, false, "35I1M", "*", -1);
        builder.addFrag("ignoredFrag", builder.getHeader().getSequenceIndex(ignoredSequence), 1, false);
        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        try (SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile)) {
            for (final SAMRecord rec : builder.getRecords()) {
                samWriter.addAlignment(rec);
            }
        }

        // Create an interval list with one ribosomal interval.
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "REF_FLAT=" + getRefFlatFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(metricsFile));

        final RnaSeqMetrics metrics = output.getMetrics().get(0);

        Assert.assertEquals(metrics.PCT_RIBOSOMAL_BASES, 0.4);
    }

    public File getRefFlatFile(String sequence, String txStart, String txEnd, String cdsStart, String cdsEnd, String exonCount, String exonStarts, String exonEnds) throws Exception {
        // Create a refFlat file containing a single gene with specified parameters.

        final String[] refFlatFields = new String[RefFlatColumns.values().length];
        refFlatFields[RefFlatColumns.GENE_NAME.ordinal()] = "myGene";
        refFlatFields[RefFlatColumns.TRANSCRIPT_NAME.ordinal()] = "myTranscript";
        refFlatFields[RefFlatColumns.CHROMOSOME.ordinal()] = sequence;
        refFlatFields[RefFlatColumns.STRAND.ordinal()] = "+";
        refFlatFields[RefFlatColumns.TX_START.ordinal()] = txStart;
        refFlatFields[RefFlatColumns.TX_END.ordinal()] = txEnd;
        refFlatFields[RefFlatColumns.CDS_START.ordinal()] = cdsStart;
        refFlatFields[RefFlatColumns.CDS_END.ordinal()] = cdsEnd;
        refFlatFields[RefFlatColumns.EXON_COUNT.ordinal()] = exonCount;
        refFlatFields[RefFlatColumns.EXON_STARTS.ordinal()] = exonStarts;
        refFlatFields[RefFlatColumns.EXON_ENDS.ordinal()] = exonEnds;

        final File refFlatFile = File.createTempFile("tmp.", ".refFlat");
        refFlatFile.deleteOnExit();
        final PrintStream refFlatStream = new PrintStream(refFlatFile);
        refFlatStream.println(StringUtil.join("\t", refFlatFields));
        refFlatStream.close();

        return refFlatFile;
    }

    public File getRefFlatFile(String sequence) throws Exception {
        // Default refFlat file with a single gene containing two exons, one of which is overlapped by the
        // ribosomal interval.

        return getRefFlatFile(sequence, "49", "500", "74", "400", "2", "49,249", "200,500");
    }
}
