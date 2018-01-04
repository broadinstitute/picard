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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.tribble.annotation.Strand;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.annotation.Feature;
import picard.annotation.GtfRecord;
import picard.annotation.RefFlatRecord;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CollectRnaSeqMetricsTest extends CommandLineProgramTest {
    private final String sequence = "chr1";
    private final String ignoredSequence = "chrM";

    public String getCommandLineProgramName() {
        return CollectRnaSeqMetrics.class.getSimpleName();
    }

    @Test
    public void testBasicWithRefFlat() throws Exception {
        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = getSamRecords();

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final File rRnaIntervalsFile = createIntervals(builder);

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

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
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

    private File createIntervals(SAMRecordSetBuilder builder) throws IOException {
        final Interval rRnaInterval = new Interval(sequence, 300, 520, true, "rRNA");
        final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
        rRnaIntervalList.add(rRnaInterval);
        final File rRnaIntervalsFile = File.createTempFile("tmp.rRna.", ".interval_list");
        rRnaIntervalsFile.deleteOnExit();
        rRnaIntervalList.write(rRnaIntervalsFile);
        return rRnaIntervalsFile;
    }

    private SAMRecordSetBuilder getSamRecords() {
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
        return builder;
    }

    @Test
    public void testBasicWithGTF() throws Exception {
        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = getSamRecords();

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final File rRnaIntervalsFile = createIntervals(builder);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "GTF=" + getGTFFile(sequence).getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + rRnaIntervalsFile.getAbsolutePath(),
                "STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND",
                "IGNORE_SEQUENCE=" + ignoredSequence
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
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

    @DataProvider(name = "rRnaIntervalsFiles")
    public static Object[][] rRnaIntervalsFiles() throws IOException {
        return new Object[][] {
                {null},
                {File.createTempFile("tmp.rRna.", ".interval_list")}
        };
    }

    @Test(dataProvider = "rRnaIntervalsFiles")
    public void testNoIntervalsNoFragPercentageWithRefFlat(final File rRnaIntervalsFile) throws Exception {
        final File samFile = getSamRecordsWithConsistentStrandedness(rRnaIntervalsFile);

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
        } catch (final Exception e) {
            if (rRnaIntervalsFile != null) {
                Assert.assertEquals(e.getCause().getClass(), PicardException.class);
            }
        }
    }

    private File getSamRecordsWithConsistentStrandedness(File rRnaIntervalsFile) throws IOException {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);

        // Add a header but no intervals
        if ( rRnaIntervalsFile != null  ) {
            final IntervalList rRnaIntervalList = new IntervalList(builder.getHeader());
            rRnaIntervalList.write(rRnaIntervalsFile);
            rRnaIntervalsFile.deleteOnExit();
        }

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
        return samFile;
    }

    @Test(dataProvider = "rRnaIntervalsFiles")
    public void testNoIntervalsNoFragPercentageWithGTF(final File rRnaIntervalsFile) throws Exception {
        final File samFile = getSamRecordsWithConsistentStrandedness(rRnaIntervalsFile);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String rRnaIntervalsPath = rRnaIntervalsFile != null ? rRnaIntervalsFile.getAbsolutePath() : null;
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "GTF=" + getGTFFile(sequence).getAbsolutePath(),
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
    public void testMultiLevelWithRefFlat() throws Exception {
        final SAMRecordSetBuilder builder = getSamRecordsWithStrandednessReads();


        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final File rRnaIntervalsFile = createIntervals(builder);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "REF_FLAT=" + getRefFlatFile(sequence).getAbsolutePath(),
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
            } else if (metrics.LIBRARY.equals("foo")) {
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
            } else if (metrics.LIBRARY.equals("bar")) {
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

    private SAMRecordSetBuilder getSamRecordsWithStrandednessReads() {
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
        return builder;
    }

    @Test
    public void testMultiLevelWithGTF() throws Exception {
        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = getSamRecordsWithStrandednessReads();

        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec: builder.getRecords()) samWriter.addAlignment(rec);
        samWriter.close();

        // Create an interval list with one ribosomal interval.
        final File rRnaIntervalsFile = createIntervals(builder);

        // Generate the metrics.
        final File metricsFile = File.createTempFile("tmp.", ".rna_metrics");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "GTF=" + getGTFFile(sequence).getAbsolutePath(),
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
    public void testTranscriptionStrandMetricsWithRefFlat() throws Exception {
        final SAMRecordSetBuilder builder = getSamRecordsWithDifferentTypesOfReads();


        final File samFile = File.createTempFile("tmp.collectRnaSeqMetrics.", ".sam");
        samFile.deleteOnExit();

        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(builder.getHeader(), false, samFile);
        for (final SAMRecord rec : builder.getRecords()) samWriter.addAlignment(rec);
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

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "REF_FLAT=" + getRefFlatFile(sequence).getAbsolutePath(),
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

    private SAMRecordSetBuilder getSamRecordsWithDifferentTypesOfReads() {
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
        return builder;
    }

    @Test
    public void testTranscriptionStrandMetricsWithGTF() throws Exception {
        // Create some alignments that hit the ribosomal sequence, various parts of the gene, and intergenic.
        final SAMRecordSetBuilder builder = getSamRecordsWithDifferentTypesOfReads();

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

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "GTF=" + getGTFFile(sequence).getAbsolutePath(),
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

    private File getRefFlatFile(String sequence) throws Exception {
        // Create a refFlat file with a single gene containing two exons, one of which is overlapped by the
        // ribosomal interval.
        final RefFlatRecord record = new RefFlatRecord("myGene", "myTranscript", sequence, Strand.POSITIVE, 49, 500, 74, 400, 2, new int[]{49, 249}, new int[]{200, 500});

        final File file = File.createTempFile("tmp.", ".refFlat");
        file.deleteOnExit();
        final PrintStream filePrintStream = new PrintStream(file);
        filePrintStream.println(record.toRow());
        filePrintStream.close();

        return file;
    }

    private File getGTFFile(String sequence) throws Exception {
        // Create a GTF file with a single gene containing two exons

        String source = "mySource";
        Map<String, String> attributes = new HashMap<>();
        attributes.put("gene_id", "myGene");
        attributes.put("transcript_id", "myTranscript");
        Strand strand = Strand.POSITIVE;

        final List<GtfRecord> records = Arrays.asList(
                new GtfRecord(sequence, source, Feature.START_CODON, 75, 77, ".", strand, ".", attributes),
                new GtfRecord(sequence, source, Feature.CDS, 75, 200, ".", strand, ".", attributes),
                new GtfRecord(sequence, source, Feature.EXON, 50, 200, ".", strand, ".", attributes),
                new GtfRecord(sequence, source, Feature.CDS, 250, 397, ".", strand, ".", attributes),
                new GtfRecord(sequence, source, Feature.STOP_CODON, 398, 400, ".", strand, ".", attributes),
                new GtfRecord(sequence, source, Feature.EXON, 250, 500, ".", strand, ".", attributes)
        );

        final File file = File.createTempFile("tmp.", ".gtf");
        file.deleteOnExit();
        final PrintStream fileStream = new PrintStream(file);
        records.forEach(record -> fileStream.println(record.toRow()));

        fileStream.close();

        return file;
    }

}

