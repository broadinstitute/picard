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

package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.util.TestNGUtil;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;

/**
 * Tests CollectAlignmentSummaryStatistics
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class CollectAlignmentSummaryMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    public String getCommandLineProgramName() {
        return CollectAlignmentSummaryMetrics.class.getSimpleName();
    }
    
    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(), 3);
        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
                case FIRST_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 7);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 3);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 59);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 19.0);
                    Assert.assertEquals(metrics.PF_READS_ALIGNED, 3);
                    Assert.assertEquals(metrics.PF_READS_IMPROPER_PAIRS, 1);
                    Assert.assertEquals(metrics.PCT_PF_READS_IMPROPER_PAIRS, 0.333333 /* 1/3 */);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 303);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*58D/303D*/0.191419);
                    Assert.assertEquals(metrics.BAD_CYCLES, 19);
                    break;
                case SECOND_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 9);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 7);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 239);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                    Assert.assertEquals(metrics.PF_READS_ALIGNED, 7);
                    Assert.assertEquals(metrics.PF_READS_IMPROPER_PAIRS, 5);
                    Assert.assertEquals(metrics.PCT_PF_READS_IMPROPER_PAIRS, 0.714286 /* 5/7 */);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 707);
                    Assert.assertEquals(metrics.PCT_READS_ALIGNED_IN_PAIRS, 0.285714 /* 2D/7 */);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*19D/707D*/0.026874);
                    Assert.assertEquals(metrics.BAD_CYCLES, 3);
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 18);
                    Assert.assertEquals(metrics.PF_READS, 16);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 10);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 298);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                    Assert.assertEquals(metrics.PF_READS_ALIGNED, 10);
                    Assert.assertEquals(metrics.PF_READS_IMPROPER_PAIRS, 6);
                    Assert.assertEquals(metrics.PCT_PF_READS_IMPROPER_PAIRS, 0.6 /* 6/10 */);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 1010);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*77D/1010D*/0.076238);
                    Assert.assertEquals(metrics.BAD_CYCLES, 22);
                    break;
                case UNPAIRED:
                default:
                    Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testBisulfite() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_bisulfite_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "IS_BISULFITE_SEQUENCED=true"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final NumberFormat format =  NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                // 19 no-calls, one potentially methylated base, one mismatch at a potentially methylated base
                Assert.assertEquals(metrics.TOTAL_READS, 1);
                Assert.assertEquals(metrics.PF_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 21.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.212121 /*21D/99D*/);
                Assert.assertEquals(metrics.BAD_CYCLES, 21);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(21/(double)99));
                break;
            case SECOND_OF_PAIR:
                // Three no-calls, two potentially methylated bases
                Assert.assertEquals(metrics.TOTAL_READS, 1);
                Assert.assertEquals(metrics.PF_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 4.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*4D/99D*/0.040404);
                Assert.assertEquals(metrics.BAD_CYCLES, 4);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(4/(double)99));
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 2);
                Assert.assertEquals(metrics.PF_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 202);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 12.5D);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 202);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.126263);// 25D/198D
                Assert.assertEquals(metrics.BAD_CYCLES, 25);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(25/(double)198));
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testBisulfiteButNot() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_bisulfite_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "IS_BISULFITE_SEQUENCED=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final NumberFormat format =  NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
                case FIRST_OF_PAIR:
                    // 19 no-calls, one potentially methylated base, one mismatch at a potentially methylated base
                    Assert.assertEquals(metrics.TOTAL_READS, 1);
                    Assert.assertEquals(metrics.PF_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 23.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 101);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.227723 /*23D/101D*/);
                    Assert.assertEquals(metrics.BAD_CYCLES, 23);
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(23/(double)101));
                    break;
                case SECOND_OF_PAIR:
                    // Three no-calls, two potentially methylated bases
                    Assert.assertEquals(metrics.TOTAL_READS, 1);
                    Assert.assertEquals(metrics.PF_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 6.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 101);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*6D/101D*/0.059406);
                    Assert.assertEquals(metrics.BAD_CYCLES, 6);
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(6/(double)101));
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 2);
                    Assert.assertEquals(metrics.PF_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 202);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 14.5D);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 202);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.143564);// 29D/202D
                    Assert.assertEquals(metrics.BAD_CYCLES, 29);
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(29/(double)202));
                    break;
                case UNPAIRED:
                default:
                    Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testNoReference() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 7);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                break;
            case SECOND_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 9);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 18);
                Assert.assertEquals(metrics.PF_READS, 16);
                Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testZeroLengthReads() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test2.sam");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));
        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            // test that it doesn't blow up
        }
    }

    @Test
    public void testMultipleLevelsOfMetrics() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test_multiple.sam");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=ALL_READS",
                "METRIC_ACCUMULATION_LEVEL=SAMPLE",
                "METRIC_ACCUMULATION_LEVEL=LIBRARY",
                "METRIC_ACCUMULATION_LEVEL=READ_GROUP",
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            if (metrics.SAMPLE == null) {
                switch (metrics.CATEGORY) {
                case FIRST_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 7);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                    Assert.assertEquals(metrics.BAD_CYCLES, 19);
                    break;
                case SECOND_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 9);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                    Assert.assertEquals(metrics.BAD_CYCLES, 3);
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 18);
                    Assert.assertEquals(metrics.PF_READS, 16);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                    Assert.assertEquals(metrics.BAD_CYCLES, 22);
                    break;
                case UNPAIRED:
                default:
                    Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                }
            }
            else if (metrics.SAMPLE.equals("Ma")) {
                // There's only one library and one read group for this sample so the metrics for
                // every level should be identical
                switch (metrics.CATEGORY) {
                    case FIRST_OF_PAIR:
                        Assert.assertEquals(metrics.TOTAL_READS, 5);
                        Assert.assertEquals(metrics.PF_READS, 3);
                        Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                        Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                        Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                        Assert.assertEquals(metrics.BAD_CYCLES, 24);
                        break;
                    case SECOND_OF_PAIR:
                        Assert.assertEquals(metrics.TOTAL_READS, 5);
                        Assert.assertEquals(metrics.PF_READS, 5);
                        Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                        Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                        Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                        Assert.assertEquals(metrics.BAD_CYCLES, 3);
                        break;
                    case PAIR:
                        Assert.assertEquals(metrics.TOTAL_READS, 10);
                        Assert.assertEquals(metrics.PF_READS, 8);
                        Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                        Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                        Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                        Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                        Assert.assertEquals(metrics.BAD_CYCLES, 27);
                        break;
                    case UNPAIRED:
                    default:
                        Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                }
            }
            else if (metrics.SAMPLE.equals("Pa")) {
                // Two libraries and three read groups for this sample
                if (metrics.LIBRARY == null) {
                    switch (metrics.CATEGORY) {
                        case FIRST_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 4);
                            Assert.assertEquals(metrics.PF_READS, 4);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                            Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                            Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                            Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 19);
                            break;
                        case SECOND_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 4);
                            Assert.assertEquals(metrics.PF_READS, 4);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                            Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                            Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                            Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 3);
                            break;
                        case PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 8);
                            Assert.assertEquals(metrics.PF_READS, 8);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                            Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                            Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                            Assert.assertEquals(metrics.PF_ALIGNED_BASES, 0);
                            Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 22);
                            break;
                        case UNPAIRED:
                        default:
                            Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                    }
                }
                else if (metrics.LIBRARY.equals("lib1")) {
                    // Only one read group in this library so library and RG metrics should be identical
                    switch (metrics.CATEGORY) {
                        case FIRST_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 2);
                            Assert.assertEquals(metrics.PF_READS, 2);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                            Assert.assertEquals(metrics.BAD_CYCLES, 19);
                            break;
                        case SECOND_OF_PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 2);
                            Assert.assertEquals(metrics.PF_READS, 2);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.BAD_CYCLES, 3);
                            break;
                        case PAIR:
                            Assert.assertEquals(metrics.TOTAL_READS, 4);
                            Assert.assertEquals(metrics.PF_READS, 4);
                            Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                            Assert.assertEquals(metrics.BAD_CYCLES, 22);
                            break;
                        case UNPAIRED:
                        default:
                            Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                    }

                }
                else if (metrics.LIBRARY.equals("lib2")) {
                    if (metrics.READ_GROUP == null) {
                        switch (metrics.CATEGORY) {
                            case FIRST_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                                break;
                            case SECOND_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                                break;
                            case PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 4);
                                Assert.assertEquals(metrics.PF_READS, 4);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                                break;
                            case UNPAIRED:
                            default:
                                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                        }
                    }
                    else if (metrics.READ_GROUP.equals("i")) {
                        switch (metrics.CATEGORY) {
                            case FIRST_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                                break;
                            case SECOND_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                                break;
                            case PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                                break;
                            case UNPAIRED:
                            default:
                                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                        }
                    }
                    else if (metrics.READ_GROUP.equals("i2")) {
                        switch (metrics.CATEGORY) {
                            case FIRST_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 27);
                                break;
                            case SECOND_OF_PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 1);
                                Assert.assertEquals(metrics.PF_READS, 1);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                                break;
                            case PAIR:
                                Assert.assertEquals(metrics.TOTAL_READS, 2);
                                Assert.assertEquals(metrics.PF_READS, 2);
                                Assert.assertEquals(metrics.PF_NOISE_READS, 0);
                                Assert.assertEquals(metrics.BAD_CYCLES, 30);
                                break;
                            case UNPAIRED:
                            default:
                                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
                        }
                    }
                    else {
                        Assert.fail("Data does not contain this read group: " + metrics.READ_GROUP);
                    }

                }
                else {
                    Assert.fail("Data does not contain this library: " + metrics.LIBRARY);
                }
            }
            else {
                Assert.fail("Data does not contain this sample: " + metrics.SAMPLE);
            }
        }
    }

    @Test
    public void testChimeras() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test_chimeras.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "MAX_INSERT_SIZE=20",
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.CATEGORY == AlignmentSummaryMetrics.Category.FIRST_OF_PAIR) {
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_CHIMERAS, 5D / 6, 0.0001);
            }
            if (metrics.CATEGORY == AlignmentSummaryMetrics.Category.SECOND_OF_PAIR) {
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_CHIMERAS, 3D / 6, 0.0001);
            }
        }
    }
}
