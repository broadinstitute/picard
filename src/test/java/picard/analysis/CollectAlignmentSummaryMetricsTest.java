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
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.util.TestNGUtil;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests CollectAlignmentSummaryStatistics
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class CollectAlignmentSummaryMetricsTest extends CommandLineProgramTest {
    public static final File TEST_DATA_DIR = new File("testdata/picard/sam/AlignmentSummaryMetrics");

    public String getCommandLineProgramName() {
        return CollectAlignmentSummaryMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile = getTempOutputFile("test", ".txt");
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

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
        final File outfile = getTempOutputFile("testBisulfite", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "IS_BISULFITE_SEQUENCED=true"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final NumberFormat format = NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

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
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(21 / (double) 99));
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
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(4 / (double) 99));
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 2);
                    Assert.assertEquals(metrics.PF_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 202);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 12.5D);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 202);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.126263);// 25D/198D
                    Assert.assertEquals(metrics.BAD_CYCLES, 25);
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(25 / (double) 198));
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
        final File outfile = getTempOutputFile("testBisulfiteButNot", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "IS_BISULFITE_SEQUENCED=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final NumberFormat format = NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

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
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(23 / (double) 101));
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
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(6 / (double) 101));
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 2);
                    Assert.assertEquals(metrics.PF_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 202);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 14.5D);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 202);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, 0.143564);// 29D/202D
                    Assert.assertEquals(metrics.BAD_CYCLES, 29);
                    Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(29 / (double) 202));
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
        final File outfile = getTempOutputFile("testNoReference", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "COLLECT_ALIGNMENT_INFORMATION=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

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
        final File outfile = getTempOutputFile("testZeroLengthReads", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "COLLECT_ALIGNMENT_INFORMATION=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }
        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            // test that it doesn't blow up
        }
    }

    @Test
    public void testMultipleLevelsOfMetrics() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test_multiple.sam");
        final File outfile = getTempOutputFile("testMultipleLevelsOfMetrics", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=ALL_READS",
                "METRIC_ACCUMULATION_LEVEL=SAMPLE",
                "METRIC_ACCUMULATION_LEVEL=LIBRARY",
                "METRIC_ACCUMULATION_LEVEL=READ_GROUP",
                "COLLECT_ALIGNMENT_INFORMATION=false"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

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
            } else if (metrics.SAMPLE.equals("Ma")) {
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
            } else if (metrics.SAMPLE.equals("Pa")) {
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
                } else if (metrics.LIBRARY.equals("lib1")) {
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

                } else if (metrics.LIBRARY.equals("lib2")) {
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
                    } else if (metrics.READ_GROUP.equals("i")) {
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
                    } else if (metrics.READ_GROUP.equals("i2")) {
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
                    } else {
                        Assert.fail("Data does not contain this read group: " + metrics.READ_GROUP);
                    }

                } else {
                    Assert.fail("Data does not contain this library: " + metrics.LIBRARY);
                }
            } else {
                Assert.fail("Data does not contain this sample: " + metrics.SAMPLE);
            }
        }
    }

    @Test
    public void testChimeras() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test_chimeras.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile = getTempOutputFile("testChimeras", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "MAX_INSERT_SIZE=20",
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.CATEGORY == AlignmentSummaryMetrics.Category.FIRST_OF_PAIR) {
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_CHIMERAS, 5D / 6, 0.0001);
            }
            if (metrics.CATEGORY == AlignmentSummaryMetrics.Category.SECOND_OF_PAIR) {
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_CHIMERAS, 3D / 6, 0.0001);
            }
        }
    }

    @Test
    public void testAdapterReads() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test_adapter_reads.sam");
        final File outfile = getTempOutputFile("testAdapterReads", ".txt");

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "MAX_INSERT_SIZE=200",
                "REFERENCE_SEQUENCE=" + CHR_M_REFERENCE.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outfile)) {
            output.read(reader);
        }

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.CATEGORY == AlignmentSummaryMetrics.Category.FIRST_OF_PAIR) {
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_ADAPTER, 0D, 0.0001);
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_PF_READS_ALIGNED, 1D, 0.0001);
            }
            if (metrics.CATEGORY == AlignmentSummaryMetrics.Category.SECOND_OF_PAIR) {
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_ADAPTER, 1D, 0.0001);
                TestNGUtil.compareDoubleWithAccuracy(metrics.PCT_PF_READS_ALIGNED, 0D, 0.0001);
            }
        }
    }


    @DataProvider
    Object[][] fileForTestReadLengthHistogram(){
        return new Object[][]{
                new Object[]{"summary_alignment_stats_test.sam"},
                new Object[]{"summary_alignment_stats_test2.sam"},
                new Object[]{"summary_alignment_stats_test3.sam"}
        };
    }

    @Test(dataProvider = "fileForTestReadLengthHistogram")
    public void testReadLengthHistogram(final String fileToUse) throws IOException {
        final File input = new File(TEST_DATA_DIR, fileToUse);
        final File outFile = getTempOutputFile("testReadLengthHistogram", ".txt");

        final List<String> argsList = new ArrayList<>();
        final File outHist = getTempOutputFile("testReadLengthHistogram", ".pdf");

        argsList.add("INPUT=" + input.getAbsolutePath());
        argsList.add("OUTPUT=" + outFile.getAbsolutePath());
        argsList.add("HISTOGRAM_FILE=" + outHist);

        Assert.assertEquals(runPicardCommandLine(argsList.toArray(new String[0])),0);

        Assert.assertTrue(outHist.exists());
    }


    @DataProvider()
    Object[][] TrueFalse() {
        return new Object[][]{
                {true},
                {false},
        };
    }

    @Test(dataProvider = "TrueFalse")
    public void testReadLengthHistogram(final boolean plotChart) throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test3.sam");
        final File outFile = getTempOutputFile("testReadLengthHistogram", ".txt");

        final List<String> argsList = new ArrayList<>();
        argsList.add("INPUT=" + input.getAbsolutePath());
        argsList.add("OUTPUT=" + outFile.getAbsolutePath());
        final File outHist;

        if (plotChart) {
            outHist = getTempOutputFile("testReadLengthHistogram", ".pdf");
            argsList.add("HISTOGRAM_FILE=" + outHist);
        } else {
            outHist = null;
        }

        final String[] args = argsList.toArray(new String[0]);

        Assert.assertEquals(runPicardCommandLine(args), 0);

        if (plotChart) {
            Assert.assertTrue(outHist.exists());
        }

        final MetricsFile<AlignmentSummaryMetrics, Integer> output = new MetricsFile<>();
        try (FileReader reader = new FileReader(outFile)) {
            output.read(reader);
        }

        Assert.assertFalse(output.getMetrics().isEmpty());


        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.CATEGORY != AlignmentSummaryMetrics.Category.SECOND_OF_PAIR) {
                Assert.assertEquals(metrics.PCT_HARDCLIP, 10 / (double) 66, 0.0001);
                Assert.assertEquals(metrics.PCT_SOFTCLIP, 23 / (double) 66, 0.0001);
                Assert.assertEquals(metrics.AVG_POS_3PRIME_SOFTCLIP_LENGTH, 22 / (double) 2, 0.0001);
            } else {
                Assert.assertEquals(metrics.PCT_HARDCLIP, 0D);
                Assert.assertEquals(metrics.PCT_SOFTCLIP, 0D);
                Assert.assertEquals(metrics.AVG_POS_3PRIME_SOFTCLIP_LENGTH, 0D);
            }
        }

        Assert.assertFalse(output.getAllHistograms().isEmpty());

        for (final Histogram<Integer> histogram : output.getAllHistograms()) {
            switch (histogram.getValueLabel()) {
                case "PAIRED_TOTAL_LENGTH_COUNT":
                case "UNPAIRED_TOTAL_LENGTH_COUNT":
                    Assert.assertEquals(histogram.getSum(), 66D); //1+2+11+10+42
                    break;
                case "PAIRED_ALIGNED_LENGTH_COUNT":
                case "UNPAIRED_ALIGNED_LENGTH_COUNT":
                    Assert.assertEquals(histogram.getSum(), 43D); //1+2+10+10+10+10;
                    break;
            }

            Assert.assertFalse(histogram.isEmpty());
            if (histogram.getValueLabel().equals("PAIRED_TOTAL_LENGTH_COUNT")) {
                for (int i = 0; i < histogram.getMax(); i++) {
                    switch (i) {
                        case 1:
                        case 2:
                        case 10:
                        case 11:
                        case 42:
                            Assert.assertEquals(histogram.get(i).getValue(), 1D);
                            break;
                        default:
                            Assert.assertTrue(
                                    histogram.get(i) == null ||
                                            histogram.get(i).getValue() == 0D);
                            break;
                    }
                }
            }
        }
    }
}
