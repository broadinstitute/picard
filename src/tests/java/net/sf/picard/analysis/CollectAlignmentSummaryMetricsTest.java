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

package net.sf.picard.analysis;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;

import net.sf.picard.analysis.CollectAlignmentSummaryMetrics;
import org.testng.Assert;
import org.testng.annotations.Test;

import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.analysis.AlignmentSummaryMetrics;

/**
 * Tests CollectAlignmentSummaryStatistics
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class CollectAlignmentSummaryMetricsTest {
    private static File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam");
    
    @Test
    public void test() throws IOException {
        CollectAlignmentSummaryMetrics program = new CollectAlignmentSummaryMetrics();
        program.INPUT = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        program.OUTPUT = File.createTempFile("alignmentMetrics", ".txt");
        program.OUTPUT.deleteOnExit();
        program.REFERENCE_SEQUENCE = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        Assert.assertEquals(program.doWork(), 0);
        
        MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<AlignmentSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(program.OUTPUT));
        
        for (AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 7);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 3);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 59);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 19.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                break;
            case SECOND_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 9);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 7);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 239);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 18);
                Assert.assertEquals(metrics.PF_READS, 16);
                Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 10);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 298);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
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
        CollectAlignmentSummaryMetrics program = new CollectAlignmentSummaryMetrics();
        program.INPUT = new File(TEST_DATA_DIR, "summary_alignment_bisulfite_test.sam");
        program.OUTPUT = File.createTempFile("alignmentMetrics", ".txt");
        program.OUTPUT.deleteOnExit();
        program.REFERENCE_SEQUENCE = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        program.IS_BISULFITE_SEQUENCED = true;
        NumberFormat format =  NumberFormat.getInstance();
        format.setMaximumFractionDigits(4);

        
        Assert.assertEquals(program.doWork(), 0);

        MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<AlignmentSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(program.OUTPUT));

        for (AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                // 19 no-calls, one potentially methylated base, one mismatch at a potentially methylated base
                Assert.assertEquals(metrics.TOTAL_READS, 1);
                Assert.assertEquals(metrics.PF_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 20.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 20);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(20/(double)100));
                break;
            case SECOND_OF_PAIR:
                // Three no-calls, two potentially methylated bases
                Assert.assertEquals(metrics.TOTAL_READS, 1);
                Assert.assertEquals(metrics.PF_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 101);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(3/(double)99));
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 2);
                Assert.assertEquals(metrics.PF_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_BASES, 202);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 23);
                Assert.assertEquals(format.format(metrics.PF_HQ_ERROR_RATE), format.format(23/(double)199));
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }


    @Test
    public void testNoReference() throws IOException {
        CollectAlignmentSummaryMetrics program = new CollectAlignmentSummaryMetrics();
        program.INPUT = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        program.OUTPUT = File.createTempFile("alignmentMetrics", ".txt");
        program.OUTPUT.deleteOnExit();
        program.REFERENCE_SEQUENCE = null;
        Assert.assertEquals(program.doWork(), 0);

        MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<AlignmentSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(program.OUTPUT));

        for (AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
            case FIRST_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 7);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 19);
                break;
            case SECOND_OF_PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 9);
                Assert.assertEquals(metrics.PF_READS, 9);
                Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 3);
                break;
            case PAIR:
                Assert.assertEquals(metrics.TOTAL_READS, 18);
                Assert.assertEquals(metrics.PF_READS, 16);
                Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 0);
                Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 0);
                Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 0.0);
                Assert.assertEquals(metrics.BAD_CYCLES, 22);
                break;
            case UNPAIRED:
            default:
                Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

}
