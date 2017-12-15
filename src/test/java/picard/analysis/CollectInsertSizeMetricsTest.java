/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.util.RExecutor;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

/**
 * Tests multi-level CollectInsertSizeMetrics
 */
public class CollectInsertSizeMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectInsertSizeMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File outfile   = File.createTempFile("test", ".insert_size_metrics");
        final File pdf   = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "Histogram_FILE=" + pdf.getAbsolutePath(),
                "LEVEL=SAMPLE",
                "LEVEL=LIBRARY",
                "LEVEL=READ_GROUP"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final InsertSizeMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.PAIR_ORIENTATION.name(), "FR");
            if (metrics.LIBRARY==null) {  // SAMPLE or ALL_READS level
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 13);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            }
            else if (metrics.LIBRARY.equals("Solexa-41753")) { // one LIBRARY and one READ_GROUP
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 44);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.READ_PAIRS, 2);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);

            }
            else if (metrics.LIBRARY.equals("Solexa-41748") && metrics.READ_GROUP == null) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 40);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 9);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);
            }
            else if (metrics.LIBRARY.equals("Solexa-41734") && metrics.READ_GROUP == null) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 38);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 2);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 7);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.7")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 37);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 4);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 9);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.6")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 38);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 5);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 9);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.5")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 1);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.3")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 36);
                Assert.assertEquals((int)metrics.MODE_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.READ_PAIRS, 1);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_95_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);
            }
            else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }

    /**
     * Histogram Width was being incorrectly set causing histograms to be trimmed an removed inappropriately.
     * See https://github.com/broadinstitute/picard/issues/253
     * Test to be sure that the right number of histograms are being output.
     */
    @Test
    public void testHistogramWidthIsSetProperly() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File outfile = File.createTempFile("test", ".insert_size_metrics");
        final File pdf = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "Histogram_FILE=" + pdf.getAbsolutePath(),
                "LEVEL=null",
                "LEVEL=READ_GROUP"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getAllHistograms().size(), 5);
    }

    @Test
    public void testMultipleOrientationsForHistogram() throws IOException {
        final File output = new File("testdata/picard/analysis/directed/CollectInsertSizeMetrics", "multiple_orientation.sam.insert_size_metrics");
        final File pdf = File.createTempFile("test", ".pdf");
        pdf.deleteOnExit();

        final int rResult;
        rResult = RExecutor.executeFromClasspath(
                CollectInsertSizeMetrics.Histogram_R_SCRIPT,
                output.getAbsolutePath(),
                pdf.getAbsolutePath(),
                "Flags of Chad and Romania");

        Assert.assertEquals(rResult, 0);
    }

    @Test
    public void testWdithOfMetrics() throws IOException {
        final File testSamFile = File.createTempFile("CollectInsertSizeMetrics", ".bam", TEST_DATA_DIR);
        testSamFile.deleteOnExit();

        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate, true, 100);
        setBuilder.setReadLength(10);

        final int insertBy = 3; // the # of bases to increase the insert by in the records below.
        int queryIndex = 0;

        // Create records such that we have 10 records in the 10th through the 90th percentiles, and 5 for the 95th and 99th percentiles.
        // WIDTH_OF_10_PERCENT through WIDTH_OF_90_PERCENT (90 pairs total))
        for (int j = 0; j < 9; j++) {
            for (int i = 0; i < 5; i++, queryIndex++) {
                setBuilder.addPair("query:" + queryIndex, 0, 1, 50 + j*insertBy, false, false, "10M", "10M", false, true, 60);
            }
            for (int i = 0; i < 5; i++, queryIndex++) {
                setBuilder.addPair("query:" + queryIndex, 0, 1, 50 - j*insertBy, false, false, "10M", "10M", false, true, 60);
            }
        }
        // WIDTH_OF_95_PERCENT through WIDTH_OF_99_PERCENT (10 pairs total)
        for (int j = 9; j < 11; j++) {
            for (int i = 0; i < 3; i++, queryIndex++) {
                setBuilder.addPair("query:" + queryIndex, 0, 1, 50 + j*insertBy, false, false, "10M", "10M", false, true, 60);
            }
            for (int i = 0; i < 2; i++, queryIndex++) {
                setBuilder.addPair("query:" + queryIndex, 0, 1, 50 - j*insertBy, false, false, "10M", "10M", false, true, 60);
            }
        }

        // Add one to make the an odd # of pairs for the median
        setBuilder.addPair("query:" + queryIndex, 0, 1, 50, false, false, "10M", "10M", false, true, 60);

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, testSamFile);
        setBuilder.forEach(writer::addAlignment);
        writer.close();


        final File outfile = File.createTempFile("test", ".insert_size_metrics");
        final File pdf     = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + testSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "Histogram_FILE=" + pdf.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));
        
        final List<InsertSizeMetrics> metrics = output.getMetrics();
        
        Assert.assertEquals(metrics.size(), 1);
        
        final InsertSizeMetrics metric = metrics.get(0);

        Assert.assertEquals(metric.PAIR_ORIENTATION.name(), "FR");
        Assert.assertEquals((int) metric.MEDIAN_INSERT_SIZE, 59);
        Assert.assertEquals((int) metric.MODE_INSERT_SIZE, 59);
        Assert.assertEquals(metric.MIN_INSERT_SIZE, 29);
        Assert.assertEquals(metric.MAX_INSERT_SIZE, 89);
        Assert.assertEquals(metric.READ_PAIRS, 101);
        Assert.assertEquals(metric.WIDTH_OF_10_PERCENT, 1);
        Assert.assertEquals(metric.WIDTH_OF_20_PERCENT, 1 + insertBy*2);
        Assert.assertEquals(metric.WIDTH_OF_30_PERCENT, 1 + insertBy*4);
        Assert.assertEquals(metric.WIDTH_OF_40_PERCENT, 1 + insertBy*6);
        Assert.assertEquals(metric.WIDTH_OF_50_PERCENT, 1 + insertBy*8);
        Assert.assertEquals(metric.WIDTH_OF_60_PERCENT, 1 + insertBy*10);
        Assert.assertEquals(metric.WIDTH_OF_70_PERCENT, 1 + insertBy*12);
        Assert.assertEquals(metric.WIDTH_OF_80_PERCENT, 1 + insertBy*14);
        Assert.assertEquals(metric.WIDTH_OF_90_PERCENT, 1 + insertBy*16);
        Assert.assertEquals(metric.WIDTH_OF_95_PERCENT, 1 + insertBy*18);
        Assert.assertEquals(metric.WIDTH_OF_99_PERCENT, 1 + insertBy*20);
    }
}
