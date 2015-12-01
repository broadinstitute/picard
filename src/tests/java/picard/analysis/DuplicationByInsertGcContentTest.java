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

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
/**
 * Created by hogstrom on 10/18/15.
 */
public class DuplicationByInsertGcContentTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return DuplicationByInsertGcContent.class.getSimpleName();
    }

    @Test
    public void testGCcontent() throws IOException {
        //final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        //final File ref = new File(TEST_DATA_DIR, "merger.fasta");
        //final File ref = new File("testdata/picard/quality/chrM.reference.fasta");
        //final File input = new File("testdata/picard/quality/chrMReadsDiffereingLengths.sam");
        //final File ref = new File("/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta");
        //final File input = new File(TEST_DATA_DIR, "test_chrm21_rg2.bam");
        final File ref = new File("testdata/picard/metrics/chrMNO.reference.fasta");
        final File input = new File("/Users/hogstrom/Documents/code/picard/testdata/picard/sam/CollectGcBiasMetrics/CollectGcBias6098159690966723109.bam");
        final File outfile = File.createTempFile("test", ".insert_GC_by_dup");
        final File pdf = File.createTempFile("test", ".pdf");
        //outfile.deleteOnExit();
        System.out.println("outfile = " +  outfile.getAbsolutePath());
        pdf.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "CHART=" + pdf.getAbsolutePath(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        /*final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getAllHistograms().size(), 5);

        for (final InsertSizeMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.PAIR_ORIENTATION.name(), "FR");
            if (metrics.LIBRARY == null) {  // SAMPLE or ALL_READS level
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 41);
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
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            }
        }*/
    }
}

