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
 * Tests CollectWgsMetricsFromSampledSites
 */
public class CollectWgsMetricsFromSampledSitesTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectWgsMetricsFromSampledSites.class.getSimpleName();
    }

    @Test
    public void testOnePos() throws IOException {
        final File input = new File(TEST_DATA_DIR, "forMetrics.sam");
        final File outfile = File.createTempFile("test", ".wgs_metrics");
        final File ref = new File(TEST_DATA_DIR, "merger.fasta");
        final File intervals = new File(TEST_DATA_DIR, "onePos.interval_list");
        final int sampleSize = 1000;
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "INTERVALS=" + intervals.getAbsolutePath(),
                "SAMPLE_SIZE=" + sampleSize
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetricsFromSampledSites.SampledWgsMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetricsFromSampledSites.SampledWgsMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetrics.WgsMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.GENOME_TERRITORY, 1);
            Assert.assertEquals(metrics.MEAN_COVERAGE, 3.0);
            Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.272727); // 3 of 11
            Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.181818); // 2 of 11
            Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.090909); // 1 of 9
            Assert.assertEquals(metrics.PCT_EXC_BASEQ, 0.090909); // 1 of 9
            Assert.assertEquals(metrics.HET_SNP_SENSITIVITY, 0.34655, .02);
        }
    }

    @Test
    public void testContiguousIntervals() throws IOException {
        final File input = new File(TEST_DATA_DIR, "forMetrics.sam");
        final File outfile = File.createTempFile("test", ".wgs_metrics");
        final File ref = new File(TEST_DATA_DIR, "merger.fasta");
        final File intervals = new File(TEST_DATA_DIR, "contiguous.interval_list");
        final int sampleSize = 1000;
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "INTERVALS=" + intervals.getAbsolutePath(),
                "SAMPLE_SIZE=" + sampleSize
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetrics.WgsMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.GENOME_TERRITORY, 5);
            Assert.assertEquals(metrics.MEAN_COVERAGE, 2.6);
            Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.0);
            Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.066667);
            Assert.assertEquals(metrics.HET_SNP_SENSITIVITY, 0.393802, .02);
        }
    }
}
