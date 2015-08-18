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
 * Tests CollectWgsMetrics
 */
public class CollectWgsMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectWgsMetrics.class.getSimpleName();
    }

    /**
     * Important note:
     *   This test is simply confirming that the correct number of sites are being evaluated by the tool.
     *   It is not testing correctness of the results.
     *   The code in general requires significant further testing of correctness.
     */
    @Test
    public void testIntervals() throws IOException {
        final File input = new File(TEST_DATA_DIR, "aligned.sam");
        final File outfile = File.createTempFile("test", ".wgs_metrics");
        final File ref = new File(TEST_DATA_DIR, "merger.fasta");
        final File intervals = new File(TEST_DATA_DIR, "several.interval_list");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "INTERVALS=" + intervals.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetrics.WgsMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.GENOME_TERRITORY, 5);
        }
    }
}
