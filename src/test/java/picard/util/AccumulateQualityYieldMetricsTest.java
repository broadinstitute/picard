/*
 * The MIT License
 *
 * Copyright (c) 2021 The Broad Institute
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

package picard.util;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.analysis.CollectQualityYieldMetrics;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class AccumulateQualityYieldMetricsTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/quality");

    public String getCommandLineProgramName() {
        return AccumulateQualityYieldMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input1 = new File(TEST_DATA_DIR, "insert_size_metrics_test.quality_yield_metrics");
        final File input2 = new File(TEST_DATA_DIR, "insert_size_metrics_test_fake.quality_yield_metrics");
        final File outfile = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input1.getAbsolutePath(),
                "INPUT=" + input2.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(), 1);

        final CollectQualityYieldMetrics.QualityYieldMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 52 * 2 - 2);
        Assert.assertEquals(metrics.PF_READS, 52 * 2 - 2);
        Assert.assertEquals(metrics.READ_LENGTH, 102);
        Assert.assertEquals(metrics.TOTAL_BASES, 5252 * 2 - 2);
        Assert.assertEquals(metrics.PF_BASES, 5252 * 2 - 2);
        Assert.assertEquals(metrics.Q20_BASES, 3532 * 2 - 2);
        Assert.assertEquals(metrics.PF_Q20_BASES, 3532 * 2 - 2);
        Assert.assertEquals(metrics.Q30_BASES, 3145 * 2 - 5);
        Assert.assertEquals(metrics.PF_Q30_BASES, 3145 * 2 - 5);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 6497 * 2 - 7);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 6497 * 2 - 7);

    }
}