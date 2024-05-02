/*
 * The MIT License
 *
 * Copyright (c) 2024 The Broad Institute
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
 * Created by dror27 on 21/12/23.
 */
public class CollectQualityYieldMetricsFlowTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectQualityYieldMetricsFlow.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test_flow.sam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics_flow_space");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetricsFlow.QualityYieldMetricsFlow, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetricsFlow.QualityYieldMetricsFlow metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 56);
        Assert.assertEquals(metrics.PF_READS, 56);
        Assert.assertEquals(metrics.MEAN_PF_READ_NUMBER_OF_FLOWS, 375);
        Assert.assertEquals(metrics.PF_FLOWS, 21053);
        Assert.assertEquals(metrics.PF_Q20_FLOWS, 20667);
        Assert.assertTrue(metrics.PCT_PF_Q20_FLOWS > 0);
        Assert.assertEquals(metrics.PF_Q30_FLOWS, 20256);
        Assert.assertTrue(metrics.PCT_PF_Q30_FLOWS > 0);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 41177);
    }
}
