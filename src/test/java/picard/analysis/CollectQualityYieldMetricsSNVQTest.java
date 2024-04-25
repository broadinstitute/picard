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
 * Created by dror27 on 24/01/2024.
 */
public class CollectQualityYieldMetricsSNVQTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");
    private static final double DOUBLE_EQUALS_EPSILON = 0.001;

    public String getCommandLineProgramName() {
        return CollectQualityYieldMetricsSNVQ.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "snvq_metrics_test.bam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetricsSNVQ.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetricsSNVQ.QualityYieldMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 26577);
        Assert.assertEquals(metrics.PF_READS, 26577);
        Assert.assertEquals(metrics.READ_LENGTH, 173);
        Assert.assertEquals(metrics.TOTAL_BASES, 4605838);
        Assert.assertEquals(metrics.PF_BASES, 4605838);
        Assert.assertEquals(metrics.Q20_BASES, 4450324);
        Assert.assertEquals(metrics.PF_Q20_BASES, 4450324);
        Assert.assertEquals(metrics.PF_Q30_BASES, 4069304);
        Assert.assertEquals(metrics.Q40_BASES, 0);
        Assert.assertEquals(metrics.PF_Q40_BASES, 0);
        Assert.assertEquals(metrics.PCT_PF_Q20_BASES, 0.966235, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_PF_Q30_BASES, 0.88351, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_PF_Q40_BASES, 0.0, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.TOTAL_SNVQ, 13817514);
        Assert.assertEquals(metrics.PF_SNVQ, 13817514);
        Assert.assertEquals(metrics.Q20_SNVQ, 13816364);
        Assert.assertEquals(metrics.PF_Q20_SNVQ, 13816364);
        Assert.assertEquals(metrics.Q30_SNVQ, 13816364);
        Assert.assertEquals(metrics.PF_Q30_SNVQ, 13816364);
        Assert.assertEquals(metrics.Q40_SNVQ, 13816364);
        Assert.assertEquals(metrics.PF_Q40_SNVQ, 13816364);
        Assert.assertEquals(metrics.PCT_Q20_SNVQ, 0.999917, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_Q30_SNVQ, 0.999917, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_Q40_SNVQ, 0.999917, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_PF_Q20_SNVQ, 0.999917, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_PF_Q30_SNVQ, 0.999917, DOUBLE_EQUALS_EPSILON);
        Assert.assertEquals(metrics.PCT_PF_Q40_SNVQ, 0.999917, DOUBLE_EQUALS_EPSILON);
    }
}
