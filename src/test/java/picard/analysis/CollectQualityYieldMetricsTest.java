/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.function.Consumer;

/**
 * Created by kbergin on 11/23/15.
 */
public class CollectQualityYieldMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectQualityYieldMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetrics metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 52);
        Assert.assertEquals(metrics.PF_READS, 52);
        Assert.assertEquals(metrics.READ_LENGTH, 101);
        Assert.assertEquals(metrics.TOTAL_BASES, 5252);
        Assert.assertEquals(metrics.PF_BASES, 5252);
        Assert.assertEquals(metrics.Q20_BASES, 3532);
        Assert.assertEquals(metrics.PF_Q20_BASES, 3532);
        Assert.assertEquals(metrics.Q30_BASES, 3145);
        Assert.assertEquals(metrics.PF_Q30_BASES, 3145);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 6497);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 6497);

    }


    @Test
    public void testFlowMode() throws IOException {
        final File input = new File(TEST_DATA_DIR, "subsample.bam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "USE_ORIGINAL_QUALITIES=false",
                "FLOW_MODE=true"
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetricsFlow, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetricsFlow metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 56);
        Assert.assertEquals(metrics.PF_READS, 56);
        Assert.assertEquals(metrics.READ_LENGTH, 285);
        Assert.assertEquals(metrics.TOTAL_BASES, 15983);
        Assert.assertEquals(metrics.PF_BASES, 15983);
        Assert.assertEquals(metrics.Q20_BASES, 15494);
        Assert.assertEquals(metrics.PF_Q20_BASES, 15494);
        Assert.assertEquals(metrics.Q30_BASES, 14786);
        Assert.assertEquals(metrics.PF_Q30_BASES, 14786);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 30589);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 30589);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30, 102);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25, 196);
    }

    @Test
    public void testFlowModeReverseReadsSameResults() throws IOException {
        final File input = new File(TEST_DATA_DIR, "subsample.bam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "USE_ORIGINAL_QUALITIES=false",
                "FLOW_MODE=true"
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final File inputReverse = new File(TEST_DATA_DIR, "subsample.reverse.bam");
        final File outfileReverse   = File.createTempFile("test", "reverse.quality_yield_metrics");
        outfileReverse.deleteOnExit();
        final String[] argsReverse = new String[] {
                "INPUT="  + inputReverse.getAbsolutePath(),
                "OUTPUT=" + outfileReverse.getAbsolutePath(),
                "USE_ORIGINAL_QUALITIES=false",
                "FLOW_MODE=true"
        };
        Assert.assertEquals(runPicardCommandLine(argsReverse), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetricsFlow, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetricsFlow, ?> outputReverse = new MetricsFile<>();
        outputReverse.read(new FileReader(outfileReverse));

        Assert.assertEquals(output.getMetrics().size(),1);
        Assert.assertEquals(outputReverse.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetricsFlow metrics = output.getMetrics().get(0);
        final CollectQualityYieldMetrics.QualityYieldMetricsFlow metricsReverse = outputReverse.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, metricsReverse.TOTAL_READS);
        Assert.assertEquals(metrics.PF_READS, metricsReverse.PF_READS);
        Assert.assertEquals(metrics.READ_LENGTH, metricsReverse.READ_LENGTH);
        Assert.assertEquals(metrics.TOTAL_BASES, metricsReverse.TOTAL_BASES);
        Assert.assertEquals(metrics.PF_BASES, metricsReverse.PF_BASES);
        Assert.assertEquals(metrics.Q20_BASES, metricsReverse.Q20_BASES);
        Assert.assertEquals(metrics.PF_Q20_BASES, metricsReverse.PF_Q20_BASES);
        Assert.assertEquals(metrics.Q30_BASES, metricsReverse.Q30_BASES);
        Assert.assertEquals(metrics.PF_Q30_BASES, metricsReverse.PF_Q30_BASES);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, metricsReverse.Q20_EQUIVALENT_YIELD);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, metricsReverse.PF_Q20_EQUIVALENT_YIELD);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30,
                metricsReverse.READ_LENGTH_AVG_Q_ABOVE_30);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25,
                metricsReverse.READ_LENGTH_AVG_Q_ABOVE_25);
    }

    @Test
    public void testMultiFlow() throws IOException {
        final File input = new File(TEST_DATA_DIR, "subsample.bam");
        final File outfile   = File.createTempFile("test", ".quality_yield_metrics");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "PROGRAM=null",
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectQualityYieldMetrics.name(),
                "EXTRA_ARGUMENT=CollectQualityYieldMetrics::FLOW_MODE=true",
                "EXTRA_ARGUMENT=CollectQualityYieldMetrics::USE_ORIGINAL_QUALITIES=false"
        };

        Assert.assertEquals(runPicardCommandLine(CollectMultipleMetrics.class.getSimpleName(), args), 0);

        final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetricsFlow, ?> output = new MetricsFile<>();
        output.read(new FileReader(outfile + ".quality_yield_metrics"));

        Assert.assertEquals(output.getMetrics().size(),1);

        final CollectQualityYieldMetrics.QualityYieldMetricsFlow metrics = output.getMetrics().get(0);
        Assert.assertEquals(metrics.TOTAL_READS, 56);
        Assert.assertEquals(metrics.PF_READS, 56);
        Assert.assertEquals(metrics.READ_LENGTH, 285);
        Assert.assertEquals(metrics.TOTAL_BASES, 15983);
        Assert.assertEquals(metrics.PF_BASES, 15983);
        Assert.assertEquals(metrics.Q20_BASES, 15494);
        Assert.assertEquals(metrics.PF_Q20_BASES, 15494);
        Assert.assertEquals(metrics.Q30_BASES, 14786);
        Assert.assertEquals(metrics.PF_Q30_BASES, 14786);
        Assert.assertEquals(metrics.Q20_EQUIVALENT_YIELD, 30589);
        Assert.assertEquals(metrics.PF_Q20_EQUIVALENT_YIELD, 30589);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_30, 102);
        Assert.assertEquals(metrics.READ_LENGTH_AVG_Q_ABOVE_25, 196);

    }

    @Test
    public void testMergeFlow() {
        CollectQualityYieldMetrics.QualityYieldMetricsFlow      m1 = createTestQualityYieldMetricsFlow();
        CollectQualityYieldMetrics.QualityYieldMetricsFlow      m2 = createTestQualityYieldMetricsFlow();

        m1.merge(m2);
        Assert.assertEquals(m1.TOTAL_READS, m2.TOTAL_READS * 2);
        Assert.assertEquals(m1.PF_READS, m2.PF_READS * 2);
        Assert.assertEquals(m1.READ_LENGTH, m2.READ_LENGTH);
        Assert.assertEquals(m1.PF_BASES, m2.PF_BASES * 2);
        Assert.assertEquals(m1.Q20_BASES, m2.Q20_BASES * 2);
        Assert.assertEquals(m1.PF_Q20_BASES, m2.PF_Q20_BASES * 2);
        Assert.assertEquals(m1.Q30_BASES, m2.Q30_BASES * 2);
        Assert.assertEquals(m1.PF_Q30_BASES, m2.PF_Q30_BASES * 2);
        Assert.assertEquals(m1.Q20_EQUIVALENT_YIELD, m2.Q20_EQUIVALENT_YIELD * 2);
        Assert.assertEquals(m1.READ_LENGTH_AVG_Q_ABOVE_25, m2.READ_LENGTH_AVG_Q_ABOVE_25);
        Assert.assertEquals(m1.READ_LENGTH_AVG_Q_ABOVE_25, 100 - m1.histogramGenerator.skipBases);

        Assert.assertEquals(m1.READ_LENGTH_AVG_Q_ABOVE_30, m2.READ_LENGTH_AVG_Q_ABOVE_30);
        Assert.assertEquals(m1.READ_LENGTH_AVG_Q_ABOVE_30, 100 - m1.histogramGenerator.skipBases);

    }

    private CollectQualityYieldMetrics.QualityYieldMetricsFlow createTestQualityYieldMetricsFlow() {

        HistogramGenerator hg  = createTestHistograms();
        CollectQualityYieldMetrics.QualityYieldMetricsFlow m = new CollectQualityYieldMetrics.QualityYieldMetricsFlow(false, hg);

        m.TOTAL_READS = 52;
        m.PF_READS = 52;
        m.READ_LENGTH = 101;
        m.TOTAL_BASES = 5252;
        m.PF_BASES = 5252;
        m.Q20_BASES = 3532;
        m.PF_Q20_BASES = 3532;
        m.Q30_BASES = 3145;
        m.PF_Q30_BASES = 3145;
        m.Q20_EQUIVALENT_YIELD = 6497;
        m.PF_Q20_EQUIVALENT_YIELD = 6497;
        m.calculateDerivedFields();
        return m;
    }

    private HistogramGenerator createTestHistograms(){
        long [] firstReadCounts = new long[100];
        Arrays.fill(firstReadCounts,100);
        long [] secondReadCounts = firstReadCounts.clone();
        double[] firstReadTotals = new double[100];
        Arrays.fill(firstReadTotals, 100*50);
        double[] secondReadsTotals = firstReadTotals.clone();

        double[] firstReadTotalProbs = new double[100];
        Arrays.fill(firstReadTotalProbs, 100*0.0001);
        double[] secondReadsTotalProbs = firstReadTotalProbs.clone();

        return new HistogramGenerator(firstReadTotals, firstReadTotalProbs, firstReadCounts,
                secondReadsTotals, secondReadsTotalProbs, secondReadCounts, 50);
    }


}
