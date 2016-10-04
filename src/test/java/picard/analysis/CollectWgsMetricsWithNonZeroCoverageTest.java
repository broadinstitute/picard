package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.analysis.CollectWgsMetricsWithNonZeroCoverage.*;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Tests for CollectWgsMetricsWithNonZeroCoverage.
 */
public class CollectWgsMetricsWithNonZeroCoverageTest  extends CommandLineProgramTest {

    private final static File TEST_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectWgsMetricsWithNonZeroCoverage.class.getSimpleName();
    }

    @Test
    public void testWithoutIntervals() throws IOException {
        final File input = new File(TEST_DIR, "forMetrics.sam");
        final File outfile = File.createTempFile("test", ".wgs_metrics");
        final File pdffile = File.createTempFile("test", ".wgs_metrics.pdf");
        final File ref = new File(TEST_DIR, "merger.fasta");
        final int sampleSize = 1000;
        outfile.deleteOnExit();
        pdffile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "SAMPLE_SIZE=" + sampleSize,
                "CHART=" + pdffile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<WgsMetricsWithNonZeroCoverage , Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final WgsMetricsWithNonZeroCoverage metrics : output.getMetrics()) {
            if (metrics.CATEGORY == WgsMetricsWithNonZeroCoverage.Category.WHOLE_GENOME) {
                Assert.assertEquals(metrics.GENOME_TERRITORY, 1210);
                Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.271403);
                Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.182149);
                Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.091075);
                Assert.assertEquals(metrics.PCT_1X, 0.107438);
            } else {
                Assert.assertEquals(metrics.GENOME_TERRITORY, 130);
                Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.271403);
                Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.182149);
                Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.091075);
                Assert.assertEquals(metrics.PCT_1X, 1.0);
            }
        }

        for (final Histogram<Integer> histogram : output.getAllHistograms()) {
            if (histogram.getValueLabel().equals("count_WHOLE_GENOME")) {
                Assert.assertEquals(histogram.get(0).getValue(), 1080d);
            } else {
                Assert.assertEquals(histogram.get(0).getValue(), 0d);
            }
            Assert.assertEquals(histogram.get(1).getValue(), 9d);
            Assert.assertEquals(histogram.get(2).getValue(), 35d);
            Assert.assertEquals(histogram.get(3).getValue(), 86d);
            Assert.assertEquals(histogram.get(4).getValue(), 0d);
        }
    }

    @Test
    public void testWithIntervals() throws IOException {
        final File input = new File(TEST_DIR, "forMetrics.sam");
        final File outfile = File.createTempFile("test", ".wgs_metrics");
        final File pdffile = File.createTempFile("test", ".wgs_metrics.pdf");
        final File ref = new File(TEST_DIR, "merger.fasta");
        final File intervals = new File(TEST_DIR, "largeIntervals.interval_list");
        final int sampleSize = 1000;
        outfile.deleteOnExit();
        pdffile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "INTERVALS=" + intervals.getAbsolutePath(),
                "SAMPLE_SIZE=" + sampleSize,
                "CHART=" + pdffile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<WgsMetricsWithNonZeroCoverage , Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final WgsMetricsWithNonZeroCoverage metrics : output.getMetrics()) {
            if (metrics.CATEGORY == WgsMetricsWithNonZeroCoverage.Category.WHOLE_GENOME) {
                Assert.assertEquals(metrics.GENOME_TERRITORY, 404);
                Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.271403);
                Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.182149);
                Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.091075);
                Assert.assertEquals(metrics.PCT_1X, 0.321782);
            } else {
                Assert.assertEquals(metrics.GENOME_TERRITORY, 130);
                Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.271403);
                Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.182149);
                Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.091075);
                Assert.assertEquals(metrics.PCT_1X, 1.0);
            }
        }

        for (final Histogram<Integer> histogram : output.getAllHistograms()) {
            if (histogram.getValueLabel().equals("count_WHOLE_GENOME")) {
                Assert.assertEquals(histogram.get(0).getValue(), 274d);
            } else {
                Assert.assertEquals(histogram.get(0).getValue(), 0d);
            }
            Assert.assertEquals(histogram.get(1).getValue(), 9d);
            Assert.assertEquals(histogram.get(2).getValue(), 35d);
            Assert.assertEquals(histogram.get(3).getValue(), 86d);
            Assert.assertEquals(histogram.get(4).getValue(), 0d);

        }
    }
}
