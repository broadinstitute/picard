package picard.analysis;

import picard.analysis.CollectWgsMetricsWithNonZeroCoverage.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CollectWgsMetricsWithNonZeroCoverageTest extends CommandLineProgramTest {
    private final static File TEST_DIR = new File("testdata/picard/sam/");
    private final static String SAMPLE = "TestSample1";
    private final static String READ_GROUP_ID = "TestReadGroup1";
    private final static String PLATFORM = "ILLUMINA";
    private final static String LIBRARY = "TestLibrary1";

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

    @Test
    public void testPoorQualityBases() throws IOException {
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final File testSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        testSamFile.deleteOnExit();

        /**
         *  Our test SAM looks as follows:
         *
         *   ----------   <- reads with great base qualities (60) ->  ----------
         *   ----------                                               ----------
         *   ----------                                               ----------
         *   **********   <- reads with poor base qualities (10) ->   **********
         *   **********                                               **********
         *   **********                                               **********
         *
         *  We exclude half of the bases because they are low quality.
         *  We do not exceed the coverage cap (3), thus none of the bases is excluded as such.
         *
         */

        final SAMRecordSetBuilder setBuilder = CollectWgsMetricsTestUtils.createTestSAMBuilder(reference, READ_GROUP_ID, SAMPLE, PLATFORM, LIBRARY);
        setBuilder.setReadLength(10);
        for (int i = 0; i < 3; i++){
            setBuilder.addPair("GreatBQRead:" + i, 0, 1, 30, false, false, "10M", "10M", false, true, 60);
        }

        for (int i = 0; i < 3; i++){
            setBuilder.addPair("PoorBQRead:" + i, 0, 1, 30, false, false, "10M", "10M", false, true, 10);
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, testSamFile);
        setBuilder.forEach(writer::addAlignment);
        writer.close();

        final File outfile = File.createTempFile("testExcludedBases-metrics", ".txt");
        outfile.deleteOnExit();

        final File chartOutFile = File.createTempFile("testExcludedBases",".pdf");
        chartOutFile.deleteOnExit();


        final String[] args = new String[] {
                "INPUT="  + testSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=true",
                "COVERAGE_CAP=3",
                "CHART_OUTPUT=" + chartOutFile.getAbsolutePath()
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);
        Assert.assertTrue(chartOutFile.exists());

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final CollectWgsMetrics.WgsMetrics metrics = output.getMetrics().get(0);
        final CollectWgsMetrics.WgsMetrics nonZeroMetrics = output.getMetrics().get(1);


        // Some metrics should not change between with and without zero
        Assert.assertEquals(nonZeroMetrics.PCT_EXC_BASEQ, metrics.PCT_EXC_BASEQ);
        Assert.assertEquals(nonZeroMetrics.PCT_EXC_CAPPED, metrics.PCT_EXC_CAPPED);

        // Other metrics change when we ignore the zero depth bin
        Assert.assertEquals(nonZeroMetrics.GENOME_TERRITORY, 20);
        Assert.assertEquals(nonZeroMetrics.MEAN_COVERAGE, 3.0);
    }

    @Test
    public void testNoCoverage() throws IOException {
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final File testSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        testSamFile.deleteOnExit();

        final SAMRecordSetBuilder setBuilder = CollectWgsMetricsTestUtils.createTestSAMBuilder(reference, READ_GROUP_ID, SAMPLE, PLATFORM, LIBRARY);
        setBuilder.setReadLength(10);
        for (int i = 0; i < 3; i++){
            setBuilder.addPair("query:" + i, 0, 1, 30, true, true, "10M", "10M", false, true, 60);
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, testSamFile);
        setBuilder.forEach(writer::addAlignment);
        writer.close();

        final File outfile = File.createTempFile("testPoorQualityBases-metrics", ".txt");
        outfile.deleteOnExit();

        final File chartOutFile = File.createTempFile("testPoorQualityBases",".pdf");
        chartOutFile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT="  + testSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=true",
                "COVERAGE_CAP=3",
                "CHART_OUTPUT=" + chartOutFile.getAbsolutePath()
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);
        Assert.assertTrue(chartOutFile.exists());

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final CollectWgsMetrics.WgsMetrics metrics = output.getMetrics().get(0);
        final CollectWgsMetrics.WgsMetrics nonZeroMetrics = output.getMetrics().get(1);

        // Some metrics should not change between with and without zero
        Assert.assertEquals(nonZeroMetrics.PCT_EXC_BASEQ, metrics.PCT_EXC_BASEQ);
        Assert.assertEquals(nonZeroMetrics.PCT_EXC_CAPPED, metrics.PCT_EXC_CAPPED);

        // Other metrics change when we ignore the zero depth bin
        Assert.assertEquals(nonZeroMetrics.GENOME_TERRITORY, 0);
        Assert.assertEquals(nonZeroMetrics.MEAN_COVERAGE, Double.NaN);
    }
}
