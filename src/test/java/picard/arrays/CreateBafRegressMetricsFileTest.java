package picard.arrays;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CreateBafRegressMetricsFileTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays");
    private static final File TEST_INPUT_FILE = new File(TEST_DATA_DIR, "BafRegress.output.txt");
    private static final File TEST_BAD_INPUT_FILE = new File(TEST_DATA_DIR, "input.vcf");

    @Test
    public void testCreateBafRegressMetricsFile() throws IOException {
        final File outputBase = File.createTempFile("testCreateBafRegressMetricsFile", "");
        outputBase.deleteOnExit();
        final File outputMetricsFile = new File(outputBase.getAbsolutePath() + "." + CreateBafRegressMetricsFile.FILE_EXTENSION);
        outputMetricsFile.deleteOnExit();

        final CreateBafRegressMetricsFile createBafRegressMetricsFile = new CreateBafRegressMetricsFile();
        createBafRegressMetricsFile.INPUT = TEST_INPUT_FILE;
        createBafRegressMetricsFile.OUTPUT = outputBase;

        Assert.assertEquals(createBafRegressMetricsFile.instanceMain(new String[0]), 0);

        final MetricsFile<BafRegressMetrics, Comparable<?>> metrics = new MetricsFile<>();
        metrics.read(new FileReader(outputMetricsFile));

        Assert.assertEquals(metrics.getMetrics().size(), 1);
        Assert.assertEquals(metrics.getMetrics().get(0).SAMPLE, "204126290052_R01C01");
        Assert.assertEquals(metrics.getMetrics().get(0).ESTIMATE, -0.000539);
        Assert.assertEquals(metrics.getMetrics().get(0).STDERR, 0.000112);
        Assert.assertEquals(metrics.getMetrics().get(0).TVAL, -4.796235);
        Assert.assertEquals(metrics.getMetrics().get(0).PVAL, 0.000002);
        Assert.assertEquals(metrics.getMetrics().get(0).LOG10_PVAL, -5.791315);
        Assert.assertEquals(metrics.getMetrics().get(0).CALL_RATE, 0.995029);
        Assert.assertEquals(metrics.getMetrics().get(0).NHOM, 1510547);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFailCreateBafRegressMetricsFile() throws IOException {
        final File outputBase = File.createTempFile("testFailCreateBafRegressMetricsFile", "");
        final File output = new File(outputBase.getAbsolutePath() + "." + CreateBafRegressMetricsFile.FILE_EXTENSION);
        output.deleteOnExit();

        final CreateBafRegressMetricsFile createBafRegressMetricsFile = new CreateBafRegressMetricsFile();
        createBafRegressMetricsFile.INPUT = TEST_BAD_INPUT_FILE;
        createBafRegressMetricsFile.OUTPUT = outputBase;

        Assert.assertEquals(createBafRegressMetricsFile.instanceMain(new String[0]), 1);
    }
}

