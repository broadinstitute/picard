package picard.arrays;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CreateVerifyIDIntensityContaminationMetricsFileTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/arrays");
    private static final File TEST_INPUT_FILE = new File(TEST_DATA_DIR, "VerifyIDIntensity.txt");
    private static final File TEST_BAD_INPUT_FILE = new File(TEST_DATA_DIR, "input.vcf");

    @Test
    public void testCreateVerifyIDIntensityContaminationMetricsFile() throws IOException {
        final File outputBase = File.createTempFile("testCreateVerifyIDIntensityContaminationMetricsFile", "");
        final File output = new File(outputBase.getAbsolutePath() + "." + CreateVerifyIDIntensityContaminationMetricsFile.FILE_EXTENSION);
        output.deleteOnExit();
        System.out.println(outputBase.getAbsolutePath());
        System.out.println(output.getAbsolutePath());

        final CreateVerifyIDIntensityContaminationMetricsFile createVerifyIDIntensityContaminationMetricsFile = new CreateVerifyIDIntensityContaminationMetricsFile();
        createVerifyIDIntensityContaminationMetricsFile.INPUT = TEST_INPUT_FILE;
        createVerifyIDIntensityContaminationMetricsFile.OUTPUT = outputBase;

        Assert.assertEquals(createVerifyIDIntensityContaminationMetricsFile.instanceMain(new String[0]), 0);

        final MetricsFile<VerifyIDIntensityContaminationMetrics, Comparable<?>> metrics = new MetricsFile<>();
        metrics.read(new FileReader(output));

        Assert.assertEquals(metrics.getMetrics().size(), 2);
        Assert.assertEquals(metrics.getMetrics().get(0).ID, 0);
        Assert.assertEquals(metrics.getMetrics().get(0).PCT_MIX, 0.214766);
        Assert.assertEquals(metrics.getMetrics().get(0).LLK, 157575);
        Assert.assertEquals(metrics.getMetrics().get(0).LLK0, 177169);

        Assert.assertEquals(metrics.getMetrics().get(1).ID, 1);
        Assert.assertEquals(metrics.getMetrics().get(1).PCT_MIX, 0.214767);
        Assert.assertEquals(metrics.getMetrics().get(1).LLK, 157576);
        Assert.assertEquals(metrics.getMetrics().get(1).LLK0, 177170);
    }

    @Test(expectedExceptions = PicardException.class)
    public void testFailCreateVerifyIDIntensityContaminationMetricsFile() throws IOException {
        final File outputBase = File.createTempFile("testFailCreateVerifyIDIntensityContaminationMetricsFile", "");
        final File output = new File(outputBase.getAbsolutePath() + "." + CreateVerifyIDIntensityContaminationMetricsFile.FILE_EXTENSION);
        output.deleteOnExit();
        System.out.println(outputBase.getAbsolutePath());
        System.out.println(output.getAbsolutePath());

        final CreateVerifyIDIntensityContaminationMetricsFile createVerifyIDIntensityContaminationMetricsFile = new CreateVerifyIDIntensityContaminationMetricsFile();
        createVerifyIDIntensityContaminationMetricsFile.INPUT = TEST_BAD_INPUT_FILE;
        createVerifyIDIntensityContaminationMetricsFile.OUTPUT = outputBase;

        Assert.assertEquals(createVerifyIDIntensityContaminationMetricsFile.instanceMain(new String[0]), 1);
    }
}
