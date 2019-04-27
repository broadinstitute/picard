package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class CalculateFingerprintMetricsTest extends CommandLineProgramTest {

    @Override
    public String getCommandLineProgramName() {
        return CalculateFingerprintMetrics.class.getSimpleName();
    }

    private static final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private static final File TEST_GENOTYPES_VCF1 = new File(TEST_DATA_DIR, "NA12892.g.vcf");
    private static final File NA12891_r1_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
    private static final File NA12892_r1_sam = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");
    private static final File na12891_fp = new File(TEST_DATA_DIR, "NA12891.fp.vcf");
    private static final File na12892_fp = new File(TEST_DATA_DIR, "NA12892.fp.vcf");

    private static final File HAPLOTYPE_MAP = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    @DataProvider
    public static Object[][] testFingerprintMetricData() {
        return new Object[][]{
                {TEST_GENOTYPES_VCF1},
                {NA12891_r1_sam},
                {NA12892_r1_sam},
                {na12891_fp},
                {na12892_fp}
        };
    }

    @Test(dataProvider = "testFingerprintMetricData")
    public void testFingerprintMetric(final File fingerprintFile) throws IOException {
        final File outputFile = File.createTempFile("testFingerprintMetric", "fingerprint_metric");
        outputFile.deleteOnExit();

        runIt(fingerprintFile, outputFile);
    }

    private void runIt(final File inputFile, final File outputFile) {
        final List<String> args = new ArrayList<>(Arrays.asList(
                "INPUT=" + inputFile.getAbsolutePath(),
                "OUTPUT=" + outputFile.getAbsolutePath(),
                "H=" + HAPLOTYPE_MAP.getAbsolutePath()));
        assertEquals(runPicardCommandLine(args),0);

        final List<FingerprintMetrics> metrics = MetricsFile.readBeans(outputFile);

        for (FingerprintMetrics metric : metrics) {
            Assert.assertTrue(metric.CHI_SQUARED_PVALUE > 0.05);
            Assert.assertTrue(metric.LOG10_CHI_SQUARED_PVALUE > -2);
            Assert.assertTrue(metric.LOD_SELF_CHECK > 0);
            Assert.assertTrue(metric.CROSS_ENTROPY_LOD < 1);
            Assert.assertTrue(metric.DEFINITE_GENOTYPES >= 0);
            Assert.assertTrue(metric.HAPLOTYPES > 4);
            Assert.assertTrue(metric.HAPLOTYPES_WITH_EVIDENCE > 0);
            Assert.assertTrue(metric.HET_CHI_SQUARED_PVALUE > 0.05);
            Assert.assertTrue(metric.LOG10_HET_CHI_SQUARED_PVALUE > -2);
            Assert.assertTrue(metric.HET_CROSS_ENTROPY_LOD < 1);
            Assert.assertTrue(metric.HOM_CHI_SQUARED_PVALUE > 0.05);
            Assert.assertTrue(metric.LOG10_HOM_CHI_SQUARED_PVALUE > -2);
            Assert.assertTrue(metric.HOM_CROSS_ENTROPY_LOD < 1);
            Assert.assertTrue(metric.NUM_HET >= 0);
            Assert.assertTrue(metric.NUM_HOM_ALLELE1 >= 0);
            Assert.assertTrue(metric.NUM_HOM_ALLELE2 >= 0);
            Assert.assertTrue(metric.EXPECTED_HET >= 0);
            Assert.assertTrue(metric.EXPECTED_HOM_ALLELE1>= 0);
            Assert.assertTrue(metric.EXPECTED_HOM_ALLELE2 >= 0);
            Assert.assertTrue(metric.DISCRIMINATORY_POWER > 0);
        }
    }
}