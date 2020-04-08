package picard.fingerprint;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.TribbleException.MalformedFeatureFile;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.TimeUnit;

public class CheckFingerprintTest extends CommandLineProgramTest {

    private final File tempFolder = new File(TEST_DATA_DIR + "/tempCheckFPDir/");
    private static final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");
    private static final File SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING =
            new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");
    private static final File SHIFTED_HAPLOTYPE_DATABASE_FOR_TESTING =
            new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.shifted.for.crams.txt");
  private static final File SHIFTED_REFERENCE =
            new File(TEST_DATA_DIR, "reference.shifted.for.crams.fasta");
    private static final String TEST_INPUT_VCF1 = new File(TEST_DATA_DIR, "NA12892.vcf").getAbsolutePath();
    private static final String TEST_INPUT_VCF_EMPTY = new File(TEST_DATA_DIR, "/tempCheckFPDir/void.vcf").getAbsolutePath();
    private static final String TEST_INPUT_VCF_NO_FILE = new File(TEST_DATA_DIR, "noFile.vcf").getAbsolutePath();
    private static final String TEST_OUTPUT = new File(TEST_DATA_DIR + "/tempCheckFPDir/otp.fp").getAbsolutePath();
    private static final String TEST_GENOTYPES_VCF1 = new File(TEST_DATA_DIR, "NA12892.g.vcf").getAbsolutePath();
    private static final String TEST_GENOTYPES_VCF_NO_FILE = new File(TEST_DATA_DIR, "noFile.g.vcf").getAbsolutePath();
    private static final File RESULT_EXAMPLE_SUMMARY =
            new File(TEST_DATA_DIR, "fingerprinting_summary_metrics.example");
    private static final File RESULT_EXAMPLE_DETAIL =
            new File(TEST_DATA_DIR, "fingerprinting_detail_metrics.example");

    private static final File HAPLOTYPE_MAP = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    private static final File NA12891_r1_sam = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r1.sam");
    private static final File NA12892_r1_sam = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r1.sam");

    private static File NA12892_1_vcf;
    private static File NA12891_named_NA12892_vcf;

    private static final File na12891_fp = new File(TEST_DATA_DIR, "NA12891.fp.vcf");
    private static final File na12892_fp = new File(TEST_DATA_DIR, "NA12892.fp.vcf");

    @BeforeClass
    private void mkTemp() {
        if (!tempFolder.exists()) {
            tempFolder.mkdir();
        }
    }

    @AfterClass
    private void rmTemp() {
        IOUtil.deleteDirectoryTree(tempFolder);
    }

    @DataProvider(name = "badData")
    public Object[][] badData() {
        return new Object[][]{
                {TEST_INPUT_VCF_EMPTY, TEST_OUTPUT,
                        TEST_GENOTYPES_VCF1, SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING},
                {TEST_INPUT_VCF_NO_FILE, TEST_OUTPUT,
                        TEST_GENOTYPES_VCF1, SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING},
                {TEST_INPUT_VCF1, TEST_OUTPUT,
                        TEST_GENOTYPES_VCF_NO_FILE, SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING},
                {TEST_INPUT_VCF_NO_FILE, TEST_OUTPUT,
                        TEST_GENOTYPES_VCF_NO_FILE, SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING},
        };
    }

    @Test
    public void testBaseOutput() {
        String[] args = new String[]{
                "I=" + TEST_INPUT_VCF1,
                "O=" + TEST_OUTPUT,
                "G=" + TEST_GENOTYPES_VCF1,
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);
        Assert.assertTrue(MetricsFile.areMetricsEqual(new File(TEST_OUTPUT + ".fingerprinting_summary_metrics"),
                RESULT_EXAMPLE_SUMMARY));
        Assert.assertTrue(MetricsFile.areMetricsEqual(new File(TEST_OUTPUT + ".fingerprinting_detail_metrics"),
                RESULT_EXAMPLE_DETAIL));
    }

    @Test
    public void testMismatchingSamples() {
        String[] args = new String[]{
                "I=" + NA12891_r1_sam,
                "O=" + TEST_OUTPUT,
                "G=" + TEST_GENOTYPES_VCF1,
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING
        };

        Assert.assertEquals(runPicardCommandLine(args), 1);
    }

    @Test
    public void testMismatchingSamples2() {
        String[] args = new String[]{
                "I=" + TEST_INPUT_VCF1,
                "O=" + TEST_OUTPUT,
                "G=" + TEST_GENOTYPES_VCF1,
                "EXPECTED_SAMPLE_ALIAS=TEST123",
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING
        };

        Assert.assertEquals(runPicardCommandLine(args), 1);
    }

    @Test
    public void testSummaryAndDetailOutputs() {
        String[] args = new String[]{
                "I=" + TEST_INPUT_VCF1,
                "S=" + TEST_DATA_DIR + "/tempCheckFPDir/summary",
                "D=" + TEST_DATA_DIR + "/tempCheckFPDir/detail",
                "G=" + TEST_GENOTYPES_VCF1,
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        Assert.assertTrue(MetricsFile.areMetricsEqual(new File(TEST_DATA_DIR + "/tempCheckFPDir/summary"),
                RESULT_EXAMPLE_SUMMARY));
        Assert.assertTrue(MetricsFile.areMetricsEqual(new File(TEST_DATA_DIR + "/tempCheckFPDir/detail"),
                RESULT_EXAMPLE_DETAIL));
    }

    @Test(dataProvider = "badData", expectedExceptions = {MalformedFeatureFile.class, SAMException.class})
    public void testBadData(final String inputVcf,
                            final String outputLoc,
                            final String genotypesFile,
                            final File haplotypeFile) {
        String[] args = new String[]{
                "I=" + inputVcf,
                "O=" + outputLoc,
                "G=" + genotypesFile,
                "H=" + haplotypeFile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
    }

    @Test
    public void testIdentifyContaminant() {

        final File mixture = new File(TEST_DATA_DIR, "NA128791_in_NA12892.25_PCT.sam");
        final File contaminant = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.r2.sam");
        final File contaminated = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.r2.sam");

        final double contamAmount = .25;
        final String sample = "NA12892";
        final FingerprintChecker checker = new FingerprintChecker(SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING);

        checker.setLocusMaxReads(100);
        final Fingerprint fpContaminant = checker.identifyContaminant(mixture.toPath(), contamAmount).get(sample);
        Assert.assertNotNull(fpContaminant);

        final Fingerprint fpContamination = checker.fingerprintFiles(Collections.singleton(contaminant.toPath()), 1, 1, TimeUnit.DAYS)
                .values().stream()
                .reduce((a, b) -> {
                    a.merge(b);
                    return a;
                })
                .orElseThrow(() -> new IllegalArgumentException("Did not find any data for contaminant"));
        Assert.assertNotNull(fpContamination);

        final Fingerprint fpContaminated = checker.fingerprintFiles(Collections.singleton(contaminated.toPath()), 1, 1, TimeUnit.DAYS)
                .values().stream()
                .reduce((a, b) -> {
                    a.merge(b);
                    return a;
                })
                .orElseThrow(() -> new IllegalArgumentException("Did not find any data for contaminated"));
        Assert.assertNotNull(fpContaminated);

        Assert.assertTrue(FingerprintChecker.calculateMatchResults(fpContaminant, fpContamination).getLOD() > 1D);
        Assert.assertTrue(FingerprintChecker.calculateMatchResults(fpContaminant, fpContaminated).getLOD() < -4D);
    }

    @Test
    public void testFPToVC() throws IOException {

        final File Na12892 = new File(TEST_DATA_DIR, "NA12892.over.fingerprints.shifted.for.crams.r1.sam");
        final File Na12891 = new File(TEST_DATA_DIR, "NA12891.over.fingerprints.shifted.for.crams.r1.sam");

       final FingerprintChecker checker = new FingerprintChecker(SHIFTED_HAPLOTYPE_DATABASE_FOR_TESTING);

        checker.setLocusMaxReads(100);
        final Fingerprint fingerprint = checker.identifyContaminant(Na12892.toPath(), 1).get("NA12892");
        Assert.assertNotNull(fingerprint);

        final File tempFile = File.createTempFile("testWriteFingerprint",".vcf");
        tempFile.deleteOnExit();

        FingerprintUtils.writeFingerPrint(fingerprint,tempFile,SHIFTED_REFERENCE,
                "NA12892",null);

        final Fingerprint NA12892FromVCF = checker.fingerprintFiles(Collections.singleton(tempFile.toPath()), 1, 1, TimeUnit.DAYS)
                .values().stream()
                .reduce((a, b) -> {
                    a.merge(b);
                    return a;
                })
                .orElseThrow(() -> new IllegalArgumentException("Did not find any data for contaminant"));
        Assert.assertNotNull(NA12892FromVCF);

        Assert.assertTrue(FingerprintChecker.calculateMatchResults(fingerprint, NA12892FromVCF).getLOD() > 1D);


        final Fingerprint NA12891_fp = checker.fingerprintFiles(Collections.singleton(Na12891.toPath()), 1, 1, TimeUnit.DAYS)
                .values().stream()
                .reduce((a, b) -> {
                    a.merge(b);
                    return a;
                })
                .orElseThrow(() -> new IllegalArgumentException("Did not find any data for contaminated"));
        Assert.assertNotNull(NA12891_fp);

        Assert.assertTrue(FingerprintChecker.calculateMatchResults(NA12891_fp, NA12892FromVCF).getLOD() < -4D);

        final Fingerprint NA12892_fp = checker.fingerprintFiles(Collections.singleton(Na12892.toPath()), 1, 1, TimeUnit.DAYS)
                .values().stream()
                .reduce((a, b) -> {
                    a.merge(b);
                    return a;
                })
                .orElseThrow(() -> new IllegalArgumentException("Did not find any data for contaminated"));
        Assert.assertNotNull(NA12892_fp);

        Assert.assertTrue(FingerprintChecker.calculateMatchResults(NA12892FromVCF, NA12892_fp).getLOD() >1D);
    }




    @Override
    public String getCommandLineProgramName() {
        return CheckFingerprint.class.getSimpleName();
    }

    @BeforeClass
    public void setup() throws IOException {
        NA12891_named_NA12892_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12891_named_NA12892.vcf"), "fingerprint");
        NA12892_1_vcf = VcfTestUtils.createTemporaryIndexedVcfFromInput(new File(TEST_DATA_DIR, "NA12892.vcf"), "fingerprint");
    }

    @DataProvider(name = "samsToFingerprint")
    Object[][] samsToFingerprint() {
        return new Object[][]{
                {NA12891_r1_sam, na12891_fp, 0},
                {NA12892_r1_sam, na12891_fp, 1},
        };
    }

    @DataProvider(name = "vcfsToFingerprint")
    Object[][] vcfsToFingerprint() {
        return new Object[][]{
                {NA12891_named_NA12892_vcf, na12892_fp, 0},
                {NA12892_1_vcf, na12892_fp, 0},
        };
    }

    @Test(dataProvider = "samsToFingerprint")
    void testCheckFingerprintSam(final File file, final File genotypes, final int expectedRetVal) throws IOException {
        tester(false, file, genotypes, expectedRetVal);
    }

    @Test(dataProvider = "vcfsToFingerprint")
    void testCheckFingerprintVcf(final File file, final File genotypes, final int expectedRetVal) throws IOException {
        tester(false, file, genotypes, expectedRetVal);
    }

    @Test(dataProvider = "samsToFingerprint")
    void testCheckFingerprintNoRg(final File file, final File genotypes, final int expectedRetVal) throws IOException {
        tester(true, file, genotypes, expectedRetVal);
    }

    private File tester(boolean ignoreRG, final File file, final File genotypes, final int expectedRetVal) throws IOException {
        final List<String> args = new ArrayList<>();
        final File outputSummary = File.createTempFile("fingerprint", "summary_metrics");
        outputSummary.deleteOnExit();
        final File outputDetail = File.createTempFile("fingerprint", "detail_metrics");
        outputSummary.deleteOnExit();

        args.add("INPUT=" + file.getAbsolutePath());
        args.add("G=" + genotypes.getAbsolutePath());
        if (ignoreRG) args.add("IGNORE_RG=true");
        args.add("H=" + HAPLOTYPE_MAP.getAbsolutePath());
        args.add("SUMMARY_OUTPUT=" + outputSummary.getAbsolutePath());
        args.add("DETAIL_OUTPUT=" + outputDetail.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args), expectedRetVal);

        Assert.assertTrue(outputSummary.exists(), "Expected output file " + outputSummary.getAbsolutePath() + " to exist.");
        Assert.assertTrue(outputDetail.exists(), "Expected output file " + outputDetail.getAbsolutePath() + " to exist.");

        return outputSummary;
    }
}