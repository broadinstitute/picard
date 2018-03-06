package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.*;
import picard.vcf.SamTestUtils;

import java.io.*;
import java.lang.reflect.Field;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Stream;

import static picard.fingerprint.FingerprintIdDetails.multipleValuesString;

/**
 * Tests for CrosscheckReadGroupFingerprints
 */
public class CrosscheckReadGroupFingerprintsTest {

    private static final File TEST_DIR = new File("testdata/picard/fingerprint/");
    private static final File HAPLOTYPE_MAP = new File(TEST_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    private static final File NA12891_r1_sam = new File(TEST_DIR, "NA12891.over.fingerprints.r1.sam");
    private static final File NA12891_r2_sam = new File(TEST_DIR, "NA12891.over.fingerprints.r2.sam");

    //this is a copy of a previous one, but with a different sample name
    private static final File NA12891_named_NA12892_r1_sam = new File(TEST_DIR, "NA12891_named_NA12892.over.fingerprints.r1.sam");

    private static final File NA12892_r1_sam = new File(TEST_DIR, "NA12892.over.fingerprints.r1.sam");
    private static final File NA12892_r2_sam = new File(TEST_DIR, "NA12892.over.fingerprints.r2.sam");

    private static File NA12891_r1, NA12891_r2, NA12891_named_NA12892_r1, NA12892_r1, NA12892_r2;

    private static final int NA12891_r1_RGs = 27;
    private static final int NA12891_r2_RGs = 26;
    private static final int NA12892_r1_RGs = 25;
    private static final int NA12892_r2_RGs = 26;

    private final Map<CrosscheckMetric.DataType, List<String>> lookupMap = new HashMap<>(4);

    @BeforeClass
    public void setup() throws IOException {
        NA12891_r1 = SamTestUtils.createIndexedBam(NA12891_r1_sam, NA12891_r1_sam);
        NA12891_r2 = SamTestUtils.createIndexedBam(NA12891_r2_sam, NA12891_r2_sam);
        NA12891_named_NA12892_r1 = SamTestUtils.createIndexedBam(NA12891_named_NA12892_r1_sam, NA12891_named_NA12892_r1_sam);
        NA12892_r1 = SamTestUtils.createIndexedBam(NA12892_r1_sam, NA12892_r1_sam);
        NA12892_r2 = SamTestUtils.createIndexedBam(NA12892_r2_sam, NA12892_r2_sam);

        lookupMap.put(CrosscheckMetric.DataType.FILE, new ArrayList<>(Arrays.asList("LEFT_FILE", "RIGHT_FILE")));

        lookupMap.put(CrosscheckMetric.DataType.SAMPLE, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.SAMPLE).addAll(Arrays.asList("LEFT_SAMPLE", "RIGHT_SAMPLE"));
        lookupMap.get(CrosscheckMetric.DataType.SAMPLE).addAll(lookupMap.get(CrosscheckMetric.DataType.FILE));

        lookupMap.put(CrosscheckMetric.DataType.LIBRARY, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.LIBRARY).addAll(Arrays.asList("LEFT_LIBRARY", "RIGHT_LIBRARY"));
        lookupMap.get(CrosscheckMetric.DataType.LIBRARY).addAll(lookupMap.get(CrosscheckMetric.DataType.SAMPLE));

        lookupMap.put(CrosscheckMetric.DataType.READGROUP, new ArrayList<>());
        lookupMap.get(CrosscheckMetric.DataType.READGROUP).addAll(Arrays.asList("LEFT_RUN_BARCODE", "LEFT_LANE",
                "LEFT_MOLECULAR_BARCODE_SEQUENCE","RIGHT_RUN_BARCODE",
                "RIGHT_LANE", "RIGHT_MOLECULAR_BARCODE_SEQUENCE"));
        lookupMap.get(CrosscheckMetric.DataType.READGROUP).addAll(lookupMap.get(CrosscheckMetric.DataType.LIBRARY));
    }

    @DataProvider(name = "bamFilesRGs")
    public Object[][] bamFilesRGs() {
        return new Object[][]{
                {NA12891_r1, NA12891_r2, false, 0, NA12891_r1_RGs + NA12891_r2_RGs },
                {NA12891_r1, NA12892_r1, false, 0, NA12891_r1_RGs + NA12892_r1_RGs },
                {NA12891_r1, NA12892_r2, false, 0, NA12891_r1_RGs + NA12892_r2_RGs },
                {NA12892_r1, NA12892_r2, false, 0, NA12892_r1_RGs + NA12892_r2_RGs },
                {NA12892_r2, NA12891_r2, false, 0, NA12892_r2_RGs + NA12891_r2_RGs },
                {NA12892_r2, NA12891_r1, false, 0, NA12892_r2_RGs + NA12891_r1_RGs },
                {NA12891_r1, NA12891_r2, true, 0, NA12891_r1_RGs + NA12891_r2_RGs },
                {NA12891_r1, NA12892_r1, true, 1, NA12891_r1_RGs + NA12892_r1_RGs },
                {NA12891_r1, NA12892_r2, true, 1, NA12891_r1_RGs + NA12892_r2_RGs },
                {NA12892_r1, NA12892_r2, true, 0, NA12892_r1_RGs + NA12892_r2_RGs },
                {NA12892_r2, NA12891_r2, true, 1, NA12892_r2_RGs + NA12891_r2_RGs },
                {NA12892_r2, NA12891_r1, true, 1, NA12892_r2_RGs + NA12891_r1_RGs}
        };
    }

    @Test(dataProvider = "bamFilesRGs")
    public void testCrossCheckRGs(final File file1, final File file2, final boolean expectAllMatch, final int expectedRetVal, final int expectedNMetrics) throws IOException, NoSuchFieldException {

        File metrics = File.createTempFile("Fingerprinting", "NA1291.RG.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -2.0,
                "EXPECT_ALL_GROUPS_TO_MATCH=" + expectAllMatch
        };

        doTest(args, metrics, expectedRetVal, expectedNMetrics * expectedNMetrics , CrosscheckMetric.DataType.READGROUP, expectAllMatch);
    }

    @DataProvider(name = "bamFilesLBs")
    public Object[][] bamFilesLBs() {

        return new Object[][]{
                {NA12891_r1, NA12891_r2, 0, true},
                {NA12891_r1, NA12892_r1, 0, false},
                {NA12892_r2, NA12891_r2, 0, false},
                {NA12892_r2, NA12891_r1, 0, false},
                {NA12891_r1, NA12891_named_NA12892_r1_sam, 1, true},
                {NA12892_r1, NA12891_named_NA12892_r1_sam, 1, false},
                {NA12891_r2, NA12891_named_NA12892_r1_sam, 1, true},
                {NA12892_r2, NA12891_named_NA12892_r1_sam, 1, false},
        };
    }

    @Test(dataProvider = "bamFilesLBs")
    public void testCrossCheckLBs(final File file1, final File file2, final int expectedRetVal, final boolean expectAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.LB.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_LIBRARIES=true"
        };
        final int numLibs = 2;
        doMatrixTest(args, metrics, expectedRetVal, numLibs, expectAllMatch);
    }

    @DataProvider(name = "bamFilesSources")
    public Object[][] bamFilesSources() {

        return new Object[][]{
                {NA12891_r1, NA12891_r2, 0},
                {NA12892_r1, NA12892_r2, 0},
                {NA12891_r1, NA12892_r1, 0},
                {NA12892_r2, NA12891_r2, 0},
                {NA12892_r2, NA12891_named_NA12892_r1, 1}, // the two files contain different samples one has wrong name
                {NA12891_r2, NA12891_named_NA12892_r1, 1}, // unexpected match
                {NA12892_r1, NA12891_named_NA12892_r1, 1}, // the two files contain different samples one has wrong name
                {NA12891_r1, NA12891_named_NA12892_r1, 1}, // unexpected match
        };
    }

    @DataProvider(name = "bamFilesSMs")
    public Object[][] bamFilesSMs() {

        return new Object[][]{
                {NA12891_r1, NA12891_r2, 0, 1, true},
                {NA12891_r1, NA12892_r1, 0, 2, false},
                {NA12892_r2, NA12891_r2, 0, 2, false},
                {NA12892_r2, NA12891_named_NA12892_r1, 0, 1, true}, // no error since only one sample in aggregate
                {NA12891_r2, NA12891_named_NA12892_r1, 1, 2, true}, // unexpected match
                {NA12892_r1, NA12891_named_NA12892_r1, 0, 1, true}, // no error since only one sample in aggregate
                {NA12891_r1, NA12891_named_NA12892_r1, 1, 2, true}, // unexpected match
        };
    }

    @Test(dataProvider = "bamFilesSMs")
    public void testCrossCheckSMs(final File file1, final File file2, final int expectedRetVal, final int numberOfSamples, final boolean expectedAllMatch) throws IOException {
        File metrics = File.createTempFile("Fingerprinting", "NA1291.SM.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_SAMPLES=true"
        };
        doMatrixTest(args, metrics, expectedRetVal, numberOfSamples, expectedAllMatch);
    }

    private void doTest(final String[] args, final File metrics, final int expectedRetVal, final int expectedNMetrics, final CrosscheckMetric.DataType expectedType, final boolean expectAllMatch) throws IOException, NoSuchFieldException {

        final CrosscheckReadGroupFingerprints crossChecker = new CrosscheckReadGroupFingerprints();
        Assert.assertEquals(crossChecker.instanceMain(args), expectedRetVal);

        final MetricsFile<CrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(metrics));

        Assert.assertTrue(metricsOutput.getMetrics().stream()
                .allMatch(m -> m.DATA_TYPE == expectedType));

        Assert.assertTrue(metricsOutput.getMetrics().stream()
                .allMatch(m -> m.LOD_SCORE_NORMAL_TUMOR != null));
        Assert.assertTrue(metricsOutput.getMetrics().stream()
                .allMatch(m -> m.LOD_SCORE != null));
        Assert.assertTrue(metricsOutput.getMetrics().stream()
                .allMatch(m -> m.LOD_SCORE_TUMOR_NORMAL != null));

        if (expectAllMatch) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isMatch() == m.LEFT_SAMPLE.equals(m.RIGHT_SAMPLE)));
        } else if (expectedRetVal == 0) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .allMatch(m -> m.RESULT == CrosscheckMetric.FingerprintResult.INCONCLUSIVE ||
                            m.RESULT.isExpected()));
        } else {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> !m.RESULT.isExpected()));
        }

        Assert.assertEquals(metricsOutput.getMetrics().size(), expectedNMetrics);

        // at the readgroup level things aren't always conclusive...
        if (!metricsOutput.getMetrics().isEmpty() && expectedType != CrosscheckMetric.DataType.READGROUP) {
            Assert.assertTrue(metricsOutput.getMetrics().stream()
                    .anyMatch(m -> m.RESULT != CrosscheckMetric.FingerprintResult.INCONCLUSIVE));
        }

        //check that fields that should have an actual value, indeed do
        for (final String fieldName : lookupMap.get(expectedType)) {
            final Field field = CrosscheckMetric.class.getField(fieldName);
            Assert.assertTrue(metricsOutput.getMetrics().stream().allMatch(m -> {
                try {
                    return field.get(m) != multipleValuesString && field.get(m) != null;
                } catch (IllegalAccessException e) {
                    e.printStackTrace();
                    return false;
                }
            }));
        }
    }

    private void doMatrixTest(final String[] args, final File metrics, final int expectedRetVal, final int expectedNMetrics, final boolean expectAllMatch) throws IOException {

        final CrosscheckReadGroupFingerprints crossChecker = new CrosscheckReadGroupFingerprints();
        Assert.assertEquals(crossChecker.instanceMain(args), expectedRetVal);
        Assert.assertTrue(metrics.canRead());

        try (Stream<String> lines = Files.lines(metrics.toPath(), Charset.defaultCharset())) {
            long numOfLines = lines.count() - 1; //subtract one to account for the column header
            Assert.assertEquals(numOfLines, expectedNMetrics);
        }

        if (expectAllMatch) {
            try (Stream<String> lines = Files.lines(metrics.toPath(), Charset.defaultCharset())) {
                lines.skip(1).forEach(s -> {
                    final List<String> strings = Arrays.asList(s.split("\t"));
                    strings.subList(1, strings.size()).forEach(str -> Assert.assertTrue(Double.parseDouble(str) > 0D, "expected positive value, found: " + str));
                });
            }
        }
    }

    @Test
    public void canWriteToDevNull() throws IOException {
        File f = new File("/dev/null");
        Assert.assertTrue(f.canRead());

        final OutputStream stream = new FileOutputStream(f);
        final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(stream));

        writer.write("Just a test");
        writer.close();

    }

    @DataProvider(name = "newParametersData")
    public Object[][] newParametersData() {

        return new Object[][]{
                {"CROSSCHECK_BY=" + CrosscheckMetric.DataType.LIBRARY},
                {"CROSSCHECK_BY=" + CrosscheckMetric.DataType.SAMPLE},
                {"CROSSCHECK_BY=" + CrosscheckMetric.DataType.FILE},
                {"MATRIX_OUTPUT=matrix.file"},
                {"SECOND_INPUT="+"/dev/stdin"}
        };
    }

    @Test(dataProvider = "newParametersData")
    public void testCannotUseNewParameters(final String extraParameter) {

        final File file1 = NA12891_r1;
        final File file2 = NA12892_r1;

        final String[] args = new String[]{
                "INPUT=" + file1.getAbsolutePath(),
                "INPUT=" + file2.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -1.0,
                "CROSSCHECK_SAMPLES=true",
                extraParameter
        };

        final CrosscheckReadGroupFingerprints crossChecker = new CrosscheckReadGroupFingerprints();

        Assert.assertEquals(crossChecker.instanceMain(args), 1);
    }
}
