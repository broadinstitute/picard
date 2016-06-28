package picard.analysis.artifacts;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class CollectSequencingArtifactMetricsTest extends CommandLineProgramTest {

    private static final File TEST_DIR = new File("testdata/picard/analysis/artifacts/CollectSequencingArtifactMetrics");
    private static final File REFERENCE = new File(TEST_DIR, "test.fasta");
    private static final File TEST_SAM = new File(TEST_DIR, "test.sam");
    private static final File DB_SNP = new File(TEST_DIR, "test.dbsnp.vcf");
    private static final File INTERVALS = new File(TEST_DIR, "test.interval_list");
    private static final File TEST_CASES = new File(TEST_DIR, "ExpectedMetricsOutput");

    private File globalTempOutputDir;

    @Override
    public String getCommandLineProgramName() {
        return CollectSequencingArtifactMetrics.class.getSimpleName();
    }

    @BeforeTest
    public void setUp() throws IOException {
        globalTempOutputDir = IOUtil.createTempDir("artifactMetrics.", ".tmp");
    }

    @AfterTest
    public void tearDown() throws IOException {
        IOUtil.deleteDirectoryTree(globalTempOutputDir);
    }

    /**
     * Run the CLP using standard arguments and maybe some additional ones, then compare to the expected results.
     *
     * @param testCase name of test case (should match one of the file sets in {@code TEST_CASES}
     * @param extraArgs extra arguments of the form {@code "KEY1=VALUE1", "KEY2=VALUE2"}, etc. These can override standard args.
     * @return the base path of the output metrics
     */
    private void runAnalysis(final String testCase, final String ... extraArgs) throws IOException {
        final File actual = new File(globalTempOutputDir, testCase);
        final File expected = new File(TEST_CASES, testCase);

        final Map<String, String> args = new HashMap<String, String>();
        args.put("INPUT", TEST_SAM.getAbsolutePath());
        args.put("OUTPUT", actual.getAbsolutePath());
        args.put("REFERENCE_SEQUENCE", REFERENCE.getAbsolutePath());
        args.put("MINIMUM_INSERT_SIZE", "30"); // test data has this insert size
        args.put("MAXIMUM_INSERT_SIZE", "30"); // test data has this insert size
        args.put("CONTEXT_SIZE", "0"); // ignore context by default, to cut down on file size

        for (final String extraArg : extraArgs) {
            final String[] kv = extraArg.split("=");
            args.put(kv[0], kv[1]);
        }

        runPicardCommandLine(args);
        assertAllFilesEqual(expected, actual);
    }

    private void assertAllFilesEqual(final File expectedBase, final File actualBase) {
        Assert.assertTrue(areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT),"Pre-Adapter summary files differ.");
        Assert.assertTrue(areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT),"Pre-Adapter details files differ.");
        Assert.assertTrue(areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT), "Bait-Bias summary files differ.");
        Assert.assertTrue(areMetricsEqual(expectedBase, actualBase, SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT), "Bait-bias details files differ.");
    }

    private boolean areMetricsEqual(final File expectedBase, final File actualBase, final String extension) {
        return MetricsFile.areMetricsEqual(new File(expectedBase + extension), new File(actualBase + extension));
    }

    @Test
    public void testContext() throws IOException {
        runAnalysis("with_context", "CONTEXT_SIZE=1");
    }

    @Test
    public void testDbSnp() throws IOException {
        runAnalysis("with_dbsnp", "DB_SNP=" + DB_SNP);
    }

    @Test
    public void testIntervalList() throws IOException {
        runAnalysis("with_intervals", "INTERVALS=" + INTERVALS);
    }

    @Test
    public void testNoBqCutoff() throws IOException {
        runAnalysis("no_bq_cutoff", "MINIMUM_QUALITY_SCORE=0");
    }

    @Test
    public void testNoMqCutoff() throws IOException {
        runAnalysis("no_mq_cutoff", "MINIMUM_MAPPING_QUALITY=0");
    }

    @Test
    public void testUnmappedMate() throws IOException {
        runAnalysis("unmapped_mate", "MINIMUM_INSERT_SIZE=0", "MAXIMUM_INSERT_SIZE=0");
    }

}
