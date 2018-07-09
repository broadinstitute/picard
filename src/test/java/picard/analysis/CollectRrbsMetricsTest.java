package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.metrics.MultilevelMetrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class CollectRrbsMetricsTest extends CommandLineProgramTest {
    private static final String CHR_M_SAM = "testdata/picard/metrics/chrMReads.sam";
    private static final String CHR_M_REFERENCE = "testdata/picard/metrics/chrM.reference.fasta";
    private static final String CHR_M_SAM_TEST_UNMAPPED = "testdata/picard/sam/collect_rrbs_metrics_test_unmapped.sam";
    private static final String CHR_M_REFERENCE_TEST_UNMAPPED = "testdata/picard/sam/summary_alignment_stats_test.fasta";

    private File pdf;
    private File summaryOutput;
    private File detailOutput;
    private File tempDir;
    private String prefix;

    @Override
    public String getCommandLineProgramName() {
        return CollectRrbsMetrics.class.getSimpleName();
    }

    @BeforeClass
    private void setUp() throws Exception {
        pdf = File.createTempFile("crmt.", ".pdf");
        pdf.deleteOnExit();
        summaryOutput = File.createTempFile("crmt.", ".rrbs_summary_metrics");
        summaryOutput.deleteOnExit();
        detailOutput = File.createTempFile("crmt.", "detail.sam");
        detailOutput.deleteOnExit();
        tempDir = Files.createTempDirectory("crmt.").toFile();
        prefix = Paths.get(tempDir.toString(), "crmt.").toString();
    }

    @AfterClass
    public void clearTempDir() throws IOException {
        IOUtil.deleteDirectoryTree(tempDir);
    }

    @Test
    public void chrMReads() throws Exception {
        Assert.assertEquals(
                runPicardCommandLine(
                        makeArgList(CHR_M_SAM, null, detailOutput.getAbsolutePath(), pdf.getAbsolutePath(), summaryOutput.getAbsolutePath(), CHR_M_REFERENCE)
                ),
                0
        );

        final RrbsSummaryMetrics metrics = getMultilevelMetrics(summaryOutput);
        Assert.assertEquals(metrics.READS_ALIGNED.intValue(), 5);
        Assert.assertEquals(metrics.NON_CPG_BASES.intValue(), 15);
        Assert.assertEquals(metrics.NON_CPG_CONVERTED_BASES.intValue(), 11);
        Assert.assertEquals(metrics.PCT_NON_CPG_BASES_CONVERTED, 0.733333);
        Assert.assertEquals(metrics.CPG_BASES_SEEN.intValue(), 5);
        Assert.assertEquals(metrics.CPG_BASES_CONVERTED.intValue(), 1);
        Assert.assertEquals(metrics.PCT_CPG_BASES_CONVERTED, 0.2);
        Assert.assertEquals(metrics.MEAN_CPG_COVERAGE, 1.666667);
        Assert.assertEquals(metrics.MEDIAN_CPG_COVERAGE.intValue(), 2);
        Assert.assertEquals(metrics.READS_WITH_NO_CPG.intValue(), 1);
        Assert.assertEquals(metrics.READS_IGNORED_SHORT.intValue(), 1);
        Assert.assertEquals(metrics.READS_IGNORED_MISMATCHES.intValue(), 1);
    }

    @Test
    public void testRrbsCpgDetailMetrics() throws Exception {
        Assert.assertEquals(
                runPicardCommandLine(
                        makeArgList(CHR_M_SAM, null, detailOutput.getAbsolutePath(), pdf.getAbsolutePath(), summaryOutput.getAbsolutePath(), CHR_M_REFERENCE)
                ),
                0
        );

        final RrbsCpgDetailMetrics metricsCpg = getMultilevelMetrics(detailOutput);
        Assert.assertEquals(metricsCpg.SEQUENCE_NAME, "chrM");
        Assert.assertEquals(metricsCpg.POSITION.intValue(), 60);
        Assert.assertEquals(metricsCpg.TOTAL_SITES.intValue(), 1);
        Assert.assertEquals(metricsCpg.CONVERTED_SITES.intValue(), 1);
        Assert.assertEquals(metricsCpg.PCT_CONVERTED.intValue(), 1);
    }

    @Test
    public void testUnmappedReads() throws Exception {
        Assert.assertEquals(
                runPicardCommandLine(
                        makeArgList(CHR_M_SAM_TEST_UNMAPPED, null, detailOutput.getAbsolutePath(), pdf.getAbsolutePath(), summaryOutput.getAbsolutePath(), CHR_M_REFERENCE_TEST_UNMAPPED)
                ),
                0
        );

        final MetricsFile<RrbsCpgDetailMetrics, ?> metricsFile = new MetricsFile<RrbsCpgDetailMetrics, Integer>();
        metricsFile.read(new FileReader(detailOutput));

        Assert.assertEquals(metricsFile.getMetrics().size(), 2); // this metric skips unmapped reads, so we get 2 reads instead of 3
    }

    @Test
    public void testRrbsCpgDetailMetricsByPrefix() throws Exception {
        Assert.assertEquals(
                runPicardCommandLine(
                        makeArgList(CHR_M_SAM, prefix, null, null, null, CHR_M_REFERENCE)
                ),
                0
        );

        final String output = prefix + CollectRrbsMetrics.DETAIL_FILE_EXTENSION;
        final RrbsCpgDetailMetrics metricsCpg = getMultilevelMetrics(output);
        Assert.assertEquals(metricsCpg.SEQUENCE_NAME, "chrM");
        Assert.assertEquals(metricsCpg.POSITION.intValue(), 60);
        Assert.assertEquals(metricsCpg.TOTAL_SITES.intValue(), 1);
        Assert.assertEquals(metricsCpg.CONVERTED_SITES.intValue(), 1);
        Assert.assertEquals(metricsCpg.PCT_CONVERTED.intValue(), 1);
    }

    private <T extends MultilevelMetrics> T getMultilevelMetrics(final File file) throws FileNotFoundException {
        final MetricsFile<T, ?> metricsFile = new MetricsFile<T, Integer>();
        metricsFile.read(new FileReader(file));
        return metricsFile.getMetrics().get(0);
    }

    private <T extends MultilevelMetrics> T getMultilevelMetrics(final String fileName) throws FileNotFoundException {
        return getMultilevelMetrics(new File(fileName));
    }

    @Test(dataProvider = "incorrectArguments", expectedExceptions = IllegalArgumentException.class)
    public void checkingArguments(final String prefix, final String detailOutput, final String chartOutput, final String summaryOutput) {
        runPicardCommandLine(makeArgList(CHR_M_SAM, prefix, detailOutput, chartOutput, summaryOutput, CHR_M_REFERENCE));
    }

    @DataProvider(name = "incorrectArguments")
    public Object[][] makeIncorrectArguments() {
        return new Object[][]{
                // METRICS_FILE_PREFIX  OUTPUT  CHART_OUTPUT    SUMMARY_OUTPUT
                {prefix, detailOutput, null, null},
                {prefix, null, pdf, null},
                {prefix, null, null, summaryOutput},
                {null, detailOutput, pdf, null},
                {null, detailOutput, null, summaryOutput},
                {null, null, pdf, summaryOutput},
        };
    }

    private List<String> makeArgList(final String input, final String prefix,
                                     final String detailOutput, final String chartOutput, final String summaryOutput,
                                     final String referenceSequence) {
        List<String> args = new ArrayList<>();
        if (input != null) {
            args.add("INPUT=" + input);
        }
        if (prefix != null) {
            args.add("METRICS_FILE_PREFIX=" + prefix);
        }
        if (detailOutput != null) {
            args.add("OUTPUT=" + detailOutput);
        }
        if (chartOutput != null) {
            args.add("CHART_OUTPUT=" + chartOutput);
        }
        if (summaryOutput != null) {
            args.add("SUMMARY_OUTPUT=" + summaryOutput);
        }
        if (referenceSequence != null) {
            args.add("R=" + referenceSequence);
        }
        return args;
    }
}