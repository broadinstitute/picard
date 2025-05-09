package picard.analysis;

import java.nio.file.Files;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.metrics.MetricsFile;
import picard.cmdline.CommandLineProgramTest;
import picard.util.RExecutor;

public class CollectRrbsMetricsTest extends CommandLineProgramTest {
    
    private static final File TEST_DATA_DIR = new File("testdata/picard");

    public String getCommandLineProgramName() {
        return CollectRrbsMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "sam/CollectRrbsMetrics/input.sam");
        final File ref = new File(TEST_DATA_DIR, "reference/test.fasta");
        final File tempdir = Files.createTempDirectory(getCommandLineProgramName()).toFile();
        final File detailFile   = new File(tempdir, "test.rrbs_detail_metrics");
        final File summaryFile   = new File(tempdir, "test.rrbs_summary_metrics");
        final File pdf   = new File(tempdir, "test.rrbs_qc.pdf");
        final String basename = detailFile.getAbsolutePath().substring(0, detailFile.getAbsolutePath().indexOf(".rrbs_detail_metrics"));
        tempdir.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "M=" + basename,
                "R=" + ref.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RrbsSummaryMetrics, Comparable<?>> summary_output = new MetricsFile<RrbsSummaryMetrics, Comparable<?>>();
        summary_output.read(new FileReader(summaryFile));

        Assert.assertEquals(summary_output.getMetrics().size(), 1);
        final RrbsSummaryMetrics summary = summary_output.getMetrics().get(0);
        Assert.assertEquals(summary.READS_ALIGNED, 6);
        Assert.assertEquals(summary.NON_CPG_BASES, 18);
        Assert.assertEquals(summary.NON_CPG_CONVERTED_BASES, 0);
        Assert.assertEquals(summary.PCT_NON_CPG_BASES_CONVERTED, 0);
        Assert.assertEquals(summary.CPG_BASES_SEEN, 4);
        Assert.assertEquals(summary.CPG_BASES_CONVERTED, 0);
        Assert.assertEquals(summary.PCT_CPG_BASES_CONVERTED, 0);
        Assert.assertEquals(summary.MEAN_CPG_COVERAGE, 4);
        Assert.assertEquals(summary.MEDIAN_CPG_COVERAGE, 4);
        Assert.assertEquals(summary.READS_WITH_NO_CPG, 2);
        Assert.assertEquals(summary.READS_IGNORED_SHORT, 0);
        Assert.assertEquals(summary.READS_IGNORED_MISMATCHES, 0);

        final MetricsFile<RrbsCpgDetailMetrics, Comparable<?>> detail_output = new MetricsFile<RrbsCpgDetailMetrics, Comparable<?>>();
        detail_output.read(new FileReader(detailFile));

        Assert.assertEquals(detail_output.getMetrics().size(), 1);
        final RrbsCpgDetailMetrics detail = detail_output.getMetrics().get(0);
        Assert.assertEquals(detail.SEQUENCE_NAME, "chr1");
        Assert.assertEquals(detail.POSITION, 20);
        Assert.assertEquals(detail.TOTAL_SITES, 4);
        Assert.assertEquals(detail.CONVERTED_SITES, 0);
        Assert.assertEquals(detail.PCT_CONVERTED, 0);

        Assert.assertTrue(pdf.length() > 0);
    }

    @Test
    public void testNoChart() throws IOException {
        final File input = new File(TEST_DATA_DIR, "sam/CollectRrbsMetrics/input.sam");
        final File ref = new File(TEST_DATA_DIR, "reference/test.fasta");
        final File tempdir = Files.createTempDirectory(getCommandLineProgramName()).toFile();
        final File detailFile   = new File(tempdir, "test.rrbs_detail_metrics");
        final File summaryFile   = new File(tempdir, "test.rrbs_summary_metrics");
        final File pdf   = new File(tempdir, "test.rrbs_qc.pdf");
        final String basename = detailFile.getAbsolutePath().substring(0, detailFile.getAbsolutePath().indexOf(".rrbs_detail_metrics"));
        tempdir.deleteOnExit();
        final String[] args = new String[] {
                "I="  + input.getAbsolutePath(),
                "M=" + basename,
                "R=" + ref.getAbsolutePath(),
                "DO_NOT_CREATE_PLOTS=true"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RrbsSummaryMetrics, Comparable<?>> summary_output = new MetricsFile<RrbsSummaryMetrics, Comparable<?>>();
        summary_output.read(new FileReader(summaryFile));

        Assert.assertEquals(summary_output.getMetrics().size(), 1);
        final RrbsSummaryMetrics summary = summary_output.getMetrics().get(0);
        Assert.assertEquals(summary.READS_ALIGNED, 6);
        Assert.assertEquals(summary.NON_CPG_BASES, 18);
        Assert.assertEquals(summary.NON_CPG_CONVERTED_BASES, 0);
        Assert.assertEquals(summary.PCT_NON_CPG_BASES_CONVERTED, 0);
        Assert.assertEquals(summary.CPG_BASES_SEEN, 4);
        Assert.assertEquals(summary.CPG_BASES_CONVERTED, 0);
        Assert.assertEquals(summary.PCT_CPG_BASES_CONVERTED, 0);
        Assert.assertEquals(summary.MEAN_CPG_COVERAGE, 4);
        Assert.assertEquals(summary.MEDIAN_CPG_COVERAGE, 4);
        Assert.assertEquals(summary.READS_WITH_NO_CPG, 2);
        Assert.assertEquals(summary.READS_IGNORED_SHORT, 0);
        Assert.assertEquals(summary.READS_IGNORED_MISMATCHES, 0);

        final MetricsFile<RrbsCpgDetailMetrics, Comparable<?>> detail_output = new MetricsFile<RrbsCpgDetailMetrics, Comparable<?>>();
        detail_output.read(new FileReader(detailFile));

        Assert.assertEquals(detail_output.getMetrics().size(), 1);
        final RrbsCpgDetailMetrics detail = detail_output.getMetrics().get(0);
        Assert.assertEquals(detail.SEQUENCE_NAME, "chr1");
        Assert.assertEquals(detail.POSITION, 20);
        Assert.assertEquals(detail.TOTAL_SITES, 4);
        Assert.assertEquals(detail.CONVERTED_SITES, 0);
        Assert.assertEquals(detail.PCT_CONVERTED, 0);

        Assert.assertEquals(pdf.length(), 0);
    }

    @Test
    public void testFailureGatkLiteDocker() throws IOException {
        final PrintStream stderr = System.err;
        final String gatkLiteDockerProperty = System.getProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR);

        try {
            final ByteArrayOutputStream stderrCapture = new ByteArrayOutputStream();
            System.setErr(new PrintStream(stderrCapture));

            System.setProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR, "true");
            final File input = new File(TEST_DATA_DIR, "sam/CollectRrbsMetrics/input.sam");
            final File ref = new File(TEST_DATA_DIR, "reference/test.fasta");
            final File tempdir = Files.createTempDirectory(getCommandLineProgramName()).toFile();
            final File detailFile   = new File(tempdir, "test.rrbs_detail_metrics");
            final File summaryFile   = new File(tempdir, "test.rrbs_summary_metrics");
            final File pdf   = new File(tempdir, "test.rrbs_qc.pdf");
            final String basename = detailFile.getAbsolutePath().substring(0, detailFile.getAbsolutePath().indexOf(".rrbs_detail_metrics"));
            tempdir.deleteOnExit();
            final String[] args = new String[] {
                    "I="  + input.getAbsolutePath(),
                    "M=" + basename,
                    "R=" + ref.getAbsolutePath()
            };
            Assert.assertEquals(runPicardCommandLine(args), 1);

            Assert.assertEquals(pdf.length(), 0);

            Assert.assertTrue(stderrCapture.toString().contains("The histogram file cannot be written because it requires R, which is not available in the GATK Lite Docker image.")); 
        }
        finally {
            System.setErr(stderr);
            if(gatkLiteDockerProperty != null) {
                System.setProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR, gatkLiteDockerProperty);
            }
            else{
                System.clearProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR);
            } 
        }
    }
}
