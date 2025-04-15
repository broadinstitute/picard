package picard.analysis;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.metrics.MetricsFile;
import picard.cmdline.CommandLineProgramTest;

public class CollectBaseDistributionByCycleTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/BaseDistributionByCycle/");

    public String getCommandLineProgramName() {
        return CollectBaseDistributionByCycle.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "input.sam");
        final File outfile   = File.createTempFile("test", ".basedistributionbycycle");
        final File pdf   = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
    
        final String[] args = new String[] {
            "INPUT="  + input.getAbsolutePath(),
            "OUTPUT=" + outfile.getAbsolutePath(),
            "CHART=" + pdf.getAbsolutePath()
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<BaseDistributionByCycleMetrics, Comparable<?>> output = new MetricsFile<BaseDistributionByCycleMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(), 9);
        
        for (final BaseDistributionByCycleMetrics metrics : output.getMetrics()) {
            switch(metrics.CYCLE) {
                case 1:
                case 4:
                case 8:
                    Assert.assertEquals(metrics.PCT_A, 0);
                    Assert.assertEquals(metrics.PCT_C, 100);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 2:
                case 3:
                    Assert.assertEquals(metrics.PCT_A, 100);
                    Assert.assertEquals(metrics.PCT_C, 0);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 5:
                    Assert.assertEquals(metrics.PCT_A, 66.666667);
                    Assert.assertEquals(metrics.PCT_C, 33.333333);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 6:
                    Assert.assertEquals(metrics.PCT_A, 33.333333);
                    Assert.assertEquals(metrics.PCT_C, 66.666667);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 7:
                    Assert.assertEquals(metrics.PCT_A, 33.333333);
                    Assert.assertEquals(metrics.PCT_C, 0);
                    Assert.assertEquals(metrics.PCT_G, 66.666667);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 9:
                    Assert.assertEquals(metrics.PCT_A, 0);
                    Assert.assertEquals(metrics.PCT_C, 0);
                    Assert.assertEquals(metrics.PCT_G, 100);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
            }
        }

        Assert.assertTrue(pdf.length() > 0, "PDF file was not created");
    }

    @Test
    public void testNoChart() throws IOException {
        final File input = new File(TEST_DATA_DIR, "input.sam");
        final File outfile   = File.createTempFile("test", ".basedistributionbycycle");
        final File pdf   = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
    
        final String[] args = new String[] {
            "INPUT="  + input.getAbsolutePath(),
            "OUTPUT=" + outfile.getAbsolutePath(),
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<BaseDistributionByCycleMetrics, Comparable<?>> output = new MetricsFile<BaseDistributionByCycleMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        Assert.assertEquals(output.getMetrics().size(), 9);
        
        for (final BaseDistributionByCycleMetrics metrics : output.getMetrics()) {
            switch(metrics.CYCLE) {
                case 1:
                case 4:
                case 8:
                    Assert.assertEquals(metrics.PCT_A, 0);
                    Assert.assertEquals(metrics.PCT_C, 100);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 2:
                case 3:
                    Assert.assertEquals(metrics.PCT_A, 100);
                    Assert.assertEquals(metrics.PCT_C, 0);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 5:
                    Assert.assertEquals(metrics.PCT_A, 66.666667);
                    Assert.assertEquals(metrics.PCT_C, 33.333333);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 6:
                    Assert.assertEquals(metrics.PCT_A, 33.333333);
                    Assert.assertEquals(metrics.PCT_C, 66.666667);
                    Assert.assertEquals(metrics.PCT_G, 0);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 7:
                    Assert.assertEquals(metrics.PCT_A, 33.333333);
                    Assert.assertEquals(metrics.PCT_C, 0);
                    Assert.assertEquals(metrics.PCT_G, 66.666667);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
                case 9:
                    Assert.assertEquals(metrics.PCT_A, 0);
                    Assert.assertEquals(metrics.PCT_C, 0);
                    Assert.assertEquals(metrics.PCT_G, 100);
                    Assert.assertEquals(metrics.PCT_T, 0);
                    Assert.assertEquals(metrics.PCT_N, 0);
                    break;
            }
        }

        Assert.assertEquals(pdf.length(), 0);
    }

    @Test
    public void testFailureGatkLiteDocker() throws IOException {
        final PrintStream stderr = System.err;
        final String gatkLiteDockerProperty = System.getProperty("IN_GATKLITE_DOCKER");

        try {
            final ByteArrayOutputStream stderrCapture = new ByteArrayOutputStream();
            System.setErr(new PrintStream(stderrCapture));

            System.setProperty("IN_GATKLITE_DOCKER", "true");

            final File input = new File(TEST_DATA_DIR, "input.sam");
            final File outfile   = File.createTempFile("test", ".basedistributionbycycle");
            final File pdf   = File.createTempFile("test", ".pdf");
            outfile.deleteOnExit();
            pdf.deleteOnExit();
            final String[] args = new String[] {
                    "INPUT="  + input.getAbsolutePath(),
                    "OUTPUT=" + outfile.getAbsolutePath(),
                    "CHART=" + pdf.getAbsolutePath()
            };
            Assert.assertEquals(runPicardCommandLine(args), 1);

            Assert.assertEquals(pdf.length(), 0);

            Assert.assertTrue(stderrCapture.toString().contains("The histogram file cannot be written because it requires R, which is not available in the GATK Lite Docker image."));      
        }
        finally {
            System.setErr(stderr);
            if(gatkLiteDockerProperty != null) {
                System.setProperty("IN_GATKLITE_DOCKER", gatkLiteDockerProperty);
            }
            else{
                System.clearProperty("IN_GATKLITE_DOCKER");
            } 
        }
    }
}
