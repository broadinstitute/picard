package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import picard.cmdline.CommandLineProgramTest;

public class MeanQualityByCycleTest extends CommandLineProgramTest {
    
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/MeanQualityByCycle");

    public String getCommandLineProgramName() {
        return MeanQualityByCycle.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "input.sam");
        final File outfile   = File.createTempFile("test", ".mean_quality_by_cycle");
        final File pdf   = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "CHART=" + pdf.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<?, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final Histogram<Integer> histogram = output.getHistogram();

        Assert.assertEquals(histogram.size(), 9);

        Assert.assertEquals(histogram.get(1).getValue(),45.166667);
        Assert.assertEquals(histogram.get(2).getValue(),42.333333);
        Assert.assertEquals(histogram.get(3).getValue(),40.833333);
        Assert.assertEquals(histogram.get(4).getValue(),38.5);
        Assert.assertEquals(histogram.get(5).getValue(),32.666667);
        Assert.assertEquals(histogram.get(6).getValue(),39.333333);
        Assert.assertEquals(histogram.get(7).getValue(),54);
        Assert.assertEquals(histogram.get(8).getValue(),41.5);
        Assert.assertEquals(histogram.get(9).getValue(),32);

        Assert.assertTrue(pdf.length() > 0);
    }

    @Test
    public void testNoChart() throws IOException {
        final File input = new File(TEST_DATA_DIR, "input.sam");
        final File outfile   = File.createTempFile("test", ".mean_quality_by_cycle");
        final File pdf   = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<?, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final Histogram<Integer> histogram = output.getHistogram();

        Assert.assertEquals(histogram.size(), 9);

        Assert.assertEquals(histogram.get(1).getValue(),45.166667);
        Assert.assertEquals(histogram.get(2).getValue(),42.333333);
        Assert.assertEquals(histogram.get(3).getValue(),40.833333);
        Assert.assertEquals(histogram.get(4).getValue(),38.5);
        Assert.assertEquals(histogram.get(5).getValue(),32.666667);
        Assert.assertEquals(histogram.get(6).getValue(),39.333333);
        Assert.assertEquals(histogram.get(7).getValue(),54);
        Assert.assertEquals(histogram.get(8).getValue(),41.5);
        Assert.assertEquals(histogram.get(9).getValue(),32);

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
            final File outfile   = File.createTempFile("test", ".mean_quality_by_cycle");
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
