package picard.analysis;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import picard.cmdline.CommandLineProgramTest;

public class QualityScoreDistributionTest  extends CommandLineProgramTest {
    
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/QualityScoreDistribution");

    public String getCommandLineProgramName() {
        return QualityScoreDistribution.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "input.sam");
        final File outfile   = File.createTempFile("test", ".quality_score_distribution");
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

        Assert.assertEquals(histogram.size(), 12);

        Assert.assertEquals(histogram.get(32).getValue(),17);
        Assert.assertEquals(histogram.get(33).getValue(),7);
        Assert.assertEquals(histogram.get(34).getValue(),2);
        Assert.assertEquals(histogram.get(35).getValue(),5);
        Assert.assertEquals(histogram.get(37).getValue(),7);
        Assert.assertEquals(histogram.get(38).getValue(),1);
        Assert.assertEquals(histogram.get(64).getValue(),3);
        Assert.assertEquals(histogram.get(68).getValue(),1);
        Assert.assertEquals(histogram.get(73).getValue(),2);
        Assert.assertEquals(histogram.get(75).getValue(),1);
        Assert.assertEquals(histogram.get(78).getValue(),1);
        Assert.assertEquals(histogram.get(82).getValue(),2);

        Assert.assertTrue(pdf.length() > 0);
    }

    @Test
    public void testNoChart() throws IOException {
        final File input = new File(TEST_DATA_DIR, "input.sam");
        final File outfile   = File.createTempFile("test", ".quality_score_distribution");
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

        Assert.assertEquals(histogram.size(), 12);

        Assert.assertEquals(histogram.get(32).getValue(),17);
        Assert.assertEquals(histogram.get(33).getValue(),7);
        Assert.assertEquals(histogram.get(34).getValue(),2);
        Assert.assertEquals(histogram.get(35).getValue(),5);
        Assert.assertEquals(histogram.get(37).getValue(),7);
        Assert.assertEquals(histogram.get(38).getValue(),1);
        Assert.assertEquals(histogram.get(64).getValue(),3);
        Assert.assertEquals(histogram.get(68).getValue(),1);
        Assert.assertEquals(histogram.get(73).getValue(),2);
        Assert.assertEquals(histogram.get(75).getValue(),1);
        Assert.assertEquals(histogram.get(78).getValue(),1);
        Assert.assertEquals(histogram.get(82).getValue(),2);

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
            final File outfile   = File.createTempFile("test", ".quality_score_distribution");
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
