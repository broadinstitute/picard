package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import static org.testng.Assert.*;

public class CollectUmiPrevalenceMetricsTest extends CommandLineProgramTest {

    @Override
    public String getCommandLineProgramName() { return CollectUmiPrevalenceMetrics.class.getSimpleName();}


    @Test
    public void integrationTest() throws IOException {
        final File TEST_DIR = new File("testdata/picard/independent_replicates/");

        final File samFile = new File(TEST_DIR, "twopairsWithManyUmis.sam");
        final File metricsFile = File.createTempFile("test", ".umi_hist");
        metricsFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "MINIMUM_MQ=20"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<?, Integer> output = new MetricsFile<>();
        output.read(new FileReader(metricsFile));

        final Histogram<Integer> hist = output.getHistogram();

        assertNull(hist.get(1));
        assertEquals(hist.get(2).getValue(),4);
        assertEquals(hist.get(3).getValue(),1);
        assertEquals(hist.get(4).getValue(),1);
        assertNull(hist.get(5));
        assertNull(hist.get(6));
        assertEquals(hist.getCount(),6); // sets
        assertEquals(hist.getSum(),15); // unique barcodes (in sets)
    }
}