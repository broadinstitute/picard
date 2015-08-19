package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.*;

/**
 * Tests for methods in CollectWgsMetricsFromQuerySorted
 *
 *
 * @author Eric Banks
 */

public class CollectWgsMetricsFromQuerySortedTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    public String getCommandLineProgramName() {
        return CollectWgsMetricsFromQuerySorted.class.getSimpleName();
    }


    @Test
    public void testMetricsFromClippedOverhangs() throws IOException {
        final File input = new File(TEST_DATA_DIR, "namesorted.test.sam");
        final File outfile   = File.createTempFile("metrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetricsFromQuerySorted.MySeqMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetricsFromQuerySorted.MySeqMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetricsFromQuerySorted.MySeqMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_BASES, 404);
            Assert.assertEquals(metrics.TOTAL_USABLE_BASES, 238);
            Assert.assertEquals(metrics.PCT_EXC_OVERLAP, 0.128713);  // 52 of 404 bases
            Assert.assertEquals(metrics.PCT_EXC_BASEQ, 0.282178);    // 114 of 404 bases
        }
    }
}