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

        final MetricsFile<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_BASES, 606);
            Assert.assertEquals(metrics.TOTAL_USABLE_BASES, 238);
            Assert.assertEquals(metrics.PCT_EXC_OVERLAP, 0.085809);  // 52 of 606 bases
            Assert.assertEquals(metrics.PCT_EXC_BASEQ, 0.188119);    // 114 of 606 bases
            Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.333333);    // 202 of 606 bases
            Assert.assertEquals(metrics.TOTAL_READ_PAIRS, 3);
            Assert.assertEquals(metrics.TOTAL_DUPE_PAIRS, 1);
            Assert.assertEquals(metrics.TOTAL_ORIENTED_PAIRS, 2);
            Assert.assertEquals(metrics.MEAN_INSERT_SIZE, 118.0);
        }
    }
}