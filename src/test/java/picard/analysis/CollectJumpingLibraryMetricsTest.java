package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CollectJumpingLibraryMetricsTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");
    private static final File SAM_FILE = new File(TEST_DATA_DIR, "forMetrics.sam");

    @Test
    public void testCollectJumpingLibraryMetrics() throws IOException {
        final File outfile = File.createTempFile("CollectJumpingLibraryMetricsTest", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + SAM_FILE.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };
        CollectJumpingLibraryMetrics collectJumpingLibraryMetrics = new CollectJumpingLibraryMetrics();
        Assert.assertEquals(collectJumpingLibraryMetrics.instanceMain(args), 0,
                "Can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        final MetricsFile<JumpingLibraryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final JumpingLibraryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.JUMP_PAIRS, 4);
            Assert.assertEquals(metrics.JUMP_DUPLICATE_PAIRS, 1);
            Assert.assertEquals(metrics.JUMP_DUPLICATE_PCT, 0.25);
            Assert.assertEquals(metrics.JUMP_LIBRARY_SIZE, 6);
            Assert.assertEquals(metrics.JUMP_MEAN_INSERT_SIZE, 176.0);
            Assert.assertEquals(metrics.JUMP_STDEV_INSERT_SIZE, 50.0);
            Assert.assertEquals(metrics.NONJUMP_PAIRS, 1);
            Assert.assertEquals(metrics.NONJUMP_DUPLICATE_PAIRS, 0);
            Assert.assertEquals(metrics.NONJUMP_DUPLICATE_PCT, 0.0);
            Assert.assertEquals(metrics.NONJUMP_LIBRARY_SIZE, 0);
            Assert.assertEquals(metrics.NONJUMP_MEAN_INSERT_SIZE, 96.0);
            Assert.assertEquals(metrics.NONJUMP_STDEV_INSERT_SIZE, Double.NaN);
            Assert.assertEquals(metrics.CHIMERIC_PAIRS, 0);
            Assert.assertEquals(metrics.FRAGMENTS, 1);
            Assert.assertEquals(metrics.PCT_JUMPS, 0.8);
            Assert.assertEquals(metrics.PCT_NONJUMPS, 0.2);
            Assert.assertEquals(metrics.PCT_CHIMERAS, 0.0);
        }
    }
}
