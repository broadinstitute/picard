package picard.illumina.quality;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import picard.illumina.quality.CollectLocalDuplicationMetrics;
import picard.illumina.quality.LocalDuplicationSummaryMetrics;
import picard.illumina.quality.LocalDuplicationDetailMetrics;

import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

public class CollectLocalDuplicationMetricsTest {
    /** Use the data already set up for testingCollectIlluminBasecallingMetrics */
    private static final File TEST_DATA_DIR = new File("testdata/picard/illumina");

    private File rootTestDir;

    @BeforeTest
    private void setUp() throws Exception {
        rootTestDir = File.createTempFile("cibm.", ".tmp");
        Assert.assertTrue(rootTestDir.delete());
        Assert.assertTrue(rootTestDir.mkdir());
        /** Test several directories.  We might want to test a single directory in testdata/picard/illumina
         * rather than steal testdata from a different tool.
         */
        for (final File source : TEST_DATA_DIR.listFiles()) {
            if (source.isDirectory() && !source.isHidden()) {
                IOUtil.copyDirectoryTree(source, new File(rootTestDir.getPath(),source.getName()));
            }
        }
    }

    @AfterTest
    private void tearDown() { IOUtil.deleteDirectoryTree(rootTestDir); }

    @Test
    public void testLocalDuplicationSummaryMetrics() throws Exception {
        final MetricsFile<LocalDuplicationSummaryMetrics, Integer> metricsFile0 = runSummaryMetrics(1, 24, "125T125T/Data/Intensities/BaseCalls",0);
        final MetricsFile<LocalDuplicationSummaryMetrics, Integer> metricsFile1000 = runSummaryMetrics(1, 24, "125T125T/Data/Intensities/BaseCalls",1000);
        final MetricsFile<LocalDuplicationSummaryMetrics, Integer> metricsFileWholeTile = runSummaryMetrics(1, 24, "125T125T/Data/Intensities/BaseCalls",Double.POSITIVE_INFINITY);

        //getMetrics().get(0) returns the first ("All") metric
        final double DELTA = 0.0001;
        final LocalDuplicationSummaryMetrics metric0 = metricsFile0.getMetrics().get(0);
        Assert.assertEquals(metric0.LANE, 1);
        Assert.assertEquals(metric0.TILE, "All");
        Assert.assertEquals(metric0.READS, 1863);
        Assert.assertEquals(metric0.LOCAL_DUPLICATES, 0);
        Assert.assertEquals(metric0.PCT_LOCAL_DUPLICATES, ((double) metric0.LOCAL_DUPLICATES) / metric0.READS, DELTA);

        final LocalDuplicationSummaryMetrics metric1000 = metricsFile1000.getMetrics().get(0);
        Assert.assertEquals(metric1000.LANE, 1);
        Assert.assertEquals(metric1000.TILE, "All");
        Assert.assertEquals(metric1000.READS, 1863);
        Assert.assertEquals(metric1000.LOCAL_DUPLICATES, 45);
        Assert.assertEquals(metric1000.PCT_LOCAL_DUPLICATES, ((double) metric1000.LOCAL_DUPLICATES) / metric1000.READS, DELTA);

        final LocalDuplicationSummaryMetrics metricWholeTile = metricsFileWholeTile.getMetrics().get(0);
        Assert.assertEquals(metricWholeTile.LANE, 1);
        Assert.assertEquals(metricWholeTile.TILE, "All");
        Assert.assertEquals(metricWholeTile.READS, 1863);
        Assert.assertEquals(metricWholeTile.LOCAL_DUPLICATES, 118);
        Assert.assertEquals(metricWholeTile.PCT_LOCAL_DUPLICATES, ((double) metricWholeTile.LOCAL_DUPLICATES) / metricWholeTile.READS, DELTA);
    }


    private MetricsFile<LocalDuplicationSummaryMetrics, Integer> runSummaryMetrics(final int lane, final int numBases, final String basecallsDirName, final double maxSeparation) throws Exception {
        final File metricsFile = File.createTempFile("ldsm", null);
        final String basename = metricsFile.getPath();
        final String extension = "." + CollectLocalDuplicationMetrics.SUMMARY_METRICS_EXTENSION;
        final File outputMetricsFile = new File(basename + extension);
        metricsFile.deleteOnExit();
        outputMetricsFile.deleteOnExit();

        File basecallsDir = new File(rootTestDir.getPath(),basecallsDirName);

        ArrayList<String> argsList = new ArrayList<String>();
        argsList.add("BASECALLS_DIR=" + basecallsDir.getPath());
        argsList.add("LANE=" + lane);
        argsList.add("OUTPUT=" + basename);
        argsList.add("NUM_BASES=" + numBases);
        argsList.add("MAX_SEPARATION=" + maxSeparation);
        argsList.add("PROB_EXPLICIT_OUTPUT=0");
        argsList.add("TILE_INDEX=0");
        argsList.add("NUM_TILES=" + CollectLocalDuplicationMetrics.TILES_PER_LANE);

        final String[] args = new String[argsList.size()];
        argsList.toArray(args);

        Assert.assertEquals(new CollectLocalDuplicationMetrics().instanceMain(args),0);

        final MetricsFile<LocalDuplicationSummaryMetrics,Integer> retval =
                new MetricsFile<LocalDuplicationSummaryMetrics,Integer>();
        retval.read(new FileReader(outputMetricsFile));
        return retval;
    }
}