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
    private static final File TEST_DATA_DIR = new File("testdata/picard/illumina/CollectIlluminaBasecallingMetrics");

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
    public void testSummaryMetrics() throws Exception {
        final MetricsFile<LocalDuplicationSummaryMetrics, Integer> metricsFile = runSummaryMetrics(1, 24, "25T8B25T/Data/Intensities/BaseCalls");

        /** Need to fill in the appropriate assertions here */
        final LocalDuplicationSummaryMetrics metric1 = metricsFile.getMetrics().get(0);
        Assert.assertEquals(metric1.LANE, "1");

        final LocalDuplicationSummaryMetrics metric2 = metricsFile.getMetrics().get(1);
        Assert.assertEquals(metric2.LANE, "1");

        final LocalDuplicationSummaryMetrics metric3 = metricsFile.getMetrics().get(2);
        Assert.assertEquals(metric3.LANE, "1");

        final LocalDuplicationSummaryMetrics metric4 = metricsFile.getMetrics().get(3);
        Assert.assertEquals(metric4.LANE, "1");

        final LocalDuplicationSummaryMetrics metric5 = metricsFile.getMetrics().get(4);
        Assert.assertEquals(metric5.LANE, "1");

        final LocalDuplicationSummaryMetrics laneMetric = metricsFile.getMetrics().get(34);
        Assert.assertEquals(laneMetric.LANE, "1");

    }

    @Test
    public void testDetailMetrics() throws Exception {
        final MetricsFile<LocalDuplicationDetailMetrics, Integer> metricsFile = runDetailMetrics(1, 24, "125T125T/Data/Intensities/BaseCalls");
        final LocalDuplicationDetailMetrics laneMetric = metricsFile.getMetrics().get(0);

    }

    private MetricsFile<LocalDuplicationSummaryMetrics, Integer> runSummaryMetrics(final int lane, final int numBases, final String basecallsDirName) throws Exception {
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

    private MetricsFile<LocalDuplicationDetailMetrics, Integer> runDetailMetrics(final int lane, final int numBases, final String basecallsDirName) throws Exception {
        final File metricsFile = File.createTempFile("lddm", null);
        final String basename = metricsFile.getPath();
        final String extension = "." + CollectLocalDuplicationMetrics.DETAILED_METRICS_EXTENSION;
        final File outputMetricsFile = new File(basename + extension);
        metricsFile.deleteOnExit();
        outputMetricsFile.deleteOnExit();

        File basecallsDir = new File(rootTestDir.getPath(),basecallsDirName);

        ArrayList<String> argsList = new ArrayList<String>();
        argsList.add("BASECALLS_DIR=" + basecallsDir.getPath());
        argsList.add("LANE=" + lane);
        argsList.add("OUTPUT=" + basename);
        argsList.add("NUM_BASES=" + numBases);
        argsList.add("PROB_EXPLICIT_OUTPUT=1.0");
        argsList.add("TILE_INDEX=0");
        argsList.add("NUM_TILES=" + CollectLocalDuplicationMetrics.TILES_PER_LANE);

        final String[] args = new String[argsList.size()];
        argsList.toArray(args);

        Assert.assertEquals(new CollectLocalDuplicationMetrics().instanceMain(args),0);

        final MetricsFile<LocalDuplicationDetailMetrics,Integer> retval =
                new MetricsFile<LocalDuplicationDetailMetrics,Integer>();
        retval.read(new FileReader(outputMetricsFile));
        return retval;
    }
}