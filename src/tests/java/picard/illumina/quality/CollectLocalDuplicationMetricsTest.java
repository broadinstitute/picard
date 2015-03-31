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
    private static final File TEST_DATA_DIR = new File("testdata/picard/illumina/CollectLocalDuplicationMetrics");

    private File rootTestDir;

    @BeforeTest
    private void setUp() throws Exception {
        rootTestDir = File.createTempFile("cibm.", ".tmp");
        Assert.assertTrue(rootTestDir.delete());
        Assert.assertTrue(rootTestDir.mkdir());
        for (final File source : TEST_DATA_DIR.listFiles()) {
            if (source.isDirectory() && !source.isHidden()) {
                IOUtil.copyDirectoryTree(source, new File(rootTestDir.getPath(),source.getName()));
            }
        }
    }

    @AfterTest
    private void tearDown() {
        IOUtil.deleteDirectoryTree(rootTestDir);
    }

    @Test
    public void testIndexedRunLane1() throws Exception {
        final MetricsFile<LocalDuplicationSummaryMetrics, Integer> metricsFile = runIt(1, "25T8B25T","25T8B25T/Data/Intensities/BaseCalls", true);

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
    public void testNonIndexedRunLane1() throws Exception {
        final MetricsFile<LocalDuplicationSummaryMetrics, Integer> metricsFile = runIt(1, "125T125T","125T125T/Data/Intensities/BaseCalls",false);
        final LocalDuplicationSummaryMetrics laneMetric = metricsFile.getMetrics().get(0);

        Assert.assertEquals(laneMetric.LANE, "1");



        Assert.assertEquals(metricsFile.getMetrics().size(),1);
    }

    private MetricsFile<LocalDuplicationSummaryMetrics, Integer> runIt(final int lane, final String readStructure, final String basecallsDirName, final boolean isIndexed) throws Exception {
        final File metricsFile = File.createTempFile("ldsm.", ".metrics");
        metricsFile.deleteOnExit();

        File basecallsDir = new File(rootTestDir.getPath(),basecallsDirName);

        ArrayList<String> argsList = new ArrayList<String>();
        argsList.add("BASECALLS_DIR=" + basecallsDir.getPath());
        argsList.add("LANE=" + lane);
        argsList.add("OUTPUT=" + metricsFile.getPath());

        if (readStructure != null) argsList.add("READ_STRUCTURE=" + readStructure);
        if (isIndexed) argsList.add("INPUT=" + new File(basecallsDir.getPath(),"barcodeData." + lane).getPath());

        final String[] args = new String[argsList.size()];
        argsList.toArray(args);

        Assert.assertEquals(new CollectLocalDuplicationMetrics().instanceMain(args),0);

        final MetricsFile<LocalDuplicationSummaryMetrics,Integer> retval =
                new MetricsFile<LocalDuplicationSummaryMetrics,Integer>();
        retval.read(new FileReader(metricsFile));
        return retval;
    }
}