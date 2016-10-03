package picard.illumina;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;
import picard.illumina.parser.ReadStructure;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;import java.lang.Exception;import java.lang.Object;import java.lang.String;

/** @author mccowan */
public class IlluminaLaneMetricsCollectorTest {
    final static File TEST_DIRECTORY = new File("testdata/picard/illumina/IlluminaLaneMetricsCollectorTest");
    final static File TILE_RUN_DIRECTORY = new File(TEST_DIRECTORY, "tileRuns");

    private static File buildOutputFile(final File directory, final String prefix, final String extension) {
        return new File(directory, String.format("%s.%s", prefix, extension));
    }

    @Test(dataProvider = "testLaneMetrics")
    public void testWriteLaneMetrics(final String testRun) throws Exception {
        final CollectIlluminaLaneMetrics clp = new CollectIlluminaLaneMetrics();
        clp.OUTPUT_DIRECTORY = IOUtil.createTempDir("illuminaLaneMetricsCollectorTest", null);
        clp.RUN_DIRECTORY = new File(TEST_DIRECTORY, testRun);
        clp.OUTPUT_PREFIX = "test";
        clp.READ_STRUCTURE = new ReadStructure("101T8B101T");
        clp.doWork();

        final File laneMetricsFile = buildOutputFile(clp.OUTPUT_DIRECTORY, clp.OUTPUT_PREFIX, IlluminaLaneMetrics.getExtension());
        final File canonicalOutputFile = buildOutputFile(TEST_DIRECTORY, testRun, IlluminaLaneMetrics.getExtension());

        IOUtil.assertFilesEqual(canonicalOutputFile, laneMetricsFile);

        IOUtil.deleteDirectoryTree(clp.OUTPUT_DIRECTORY);
    }

    @DataProvider(name = "testLaneMetrics")
    public Object[][] testLaneMetricsDataProvider() {
        return new Object[][] {
                {"130321_SL-MAK_0035_FC000000000-A306B"},
                {"130318_SL-HBB_0226_BFCC1WYMACXX"},
                {"130401_SL-HAC_0022_BH07PBADXX"}
        };
    }

    @Test(dataProvider = "testCollectIlluminaLaneMetrics")
    public void testCollectIlluminaLaneMetrics(final String testRun, final ReadStructure readStructure) throws Exception {
        final File runDirectory = new File(TILE_RUN_DIRECTORY, testRun);
        final CollectIlluminaLaneMetrics clp = new CollectIlluminaLaneMetrics();
        clp.OUTPUT_DIRECTORY = IOUtil.createTempDir("illuminaLaneMetricsCollectorTest", null);
        clp.RUN_DIRECTORY = runDirectory;
        clp.OUTPUT_PREFIX = "test";
        clp.READ_STRUCTURE = readStructure;
        clp.doWork();

        final File phasingMetricsPhile = buildOutputFile(clp.OUTPUT_DIRECTORY, clp.OUTPUT_PREFIX, IlluminaPhasingMetrics.getExtension());
        final File canonicalPhasingPhile = buildOutputFile(runDirectory, testRun, IlluminaPhasingMetrics.getExtension());
        IOUtil.assertFilesEqual(canonicalPhasingPhile, phasingMetricsPhile);

        final File laneMetricsFile = buildOutputFile(clp.OUTPUT_DIRECTORY, clp.OUTPUT_PREFIX, IlluminaLaneMetrics.getExtension());
        final File canonicalLaneFile = buildOutputFile(runDirectory, testRun, IlluminaLaneMetrics.getExtension());
        IOUtil.assertFilesEqual(canonicalLaneFile, laneMetricsFile);
        IOUtil.deleteDirectoryTree(clp.OUTPUT_DIRECTORY);
    }

    @DataProvider(name = "testCollectIlluminaLaneMetrics")
    public Object[][] testCollectIlluminaLaneMetricsDataProvider() {
        return new Object[][] {
                {"A7LE0", new ReadStructure("25T8B8B25T")},
                {"C2MFAACXX", new ReadStructure("95T101T")},
                {"H7BATADXX", new ReadStructure("76T8B76T")},
                {"H7H7RADXX", new ReadStructure("101T8B8B101T")},
                {"A67HY", new ReadStructure("8B8B")}
        };
    }
}
