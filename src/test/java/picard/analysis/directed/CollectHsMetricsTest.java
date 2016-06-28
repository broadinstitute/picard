package picard.analysis.directed;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class CollectHsMetricsTest extends CommandLineProgramTest {
    private final static File TEST_DIR = new File("testdata/picard/analysis/directed/CollectHsMetrics");

    @Override
    public String getCommandLineProgramName() {
        return CollectHsMetrics.class.getSimpleName();
    }

    @DataProvider(name = "collectHsMetricsDataProvider")
    public Object[][] targetedIntervalDataProvider() {
        final String referenceFile = TEST_DIR + "/chrM.fasta";
        final String intervals = TEST_DIR + "/chrM.interval_list";

        return new Object[][] {
                // test that all bases (read 2) with base quality 1 are filtered out
                {TEST_DIR + "/lowbaseq.sam",    referenceFile, intervals, "NONE", 1, 1, true,  2, 202, 0.5, 0.0, 0.505, 0.0,   1000},
                // test that read 2 (with mapping quality 1) is filtered out with minimum mapping quality 2
                {TEST_DIR + "/lowmapq.sam",     referenceFile, intervals, "NONE", 2, 0, true,  2, 202, 0,   0.0, 0.505, 0.0,   1000},
                // test that we clip overlapping bases
                {TEST_DIR + "/overlapping.sam", referenceFile, intervals, "NONE", 0, 0, true,  2, 202, 0,   0.5, 0.505, 0.505, 1000},
                // test that we do not clip overlapping bases
                {TEST_DIR + "/overlapping.sam", referenceFile, intervals, "NONE", 0, 0, false, 2, 202, 0,   0.0, 0.505, 0.505, 1000}
        };
    }

    @Test(dataProvider = "collectHsMetricsDataProvider")
    public void runCollectTargetedMetricsTest(final String input,
                                              final String referenceFile,
                                              final String targetIntervals,
                                              final String metricsFile,
                                              final int minimumMappingQuality,
                                              final int minimumBaseQuality,
                                              final boolean clipOverlappingReads,
                                              final int totalReads,
                                              final int pfUqBasesAligned,
                                              final double pctExcBaseq,
                                              final double pctExcOverlap,
                                              final double pctTargetBases1x,
                                              final double pctTargetBases2x,
                                              final int sampleSize) throws IOException {

        final File outfile = File.createTempFile("CollectHsMetrics", ".hs_metrics", TEST_DIR);
        outfile.deleteOnExit();

        final String[] args = new String[] {
                "TARGET_INTERVALS=" + targetIntervals,
                "BAIT_INTERVALS=" + targetIntervals,
                "INPUT=" + input,
                "OUTPUT=" + outfile,
                "MINIMUM_MAPPING_QUALITY=" + minimumMappingQuality,
                "MINIMUM_BASE_QUALITY=" + minimumBaseQuality,
                "CLIP_OVERLAPPING_READS=" + clipOverlappingReads,
                "SAMPLE_SIZE=" + sampleSize
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<HsMetrics, Comparable<?>> output = new MetricsFile<HsMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final HsMetrics metrics : output.getMetrics()) {
            // overlap
            Assert.assertEquals(metrics.TOTAL_READS, totalReads);
            Assert.assertEquals(metrics.PF_UQ_BASES_ALIGNED, pfUqBasesAligned);
            Assert.assertEquals(metrics.PCT_EXC_BASEQ, pctExcBaseq);
            Assert.assertEquals(metrics.PCT_EXC_OVERLAP, pctExcOverlap);
            Assert.assertEquals(metrics.PCT_TARGET_BASES_1X, pctTargetBases1x);
            Assert.assertEquals(metrics.PCT_TARGET_BASES_2X, pctTargetBases2x);
        }
    }
}
