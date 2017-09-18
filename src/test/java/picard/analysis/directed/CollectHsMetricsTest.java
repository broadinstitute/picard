package picard.analysis.directed;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
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
        final String twoSmallIntervals = TEST_DIR + "/two-small.interval_list";

        return new Object[][] {
                // two reads, each has 100 bases. bases in one read are medium quality (20), in the other read poor quality (10).
                // test that we exclude half of the bases
                {TEST_DIR + "/lowbaseq.sam",    intervals, 1, 10, true,  2, 200, 0.5, 0.0, 0.50, 0.0,  1, 1000},
                // test that read 2 (with mapping quality 1) is filtered out with minimum mapping quality 2
                {TEST_DIR + "/lowmapq.sam",     intervals, 2, 0, true,  2, 202, 0,   0.0, 0.505, 0.0,   1, 1000},
                // test that we clip overlapping bases
                {TEST_DIR + "/overlapping.sam", intervals, 0, 0, true,  2, 202, 0,   0.5, 0.505, 0, 1, 1000},
                // test that we do not clip overlapping bases
                {TEST_DIR + "/overlapping.sam", intervals, 0, 0, false, 2, 202, 0,   0.0, 0.505, 0.505, 2, 1000},
                // A read 10 base pairs long. two intervals: one maps identically to the read, other does not overlap at all
                {TEST_DIR + "/single-short-read.sam", twoSmallIntervals, 20, 20, true, 1, 10, 0.0, 0.0, 0.5, 0.0, 1, 1000 }

        };
    }

    @Test(dataProvider = "collectHsMetricsDataProvider")
    public void runCollectHsMetricsTest(final String input,
                                              final String targetIntervals,
                                              final int minimumMappingQuality,
                                              final int minimumBaseQuality,
                                              final boolean clipOverlappingReads,
                                              final int totalReads,
                                              final int pfUqBasesAligned,
                                              final double pctExcBaseq,
                                              final double pctExcOverlap,
                                              final double pctTargetBases1x,
                                              final double pctTargetBases2x,
                                              final long maxTargetCoverage,
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
            Assert.assertEquals(metrics.MAX_TARGET_COVERAGE, maxTargetCoverage);
        }
    }

    @Test
    public void testCoverageHistogram() throws IOException {

        /**
         *  We have a read 10 base pairs long and two intervals: one maps identically to the read, other does not overlap with the read
         *
         *  intervals:    [----------]          [----------]
         *  read:          xxxxxxxxxx
         *
         *  Test that the depth histogram is [10,10,0,...,0]
         */

        final String input = TEST_DIR + "/single-short-read.sam";
        final String targetIntervals = TEST_DIR + "/two-small.interval_list";
        final int minimumMappingQuality = 20;
        final int minimumBaseQuality = 20;
        final boolean clipOverlappingReads = true;
        final int sampleSize = 10;

        final File outfile = File.createTempFile("testCoverageHistogram", ".hs_metrics", TEST_DIR);
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

        final MetricsFile<HsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));
        final Histogram<Integer> coverageHistogram = output.getAllHistograms().get(0);
        Assert.assertEquals(coverageHistogram.get(0).getValue(), 10.0);
        Assert.assertEquals(coverageHistogram.get(1).getValue(), 10.0);
    }
}
