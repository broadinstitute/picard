package picard.analysis.directed;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;

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
        final String halfIntervals = TEST_DIR + "/chrM_100bp.interval_list";
        final String twoSmallIntervals = TEST_DIR + "/two-small.interval_list";

        return new Object[][] {
                // two reads, each has 100 bases. bases in one read are medium quality (20), in the other read poor quality (10).
                // test that we exclude half of the bases
                {TEST_DIR + "/lowbaseq.sam",    intervals, 1, 10, true,  2, 200, 0.5, 0.0, 0.50, 0.0,  1, 0, 200, 1000},
                // test that read 2 (with mapping quality 1) is filtered out with minimum mapping quality 2
                {TEST_DIR + "/lowmapq.sam",     intervals, 2, 0, true,  2, 202, 0,   0.0, 0.505, 0.0,   1, 0, 202, 1000},
                // test that we clip overlapping bases
                {TEST_DIR + "/overlapping.sam", intervals, 0, 0, true,  2, 202, 0,   0.5, 0.505, 0, 1, 0, 202, 1000},
                // test that we do not clip overlapping bases
                {TEST_DIR + "/overlapping.sam", intervals, 0, 0, false, 2, 202, 0,   0.0, 0.505, 0.505, 2, 0, 202, 1000},
                // test that we exclude half of the bases (due to poor quality) with an interval that is completely covered
                {TEST_DIR + "/lowbaseq.sam",    halfIntervals, 1, 10, true,  2, 200, 0.5, 0.0, 1.0, 0.0,  1, 1, 200, 1000},
                // test that read 2 (with mapping quality 1) is filtered out with minimum mapping quality 2 with an interval that is completely covered
                {TEST_DIR + "/lowmapq.sam",     halfIntervals, 2, 0, true,  2, 202, 0,   0.0, 1.0, 0.0,   1, 1, 202, 1000},
                // test that we clip overlapping bases with an interval that is completely covered
                {TEST_DIR + "/overlapping.sam", halfIntervals, 0, 0, true,  2, 202, 0,   0.5, 1.0, 0, 1, 1, 202, 1000},
                // test that we do not clip overlapping bases with an interval that is completely covered
                {TEST_DIR + "/overlapping.sam", halfIntervals, 0, 0, false, 2, 202, 0,   0.0, 1.0, 1.0, 2, 2, 202, 1000},
                // A read 10 base pairs long. two intervals: one maps identically to the read, other does not overlap at all
                {TEST_DIR + "/single-short-read.sam", twoSmallIntervals, 20, 20, true, 1, 10, 0.0, 0.0, 0.5, 0.0, 1, 0, 10, 1000 }
        };
    }

    /** Read back the first metrics record in an hs metrics file. */
    private HsMetrics readMetrics(final File f) {
        try {
            final MetricsFile<HsMetrics, Comparable<?>> mFile = new MetricsFile<HsMetrics, Comparable<?>>();
            mFile.read(new FileReader(f));
            return mFile.getMetrics().get(0);
        }
        catch (IOException ioe) {
             throw new RuntimeIOException(ioe);
        }
    }

    /** Writes the contents of a SAMRecordSetBuilder out to a file. */
    private File writeBam(final SAMRecordSetBuilder builder, final File f) {
        try (final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(builder.getHeader(), false, f)) {
            builder.forEach(out::addAlignment);
        }
        return f;
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
                                              final long minTargetCoverage,
                                              final long pfBases,
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

        final HsMetrics metrics = readMetrics(outfile);
        Assert.assertEquals(metrics.TOTAL_READS, totalReads);
        Assert.assertEquals(metrics.PF_UQ_BASES_ALIGNED, pfUqBasesAligned);
        Assert.assertEquals(metrics.PCT_EXC_BASEQ, pctExcBaseq);
        Assert.assertEquals(metrics.PCT_EXC_OVERLAP, pctExcOverlap);
        Assert.assertEquals(metrics.PCT_TARGET_BASES_1X, pctTargetBases1x);
        Assert.assertEquals(metrics.PCT_TARGET_BASES_2X, pctTargetBases2x);
        Assert.assertEquals(metrics.MAX_TARGET_COVERAGE, maxTargetCoverage);
        Assert.assertEquals(metrics.MIN_TARGET_COVERAGE, minTargetCoverage);
        Assert.assertEquals(metrics.PF_BASES, pfBases);
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

    @Test
    public void testHsMetricsHandlesIndelsAppropriately() throws IOException {
        final SAMRecordSetBuilder withDeletions = new SAMRecordSetBuilder(true, SortOrder.coordinate);
        final SAMRecordSetBuilder withInsertions = new SAMRecordSetBuilder(true, SortOrder.coordinate);
        final IntervalList targets = new IntervalList(withDeletions.getHeader());
        final IntervalList baits   = new IntervalList(withDeletions.getHeader());
        targets.add(new Interval("chr1", 1000, 1199, false, "t1"));
        baits.add(new Interval("chr1", 950,  1049, false, "b1"));
        baits.add(new Interval("chr1", 1050, 1149, false, "b2"));
        baits.add(new Interval("chr1", 1150, 1249, false, "b3"));

        // Generate 100 reads that fully cover the the target in each BAM
        for (int i=0; i<100; ++i) {
            withDeletions.addFrag( "d" + i, 0, 1000, false, false, "100M20D80M", null, 30);
            withInsertions.addFrag("i" + i, 0, 1000, false, false, "100M50I100M", null, 30);
        }

        // Write things out to file
        final File dir = IOUtil.createTempDir("hsmetrics.", ".test");
        final File bs = new File(dir, "baits.interval_list").getAbsoluteFile();
        final File ts = new File(dir, "targets.interval_list").getAbsoluteFile();
        baits.write(bs);
        targets.write(ts);
        final File withDelBam = writeBam(withDeletions,  new File(dir, "with_del.bam"));
        final File withInsBam = writeBam(withInsertions, new File(dir, "with_ins.bam"));

        // Now run CollectHsMetrics four times
        final File out = Files.createTempFile("hsmetrics.", ".txt").toFile();
        runPicardCommandLine(Arrays.asList("INCLUDE_INDELS=false", "SAMPLE_SIZE=0", "TI="+ts.getPath(), "BI="+bs.getPath(), "O="+out.getPath(), "I="+withDelBam.getAbsolutePath()));
        final HsMetrics delsWithoutIndelHandling = readMetrics(out);
        runPicardCommandLine(Arrays.asList("INCLUDE_INDELS=true", "SAMPLE_SIZE=0", "TI="+ts.getPath(), "BI="+bs.getPath(), "O="+out.getPath(), "I="+withDelBam.getAbsolutePath()));
        final HsMetrics delsWithIndelHandling = readMetrics(out);
        runPicardCommandLine(Arrays.asList("INCLUDE_INDELS=false", "SAMPLE_SIZE=0", "TI="+ts.getPath(), "BI="+bs.getPath(), "O="+out.getPath(), "I="+withInsBam.getAbsolutePath()));
        final HsMetrics insWithoutIndelHandling = readMetrics(out);
        runPicardCommandLine(Arrays.asList("INCLUDE_INDELS=true", "SAMPLE_SIZE=0", "TI="+ts.getPath(), "BI="+bs.getPath(), "O="+out.getPath(), "I="+withInsBam.getAbsolutePath()));
        final HsMetrics insWithIndelHandling = readMetrics(out);

        IOUtil.deleteDirectoryTree(dir);

        Assert.assertEquals(delsWithoutIndelHandling.MEAN_TARGET_COVERAGE, 90.0);  // 100X over 180/200 bases due to deletion
        Assert.assertEquals(delsWithIndelHandling.MEAN_TARGET_COVERAGE, 100.0);    // 100X with counting the deletion

        Assert.assertEquals(insWithoutIndelHandling.PCT_USABLE_BASES_ON_TARGET, 200/250d); // 50/250 inserted bases are not counted as on target
        Assert.assertEquals(insWithIndelHandling.PCT_USABLE_BASES_ON_TARGET,   1.0d);      // inserted bases are counted as on target
    }
}
