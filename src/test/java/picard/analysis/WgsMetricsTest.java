package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

/**
 * Tests for WgsMetrics.
 */
public class WgsMetricsTest {

    private Histogram<Integer> emptyDepthHistogram() { return new Histogram<>(); }
    private Histogram<Integer> singleDepthHistogram(final int depth, final int count) {
        final Histogram<Integer> histogram = new Histogram<>();
        histogram.increment(depth, count);
        return histogram;
    }
    private Histogram<Integer> twoSiteDepthHistogram(final int depth1, final int count1, final int depth2, final int count2) {
        final Histogram<Integer> histogram = new Histogram<>();
        if (0 < depth1) histogram.increment(depth1, count1);
        if (0 < depth2) histogram.increment(depth2, count2);
        return histogram;
    }

    private IntervalList buildIntervalList(final int start, final int end) {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new SAMSequenceRecord("CONTIG", 100000000));
        final IntervalList intervals = new IntervalList(header);
        if (0 < start) intervals.add(new Interval("CONTIG", start, end));
        return intervals;
    }

    private CollectWgsMetrics.WgsMetrics emptyMetrics() {
        return new CollectWgsMetrics.WgsMetrics(
                buildIntervalList(-1, -1),
                emptyDepthHistogram(), emptyDepthHistogram(),
                0, 0, 0, 0, 0, 0, 0, 1000000,
                null, -1
        );
    }

    private CollectWgsMetrics.WgsMetrics singleDepthMetrics(final int depth, final int countScale, final int start) {
        final int count = 100000 * countScale;
        final int totalExcluded = (10 + 20 + 30 + 40 + 50 + 60) * countScale;
        return new CollectWgsMetrics.WgsMetrics(
                buildIntervalList(start, start),
                singleDepthHistogram(depth, count),
                singleDepthHistogram(depth, count),
                10d * countScale / count, 20d * countScale / count, 30d * countScale / count,
                40d * countScale / count, 50d * countScale / count, 60d * countScale / count,
                totalExcluded / (double) (count + totalExcluded),
                1000000,
                null, -1
        );
    }

    private CollectWgsMetrics.WgsMetrics twoSiteDepthMetrics(final int depth1, final int countScale1,
                                                             final int depth2, final int countScale2,
                                                             final int start) {
        final int count1 = 100000 * countScale1;
        final int count2 = 100000 * countScale2;
        final int count  = count1 + count2;
        final int countScale = countScale1 + countScale2;
        final int totalExcluded = (10 + 20 + 30 + 40 + 50 + 60) * countScale;
        return new CollectWgsMetrics.WgsMetrics(
                buildIntervalList(start, start+1),
                twoSiteDepthHistogram(depth1, count1, depth2, count2),
                twoSiteDepthHistogram(depth1, count1, depth2, count2),
                10d * countScale / count, 20d * countScale / count, 30d * countScale / count,
                40d * countScale / count, 50d * countScale / count, 60d * countScale / count,
                totalExcluded / (double) (count + totalExcluded),
                100000,
                null, -1
        );
    }

    @Test(dataProvider = "testWgsMetricsMergeDataProvider")
    public void testWgsMetricsMerge(final CollectWgsMetrics.WgsMetrics left,
                                    final CollectWgsMetrics.WgsMetrics right,
                                    final CollectWgsMetrics.WgsMetrics expected) {
        left.merge(right);
        left.calculateDerivedFields();
        Assert.assertTrue(left.equals(expected));
    }

    @DataProvider(name = "testWgsMetricsMergeDataProvider")
    public Object[][] testWgsMetricsMergeDataProvider() {
        return new Object[][] {
                {emptyMetrics(), emptyMetrics(), emptyMetrics()},
                {emptyMetrics(), singleDepthMetrics(1, 1, 1), singleDepthMetrics(1, 1, 1)},
                {singleDepthMetrics(1, 1, 1), emptyMetrics(), singleDepthMetrics(1, 1, 1)},
                {singleDepthMetrics(1, 1, 1), singleDepthMetrics(1, 1, 2), twoSiteDepthMetrics(1, 2, 0, 0, 1)},
                {singleDepthMetrics(1, 1, 1), singleDepthMetrics(1, 2, 2), twoSiteDepthMetrics(1, 3, 0, 0, 1)},
                {singleDepthMetrics(1, 1, 1), singleDepthMetrics(1, 1, 2), twoSiteDepthMetrics(1, 2, 0, 0, 1)},
                {singleDepthMetrics(1, 1, 1), singleDepthMetrics(2, 1, 2), twoSiteDepthMetrics(1, 1, 2, 1, 1)},
                {twoSiteDepthMetrics(1, 1, 2, 1, 1), twoSiteDepthMetrics(1, 1, 2, 1, 3), twoSiteDepthMetrics(1, 2, 2, 2, 1)}
        };
    }

    @Test(expectedExceptions = {PicardException.class})
    public void testMergeOverlappingIntervals() {
        singleDepthMetrics(1, 1, 1).merge(singleDepthMetrics(1, 1, 1));
    }
}
