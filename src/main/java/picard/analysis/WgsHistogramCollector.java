package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;

public class WgsHistogramCollector {
    private long[] highQualityDepthHistogramArray;
    protected final long[] unfilteredDepthHistogramArray;
    protected final long[] unfilteredBaseQHistogramArray;

    public WgsHistogramCollector(final int coverageCap) {
        unfilteredDepthHistogramArray = new long[coverageCap + 1];
        highQualityDepthHistogramArray = new long[coverageCap + 1];
        unfilteredBaseQHistogramArray = new long[Byte.MAX_VALUE];
    }

    protected void addBaseQHistogram(final MetricsFile<WgsMetrics, Integer> file) {
        file.addHistogram(getUnfilteredBaseQHistogram());
    }

    protected Histogram<Integer> getHighQualityDepthHistogram() {
        return getHistogram(highQualityDepthHistogramArray, "coverage", "high_quality_coverage_count");
    }

    protected Histogram<Integer> getUnfilteredDepthHistogram() {
        return getHistogram(unfilteredDepthHistogramArray, "coverage", "unfiltered_coverage_count");
    }

    protected Histogram<Integer> getUnfilteredBaseQHistogram() {
        return getHistogram(unfilteredBaseQHistogramArray, "baseq", "unfiltered_baseq_count");
    }

    protected Histogram<Integer> getHistogram(final long[] array, final String binLabel, final String valueLabel) {
        final Histogram<Integer> histogram = new Histogram<>(binLabel, valueLabel);
        for (int i = 0; i < array.length; ++i) {
            histogram.increment(i, array[i]);
        }
        return histogram;
    }
}
