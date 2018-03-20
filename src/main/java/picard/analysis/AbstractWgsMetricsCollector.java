/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SequenceUtil;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;

/**
 * Class for collecting data on reference coverage, base qualities and excluded bases from one AbstractLocusInfo object for
 * CollectWgsMetrics.
 * <p>
 * The shared code for forming result for CollectWgsMetrics is abstracted into this class.
 * Classes that extend this collector implement their logic in addInfo() method.
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */

public abstract class AbstractWgsMetricsCollector<T extends AbstractRecordAndOffset> {

    /**
     * The source CollectWgsMetrics object
     */
    final CollectWgsMetrics collectWgsMetrics;
    /** Count of sites with a given depth of coverage. Includes all but quality 2 bases.
     * We draw depths from this histogram when we calculate the theoretical het sensitivity.
     */
    protected final long[] unfilteredDepthHistogramArray;
    /** Count of bases observed with a given base quality. Includes all but quality 2 bases.
     * We draw bases from this histogram when we calculate the theoretical het sensitivity.
     */
    protected final long[] unfilteredBaseQHistogramArray;
    /**
     * Count of sites with a given depth of coverage.
     * Excludes bases with quality below MINIMUM_BASE_QUALITY (default 20).
     */
    protected final long[] highQualityDepthHistogramArray;
    /**
     * Number of aligned bases that were filtered out because they were of low base quality (default is < 20).
     */
    long basesExcludedByBaseq = 0;
    /**
     * Number of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
     */
    long basesExcludedByOverlap = 0;
    /**
     * Number of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x).
     */
    long basesExcludedByCapping = 0;
    /**
     * Positions with coverage exceeding this value are treated as if they had coverage at this value
     */
    protected final int coverageCap;

    protected final IntervalList intervals;
    /**
     * This value indicates that processing will stop after specified int the metric amount of genomic bases.
     */
    private final boolean usingStopAfter;
    /**
     * The number of processed genomic bases
     */
    protected long counter = 0;

    /**
     * Creates a collector and initializes the inner data structures
     *
     * @param collectWgsMetrics CollectWgsMetrics, that creates this collector
     * @param coverageCap       coverage cap
     */
    AbstractWgsMetricsCollector(CollectWgsMetrics collectWgsMetrics, final int coverageCap, final IntervalList intervals) {
        if (coverageCap <= 0) {
            throw new IllegalArgumentException("Coverage cap must be positive.");
        }
        this.collectWgsMetrics = collectWgsMetrics;
        unfilteredDepthHistogramArray = new long[coverageCap + 1];
        highQualityDepthHistogramArray = new long[coverageCap + 1];
        unfilteredBaseQHistogramArray = new long[Byte.MAX_VALUE];
        this.coverageCap    = coverageCap;
        this.intervals      = intervals;
        this.usingStopAfter = collectWgsMetrics.STOP_AFTER > 0;
    }

    /**
     * Accumulates the data from AbstractLocusInfo in inner structures
     * @param info {@link htsjdk.samtools.util.AbstractLocusInfo} with aligned to reference position reads
     * @param ref  {@link htsjdk.samtools.reference.ReferenceSequence}
     * @param referenceBaseN true if current the value of reference base represents a no call
     */
    public abstract void addInfo(final AbstractLocusInfo<T> info, final ReferenceSequence ref, boolean referenceBaseN);

    /**
     * Adds collected metrics and depth histogram to file
     * @param file MetricsFile for result of collector's work
     * @param dupeFilter         counting filter for duplicate reads
     * @param mapqFilter         counting filter for mapping quality
     * @param pairFilter         counting filter for reads without a mapped mate pair
     */
    public void addToMetricsFile(final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> file,
            final boolean includeBQHistogram,
            final CountingFilter dupeFilter,
            final CountingFilter mapqFilter,
            final CountingPairedFilter pairFilter) {
        final CollectWgsMetrics.WgsMetrics
                metrics = getMetrics(dupeFilter, mapqFilter, pairFilter);

        // add them to the file
        file.addMetric(metrics);
        file.addHistogram(getHighQualityDepthHistogram());
        if (includeBQHistogram) addBaseQHistogram(file);
    }

    protected void addBaseQHistogram(final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> file) {
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

    /**
     * Creates CollectWgsMetrics.WgsMetrics - the object holding the result of CollectWgsMetrics
     *
     * @param dupeFilter     counting filter for duplicate reads
     * @param mapqFilter     counting filter for mapping quality
     * @param pairFilter     counting filter for reads without a mapped mate pair
     * @return CollectWgsMetrics.WgsMetrics with set fields
     */
    protected CollectWgsMetrics.WgsMetrics getMetrics(final CountingFilter dupeFilter,
            final CountingFilter mapqFilter,
            final CountingPairedFilter pairFilter) {
        return collectWgsMetrics.generateWgsMetrics(
                this.intervals,
                getHighQualityDepthHistogram(),
                getUnfilteredDepthHistogram(),
                collectWgsMetrics.getBasesExcludedBy(mapqFilter),
                collectWgsMetrics.getBasesExcludedBy(dupeFilter),
                collectWgsMetrics.getBasesExcludedBy(pairFilter),
                basesExcludedByBaseq,
                basesExcludedByOverlap,
                basesExcludedByCapping,
                coverageCap,
                getUnfilteredBaseQHistogram(),
                collectWgsMetrics.SAMPLE_SIZE);
    }

    /**
     * @return true, of number of processed loci exceeded the threshold, otherwise false
     */
    boolean isTimeToStop(final long processedLoci) {
        return usingStopAfter && processedLoci > collectWgsMetrics.STOP_AFTER - 1;
    }

    /**
     * Sets the counter to the current number of processed loci. Counter, must be updated
     * from outside, since we are skipping a no call reference positions outside of the collector
     *
     * @param counter number of processed loci
     */
    public void setCounter(long counter) {
        this.counter = counter;
    }

    /**
     * Checks if reference base at given position is unknown.
     *
     * @param position to check the base
     * @param ref      reference sequence
     * @return true if reference base at position represents a no call, otherwise false
     */
    boolean isReferenceBaseN(final int position, final ReferenceSequence ref) {
        final byte base = ref.getBases()[position - 1];
        return SequenceUtil.isNoCall(base);
    }
}
