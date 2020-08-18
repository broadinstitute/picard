/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.QualityUtil;
import picard.PicardException;
import picard.util.MathUtil;

/** Metrics for evaluating the performance of whole genome sequencing experiments. */
public class WgsMetrics extends MergeableMetricBase {
    private static final double LOG_ODDS_THRESHOLD = 3.0;


    /** The intervals over which this metric was computed. */
    @MergingIsManual
    protected IntervalList intervals;

    /** Count of sites with a given depth of coverage. Excludes bases with quality below MINIMUM_BASE_QUALITY*/
    @MergingIsManual
    protected final Histogram<Integer> highQualityDepthHistogram;

    /** Count of sites with a given depth of coverage. Includes all but quality 2 bases */
    @MergingIsManual
    protected final Histogram<Integer> unfilteredDepthHistogram;

    /** The count of bases observed with a given base quality. */
    @MergingIsManual
    protected final Histogram<Integer> unfilteredBaseQHistogram;

    /** The maximum depth/coverage to consider. */
    @MergeByAssertEquals
    protected final int coverageCap;

    /** The sample size used for theoretical het sensitivity. */
    @NoMergingKeepsValue
    protected final int theoreticalHetSensitivitySampleSize;

    /**
     * Create an instance of this metric that is not mergeable.
     */
    public WgsMetrics() {
        intervals                           = null;
        highQualityDepthHistogram           = null;
        unfilteredDepthHistogram            = null;
        unfilteredBaseQHistogram            = null;
        theoreticalHetSensitivitySampleSize = -1;
        coverageCap                         = -1;
    }

    /**
     * Create an instance of this metric that is mergeable.
     *
     * @param highQualityDepthHistogram the count of genomic positions observed for each observed depth. Excludes bases with quality below MINIMUM_BASE_QUALITY.
     * @param unfilteredDepthHistogram the depth histogram that includes all but quality 2 bases.
     * @param pctExcludedByAdapter the fraction of aligned bases that were filtered out because they were in reads with 0 mapping quality that looked like adapter sequence.
     * @param pctExcludedByMapq the fraction of aligned bases that were filtered out because they were in reads with low mapping quality.
     * @param pctExcludedByDupes the fraction of aligned bases that were filtered out because they were in reads marked as duplicates.
     * @param pctExcludedByPairing the fraction of bases that were filtered out because they were in reads without a mapped mate pair.
     * @param pctExcludedByBaseq the fraction of aligned bases that were filtered out because they were of low base quality.
     * @param pctExcludedByOverlap the fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
     * @param pctExcludedByCapping the fraction of aligned bases that were filtered out because they would have raised coverage above the capped value.
     * @param pctExcludeTotal the fraction of bases excluded across all filters.
     * @param coverageCap Treat positions with coverage exceeding this value as if they had coverage at this value.
     * @param unfilteredBaseQHistogram the count of bases observed with a given quality. Includes all but quality 2 bases.
     * @param theoreticalHetSensitivitySampleSize the sample size used for theoretical het sensitivity sampling.
     */
    public WgsMetrics(final IntervalList intervals,
                      final Histogram<Integer> highQualityDepthHistogram,
                      final Histogram<Integer> unfilteredDepthHistogram,
                      final double pctExcludedByAdapter,
                      final double pctExcludedByMapq,
                      final double pctExcludedByDupes,
                      final double pctExcludedByPairing,
                      final double pctExcludedByBaseq,
                      final double pctExcludedByOverlap,
                      final double pctExcludedByCapping,
                      final double pctExcludeTotal,
                      final int coverageCap,
                      final Histogram<Integer> unfilteredBaseQHistogram,
                      final int theoreticalHetSensitivitySampleSize) {
        this.intervals      = intervals.uniqued();
        this.highQualityDepthHistogram = highQualityDepthHistogram;
        this.unfilteredDepthHistogram = unfilteredDepthHistogram;
        this.unfilteredBaseQHistogram = unfilteredBaseQHistogram;
        this.coverageCap    = coverageCap;
        this.theoreticalHetSensitivitySampleSize = theoreticalHetSensitivitySampleSize;

        PCT_EXC_ADAPTER  = pctExcludedByAdapter;
        PCT_EXC_MAPQ     = pctExcludedByMapq;
        PCT_EXC_DUPE     = pctExcludedByDupes;
        PCT_EXC_UNPAIRED = pctExcludedByPairing;
        PCT_EXC_BASEQ    = pctExcludedByBaseq;
        PCT_EXC_OVERLAP  = pctExcludedByOverlap;
        PCT_EXC_CAPPED   = pctExcludedByCapping;
        PCT_EXC_TOTAL    = pctExcludeTotal;

        calculateDerivedFields();
    }

    /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
    @NoMergingIsDerived
    public long GENOME_TERRITORY;
    /** The mean coverage in bases of the genome territory, after all filters are applied. */
    @NoMergingIsDerived
    public double MEAN_COVERAGE;
    /** The standard deviation of coverage of the genome after all filters are applied. */
    @NoMergingIsDerived
    public double SD_COVERAGE;
    /** The median coverage in bases of the genome territory, after all filters are applied. */
    @NoMergingIsDerived
    public double MEDIAN_COVERAGE;
    /** The median absolute deviation of coverage of the genome after all filters are applied. */
    @NoMergingIsDerived
    public double MAD_COVERAGE;

    /** The fraction of aligned bases that were filtered out because they were in reads with mapping quality 0 and the looked like adapter reads. */
    @NoMergingIsDerived
    public double PCT_EXC_ADAPTER;
    /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (lower than MIN_MAPPING_QUALITY). */
    @NoMergingIsDerived
    public double PCT_EXC_MAPQ;
    /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
    @NoMergingIsDerived
    public double PCT_EXC_DUPE;
    /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
    @NoMergingIsDerived
    public double PCT_EXC_UNPAIRED;
    /** The fraction of aligned bases that were filtered out because they were of low base quality (lower than MIN_BASE_QUALITY). */
    @NoMergingIsDerived
    public double PCT_EXC_BASEQ;
    /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
    @NoMergingIsDerived
    public double PCT_EXC_OVERLAP;
    /** The fraction of aligned bases that were filtered out because they would have raised coverage above COVERAGE_CAP. */
    @NoMergingIsDerived
    public double PCT_EXC_CAPPED;
    /** The total fraction of aligned bases excluded due to all filters. */
    @NoMergingIsDerived
    public double PCT_EXC_TOTAL;

    /** The fraction of bases that attained at least 1X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_1X;
    /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_5X;
    /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_10X;
    /** The fraction of bases that attained at least 15X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_15X;
    /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_20X;
    /** The fraction of bases that attained at least 25X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_25X;
    /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_30X;
    /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_40X;
    /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_50X;
    /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_60X;
    /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_70X;
    /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_80X;
    /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_90X;
    /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
    @NoMergingIsDerived
    public double PCT_100X;

    /** The fold over-coverage necessary to raise 80% of bases to the mean coverage level. */
    @NoMergingIsDerived
    public double FOLD_80_BASE_PENALTY;

    /** The fold over-coverage necessary to raise 90% of bases to the mean coverage level. */
    @NoMergingIsDerived
    public double FOLD_90_BASE_PENALTY;

    /** The fold over-coverage necessary to raise 95% of bases to the mean coverage level. */
    @NoMergingIsDerived
    public double FOLD_95_BASE_PENALTY;

    /** The theoretical HET SNP sensitivity. */
    @NoMergingIsDerived
    public double HET_SNP_SENSITIVITY;

    /** The Phred Scaled Q Score of the theoretical HET SNP sensitivity. */
    @NoMergingIsDerived
    public double HET_SNP_Q;

    /**
     * Merges the various PCT_EXC_* metrics.
     * @param other metric to merge into this one.
     *
     * @return result of merging, also known as "this"
     */
    @Override
    public MergeableMetricBase merge(final MergeableMetricBase other) {
        final WgsMetrics otherMetric = (WgsMetrics) other;

        if (highQualityDepthHistogram == null || otherMetric.highQualityDepthHistogram == null ||
                unfilteredDepthHistogram == null || otherMetric.unfilteredDepthHistogram == null) {
            throw new PicardException("Depth histogram is required when deriving metrics.");
        }

        // Union the intervals over which bases are called.  They should have no overlaps!
        // NB: interval lists are already uniqued.
        final long genomeTerritory = this.intervals.getBaseCount() + otherMetric.intervals.getBaseCount();
        this.intervals.addall(otherMetric.intervals.getIntervals());
        this.intervals = this.intervals.uniqued();
        if (this.intervals.getBaseCount() != genomeTerritory) {
            throw new PicardException("Trying to merge WgsMetrics calculated on intervals that overlap.");
        }

        // NB:
        // Since: PCT_EXC_TOTAL     = (totalWithExcludes - thisMetricTotal) / totalWithExcludes;
        // Thus:  totalWithExcludes = total / (1 - PCT_EXC_TOTAL)
        // Proof: Exercise is left to the reader.
        final long thisMetricTotal        = (long) highQualityDepthHistogram.getSum();
        final long otherMetricTotal       = (long) otherMetric.highQualityDepthHistogram.getSum();
        final long total                  = thisMetricTotal + otherMetricTotal;
        final long thisTotalWithExcludes  = (long) (thisMetricTotal / (1.0 - PCT_EXC_TOTAL));
        final long otherTotalWithExcludes = (long) (otherMetricTotal / (1.0 - otherMetric.PCT_EXC_TOTAL));
        final double totalWithExcludes    = thisTotalWithExcludes + otherTotalWithExcludes;

        if (0 < totalWithExcludes) {
            PCT_EXC_DUPE     = (PCT_EXC_DUPE * thisTotalWithExcludes + otherMetric.PCT_EXC_DUPE * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_ADAPTER  = (PCT_EXC_ADAPTER * thisTotalWithExcludes + otherMetric.PCT_EXC_ADAPTER * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_MAPQ     = (PCT_EXC_MAPQ * thisTotalWithExcludes + otherMetric.PCT_EXC_MAPQ * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_UNPAIRED = (PCT_EXC_UNPAIRED * thisTotalWithExcludes + otherMetric.PCT_EXC_UNPAIRED * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_BASEQ    = (PCT_EXC_BASEQ * thisTotalWithExcludes + otherMetric.PCT_EXC_BASEQ * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_OVERLAP  = (PCT_EXC_OVERLAP * thisTotalWithExcludes + otherMetric.PCT_EXC_OVERLAP * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_CAPPED   = (PCT_EXC_CAPPED * thisTotalWithExcludes + otherMetric.PCT_EXC_CAPPED * otherTotalWithExcludes) / totalWithExcludes;
            PCT_EXC_TOTAL    = (totalWithExcludes - total) / totalWithExcludes;
        }

        // do any merging that are dictated by the annotations.
        super.merge(other);

        // merge the histograms
        highQualityDepthHistogram.addHistogram(otherMetric.highQualityDepthHistogram);
        unfilteredDepthHistogram.addHistogram(otherMetric.unfilteredDepthHistogram);
        if (unfilteredBaseQHistogram != null && otherMetric.unfilteredBaseQHistogram != null)
            unfilteredBaseQHistogram.addHistogram(otherMetric.unfilteredBaseQHistogram);
        return this;
    }

    @Override
    public void calculateDerivedFields() {
        if (highQualityDepthHistogram == null || unfilteredDepthHistogram == null) throw new PicardException("Depth histogram is required when deriving metrics.");
        if (unfilteredBaseQHistogram != null && theoreticalHetSensitivitySampleSize <= 0) {
            throw new PicardException("Sample size is required when a baseQ histogram is given when deriving metrics.");
        }

        final long[] depthHistogramArray = new long[coverageCap + 1];

        for (final Histogram.Bin<Integer> bin : highQualityDepthHistogram.values()) {
            final int depth = Math.min((int) bin.getIdValue(), coverageCap);
            depthHistogramArray[depth] += bin.getValue();
        }

        GENOME_TERRITORY = (long) highQualityDepthHistogram.getSumOfValues();
        MEAN_COVERAGE    = highQualityDepthHistogram.getMean();
        SD_COVERAGE      = highQualityDepthHistogram.getStandardDeviation();
        MEDIAN_COVERAGE  = highQualityDepthHistogram.getMedian();
        MAD_COVERAGE     = highQualityDepthHistogram.getMedianAbsoluteDeviation();

        PCT_1X   = MathUtil.sum(depthHistogramArray, 1, depthHistogramArray.length)   / (double) GENOME_TERRITORY;
        PCT_5X   = MathUtil.sum(depthHistogramArray, 5, depthHistogramArray.length)   / (double) GENOME_TERRITORY;
        PCT_10X  = MathUtil.sum(depthHistogramArray, 10, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_15X  = MathUtil.sum(depthHistogramArray, 15, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_20X  = MathUtil.sum(depthHistogramArray, 20, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_25X  = MathUtil.sum(depthHistogramArray, 25, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_30X  = MathUtil.sum(depthHistogramArray, 30, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_40X  = MathUtil.sum(depthHistogramArray, 40, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_50X  = MathUtil.sum(depthHistogramArray, 50, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_60X  = MathUtil.sum(depthHistogramArray, 60, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_70X  = MathUtil.sum(depthHistogramArray, 70, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_80X  = MathUtil.sum(depthHistogramArray, 80, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_90X  = MathUtil.sum(depthHistogramArray, 90, depthHistogramArray.length)  / (double) GENOME_TERRITORY;
        PCT_100X = MathUtil.sum(depthHistogramArray, 100, depthHistogramArray.length) / (double) GENOME_TERRITORY;


        // This roughly measures by how much we must over-sequence so that xx% of bases have coverage at least as deep as the current mean coverage:
        if (highQualityDepthHistogram.getCount() > 0) {
            FOLD_80_BASE_PENALTY = MEAN_COVERAGE / highQualityDepthHistogram.getPercentile(0.2);
            FOLD_90_BASE_PENALTY = MEAN_COVERAGE / highQualityDepthHistogram.getPercentile(0.1);
            FOLD_95_BASE_PENALTY = MEAN_COVERAGE / highQualityDepthHistogram.getPercentile(0.05);
        } else {
            FOLD_80_BASE_PENALTY = 0;
            FOLD_90_BASE_PENALTY = 0;
            FOLD_95_BASE_PENALTY = 0;
        }

        // Get Theoretical Het SNP Sensitivity
        if (unfilteredBaseQHistogram != null && unfilteredDepthHistogram != null) {
            final double[] depthDoubleArray = TheoreticalSensitivity.normalizeHistogram(unfilteredDepthHistogram);
            final double[] baseQDoubleArray = TheoreticalSensitivity.normalizeHistogram(unfilteredBaseQHistogram);
            HET_SNP_SENSITIVITY = TheoreticalSensitivity.hetSNPSensitivity(depthDoubleArray, baseQDoubleArray, theoreticalHetSensitivitySampleSize, LOG_ODDS_THRESHOLD);
            HET_SNP_Q = QualityUtil.getPhredScoreFromErrorProbability((1 - HET_SNP_SENSITIVITY));
        }
    }
}