package net.sf.picard.analysis;

import net.sf.picard.metrics.MetricBase;

/** Metrics for evaluating the performance of whole genome sequencing experiments. */
public class WgsMetrics extends MetricBase {
    /** The number of non-N bases in the genome reference over which coverage will be evaluated. */
    public long GENOME_TERRITORY;
    /** The mean coverage in bases of the genome territory, after all filters are applied. */
    public double MEAN_COVERAGE;
    /** The standard deviation of coverage of the genome after all filters are applied. */
    public double SD_COVERAGE;
    /** The median coverage in bases of the genome territory, after all filters are applied. */
    public double MEDIAN_COVERAGE;
    /** The median absolute deviation of coverage of the genome after all filters are applied. */
    public double MAD_COVERAGE;

    /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20). */
    public double PCT_EXC_MAPQ;
    /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
    public double PCT_EXC_DUPE;
    /** The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair. */
    public double PCT_EXC_UNPAIRED;
    /** The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20). */
    public double PCT_EXC_BASEQ;
    /** The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads. */
    public double PCT_EXC_OVERLAP;
    /** The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x). */
    public double PCT_EXC_CAPPED;
    /** The total fraction of aligned bases excluded due to all filters. */
    public double PCT_EXC_TOTAL;

    /** The fraction of bases that attained at least 5X sequence coverage in post-filtering bases. */
    public double PCT_5X;
    /** The fraction of bases that attained at least 10X sequence coverage in post-filtering bases. */
    public double PCT_10X;
    /** The fraction of bases that attained at least 20X sequence coverage in post-filtering bases. */
    public double PCT_20X;
    /** The fraction of bases that attained at least 30X sequence coverage in post-filtering bases. */
    public double PCT_30X;
    /** The fraction of bases that attained at least 40X sequence coverage in post-filtering bases. */
    public double PCT_40X;
    /** The fraction of bases that attained at least 50X sequence coverage in post-filtering bases. */
    public double PCT_50X;
    /** The fraction of bases that attained at least 60X sequence coverage in post-filtering bases. */
    public double PCT_60X;
    /** The fraction of bases that attained at least 70X sequence coverage in post-filtering bases. */
    public double PCT_70X;
    /** The fraction of bases that attained at least 80X sequence coverage in post-filtering bases. */
    public double PCT_80X;
    /** The fraction of bases that attained at least 90X sequence coverage in post-filtering bases. */
    public double PCT_90X;
    /** The fraction of bases that attained at least 100X sequence coverage in post-filtering bases. */
    public double PCT_100X;
}
