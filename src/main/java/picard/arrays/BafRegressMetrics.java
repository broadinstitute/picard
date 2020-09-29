package picard.arrays;

import htsjdk.samtools.metrics.MetricBase;

public class BafRegressMetrics extends MetricBase {
    /** The sample name */
    public String SAMPLE;

    /** The estimate of contamination from the model (on the 0.0-1.0 scale) */
    public double ESTIMATE;

    /** The standard error of the estimate */
    public double STDERR;

    /** The test statistic for the estimate */
    public double TVAL;

    /** The p-value of the estimate */
    public double PVAL;

    /** The call rate of the sample (number of non-missing genotypes) */
    public double CALL_RATE;

    /** The number of homozygous genotypes used to fit the model */
    public int NHOM;
}
