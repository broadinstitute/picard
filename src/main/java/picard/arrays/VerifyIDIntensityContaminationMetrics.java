package picard.arrays;

import htsjdk.samtools.metrics.MetricBase;

public class VerifyIDIntensityContaminationMetrics extends MetricBase {
    /** The ID of this entry */
    public int ID;

    /** The percent mixture (contamination) of the sample for ID */
    public double PCT_MIX;

    /** The log likelihood */
    public double LLK;

    /** The log likelihood 0 */
    public double LLK0;
}
