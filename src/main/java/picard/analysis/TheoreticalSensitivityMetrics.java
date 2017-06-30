package picard.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Created by fleharty on 6/30/17.
 */
public class TheoreticalSensitivityMetrics extends MetricBase {
    /** Theoretical sensitivity at 0.1% allele fraction */
    public double SENSITIVITY_AT_0_1 = 0.0;
    /** Theoretical sensitivity at 0.1% allele fraction */

    /** Theoretical sensitivity at 0.5% allele fraction */
    public double SENSITIVITY_AT_0_5 = 0.0;

    /** Theoretical sensitivity at 1% allele fraction */
    public double SENSITIVITY_AT_01 = 0.0;

    /** Theoretical sensitivity at 2% allele fraction */
    public double SENSITIVITY_AT_02 = 0.0;

    /** Theoretical sensitivity at 5% allele fraction */
    public double SENSITIVITY_AT_05 = 0.0;

    /** Theoretical sensitivity at 10% allele fraction */
    public double SENSITIVITY_AT_10 = 0.0;

    /** Theoretical sensitivity at 30% allele fraction */
    public double SENSITIVITY_AT_30 = 0.0;

    /** Theoretical sensitivity at 50% allele fraction */
    public double SENSITIVITY_AT_50 = 0.0;
}