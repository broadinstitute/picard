package picard.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Created by fleharty on 6/30/17.
 */
public class TheoreticalSensitivityMetrics extends MetricBase {
    /** Theoretical Het Sensitivity for log-odss of 6.2 **/
    public double HET_SENSITIVITY_LOD6_2;

    /** Theoretical Het Sensitivity for log-odss of 3.0 **/
    public double HET_SENSITIVITY_LOD3_0;
}