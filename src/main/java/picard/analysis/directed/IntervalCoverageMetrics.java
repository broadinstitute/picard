package picard.analysis.directed;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Created by farjoun on 8/22/17.
 */
public class IntervalCoverageMetrics extends MetricBase{

    public int PF_UNIQUE_TANDEM_READS = 0;
    public int PF_UNIQUE_OUTIE_READS = 0;
    public int PF_UNIQUE_INNI_READS = 0;
    public int PF_UNIQUE_READS = 0;
    public int PF_READS = 0;
    public int TOTAL_READS = 0;
}
