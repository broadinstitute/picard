package picard.arrays;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
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
