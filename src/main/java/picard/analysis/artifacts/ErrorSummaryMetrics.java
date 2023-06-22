package picard.analysis.artifacts;

import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.MergeableMetricBase;
import picard.util.help.HelpConstants;

/**
 * Summary metrics produced by {@link CollectSequencingArtifactMetrics} as a roll up of the
 * context-specific error rates, to provide global error rates per type of base substitution.
 *
 * Errors are normalized to the lexically lower reference base and summarized together. E.g.
 * G>T is converted to C>A and merged with data from C>A for reporting.
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class ErrorSummaryMetrics extends MergeableMetricBase {
    /** The reference base (or it's complement). */
    @MergeByAssertEquals public char REF_BASE;

    /** The alternative base (or it's complement). */
    @MergeByAssertEquals public char ALT_BASE;

    /** A single string representing the substition from REF_BASE to ALT_BASE for convenience. */
    @MergeByAssertEquals public String SUBSTITUTION;

    /** The number of reference bases observed. */
    @MergeByAdding public long REF_COUNT;

    /** The number of alt bases observed. */
    @MergeByAdding public long ALT_COUNT;

    /** The rate of the substitution in question. */
    @NoMergingIsDerived public double SUBSTITUTION_RATE;

    @Override
    public void calculateDerivedFields() {
        final double total = REF_COUNT + ALT_COUNT;
        this.SUBSTITUTION_RATE = (total == 0) ? 0 : ALT_COUNT / total;
    }
}
