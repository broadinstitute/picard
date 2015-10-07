package picard.vcf;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.variant.variantcontext.VariantContext;
import picard.vcf.GenotypeConcordanceStates.*;

import java.util.Map;

/**
 * Class that holds metrics about the Genotype Concordance contingency tables.
 */
public class GenotypeConcordanceContingencyMetrics extends MetricBase {
    /**
     * Empty constructor - needed for unit tests
     */
    public GenotypeConcordanceContingencyMetrics() {
    }

    GenotypeConcordanceContingencyMetrics(final VariantContext.Type variantType, final GenotypeConcordanceCounts concordanceCounts,
                                          final String truthSample, final String callSample, final boolean missingSitesFlag) {
        this.VARIANT_TYPE = variantType;
        this.TRUTH_SAMPLE = truthSample;
        this.CALL_SAMPLE = callSample;

        final GenotypeConcordanceSchemeFactory schemeFactory = new GenotypeConcordanceSchemeFactory();
        final GenotypeConcordanceScheme scheme = schemeFactory.getScheme(missingSitesFlag);
        scheme.validateScheme();
        concordanceCounts.validateCountsAgainstScheme(scheme);

        Map<ContingencyState, Long> counts = concordanceCounts.getContingencyStateCounts(scheme);
        this.TP_COUNT = counts.get(ContingencyState.TP);
        this.TN_COUNT = counts.get(ContingencyState.TN);
        this.FP_COUNT = counts.get(ContingencyState.FP);
        this.FN_COUNT = counts.get(ContingencyState.FN);
        this.EMPTY_COUNT = counts.get(ContingencyState.EMPTY);
    }

    /** The type of the event (i.e. either SNP or INDEL) */
    public VariantContext.Type VARIANT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE;

    /** The TP (true positive) count across all variants */
    public long TP_COUNT;

    /** The TN (true negative) count across all variants */
    public long TN_COUNT;

    /** The FP (false positive) count across all variants */
    public long FP_COUNT;

    /** The FN (false negative) count across all variants */
    public long FN_COUNT;

    /** The empty (no contingency info) count across all variants */
    public long EMPTY_COUNT;
}
