package picard.vcf;

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * Class that holds detail metrics about Genotype Concordance
 *
 * @author George Grant
 */
public class GenotypeConcordanceDetailMetrics extends MetricBase {
    /** The type of the event (i.e. either SNP or INDEL) */
    public VariantContext.Type VARIANT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE;

    /** The state of the 'truth' sample (i.e. HOM_REF, HET_REF_VAR1, HET_VAR1_VAR2...) */
    public GenotypeConcordanceStates.TruthState TRUTH_STATE;

    /** The state of the 'call' sample (i.e. HOM_REF, HET_REF_VAR1...) */
    public GenotypeConcordanceStates.CallState CALL_STATE;

    /** The number of events of type TRUTH_STATE and CALL_STATE for the EVENT_TYPE and SAMPLEs */
    public long COUNT;
}
