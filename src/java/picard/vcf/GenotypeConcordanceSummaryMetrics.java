package picard.vcf;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Class that holds summary metrics about Genotype Concordance
 *
 * @author George Grant
 */
public class GenotypeConcordanceSummaryMetrics extends MetricBase {
    GenotypeConcordanceSummaryMetrics(final String eventType, final ConcordanceResults concordanceResults) {
        this.EVENT_TYPE = eventType;
        this.TRUTH_SAMPLE_NAME = concordanceResults.getTruthSample();
        this.CALL_SAMPLE_NAME = concordanceResults.getCallSample();
        this.HET_SENSITIVITY = concordanceResults.hetSensitivity();
        this.HET_PPV = concordanceResults.hetPpv();
        this.HOMVAR_SENSITIVITY = concordanceResults.homVarSensitivity();
        this.HOMVAR_PPV = concordanceResults.homVarPpv();
        this.VAR_SENSITIVITY = concordanceResults.varSensitivity();
        this.VAR_PPV = concordanceResults.varPpv();
        this.FALSE_POS_PER_MB = concordanceResults.numFalsePositives();
    }

    /** The type of event SNP/Indel */
    public String EVENT_TYPE;

    /** The name of the 'truth' sample */
    public String TRUTH_SAMPLE_NAME;

    /** The name of the 'call' sample */
    public String CALL_SAMPLE_NAME;

    /** The het sensitivity */
    public double HET_SENSITIVITY;

    /** The het ppv (positive predictive value) */
    public double HET_PPV;

    /** The homVar sensitivity */
    public double HOMVAR_SENSITIVITY;

    /** The homVar ppv (positive predictive value) */
    public double HOMVAR_PPV;

    /** The var sensitivity */
    public double VAR_SENSITIVITY;

    /** The var ppv (positive predictive value) */
    public double VAR_PPV;

    /** The number of false positives per megabase of sequence */
    public double FALSE_POS_PER_MB;
}
