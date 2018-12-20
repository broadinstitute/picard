package picard.fingerprint;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Class for holding
 */
public class FingerprintMetrics extends MetricBase {
    /* The Sample alias taken from RG header or #CROME line */
    public String SAMPLE_ALIAS;
    /* the originating file (if available) for this fingerprint */
    public String SOURCE;
    /* additional information about the fingerprint */
    public String INFO;


    public long HAPLOTYPES;
    public long HAPLOTYPES_WITH_EVIDENCE;
    public long DEFINITE_GENOTYPES;
    public long NUM_HOM_REF;
    public long NUM_HET;
    public long NUM_HOM_VAR;
    public double CHI_SQUARED_PVALUE;
    public double LOG10_CHI_SQUARED_PVALUE;
    public double CROSS_ENTROPY_LOD;
    public double HET_CHI_SQUARED_PVALUE;
    public double HET_LOG10_CHI_SQUARED_PVALUE;
    public double HET_CROSS_ENTROPY_LOD;
    public double HOM_CHI_SQUARED_PVALUE;
    public double HOM_LOG10_CHI_SQUARED_PVALUE;
    public double HOM_CROSS_ENTROPY_LOD;
    public double DISCRIMINATORY_POWER;
    public double LOD_SELF_CHECK;

}
