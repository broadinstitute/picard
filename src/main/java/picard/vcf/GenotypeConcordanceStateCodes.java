package picard.vcf;

/**
 * Created by kbergin on 7/30/15.
 */
public enum GenotypeConcordanceStateCodes {
    /**
     * These codes allow for constants in the truth state enum and the call state enum to be compared.
     * The final INCOMPARABLE_CODE is for Call States that are not reflected in the Truth State, so they can never match.
     */
        MISSING_CODE,
        HOM_REF_CODE,
        HET_REF_VAR1_CODE,
        HET_VAR1_VAR2_CODE,
        HOM_VAR1_CODE,
        NO_CALL_CODE,
        LOW_GQ_CODE,
        LOW_DP_CODE,
        VC_FILTERED_CODE,
        GT_FILTERED_CODE,
        IS_MIXED_CODE,
        INCOMPARABLE_CODE;
}
