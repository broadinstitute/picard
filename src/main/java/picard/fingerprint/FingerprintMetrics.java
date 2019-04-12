package picard.fingerprint;

import htsjdk.samtools.metrics.MetricBase;
import picard.analysis.FingerprintingDetailMetrics;
import picard.analysis.FingerprintingSummaryMetrics;

/**
 * Class for holding metrics on a single fingerprint.
 * Note: this is distinct from {@link FingerprintingDetailMetrics} and
 * {@link FingerprintingSummaryMetrics} in that it is calculated on a single fingerprint,
 * and attempts to describe how likely that fingerprint is to have arisen from an actual sample
 * as opposed to having artifacts such as contamination, or strong bias towards homozygous genotypes.
 *
 */
public class FingerprintMetrics extends MetricBase {
    /* The Sample alias taken from RG header or #CROME line */
    public String SAMPLE_ALIAS;
    /* the originating file (if available) for this fingerprint */
    public String SOURCE;
    /* additional information about the fingerprint */
    public String INFO;

    /* Number of haplotypes examined */
    public long HAPLOTYPES;

    /* Number of haplotypes that had evidence in the source file */
    public long HAPLOTYPES_WITH_EVIDENCE;

    /* Number of haplotypes that had enough evidence to make a definite genotype call */
    public long DEFINITE_GENOTYPES;

    /* Number of major allele homozygous calls*/
    public long NUM_HOM_ALLELE1;

    /* Number of definite heterozygous calls*/
    public long NUM_HET;

    /* Number of definite minor allele homozygous calls*/
    public long NUM_HOM_ALLELE2;

    /* The Chi-squared pvalue of the counts vector relative to the expected counts, (3x2 table) */
    public double CHI_SQUARED_PVALUE;

    /* The log10 of the Chi-squared pvalue*/
    public double LOG10_CHI_SQUARED_PVALUE;

    /* The categorical cross entropy of the counts of genotypes relative to expected */
    public double CROSS_ENTROPY_LOD;

    /* The Chi-squared pvalue for the number of HETs and HOMs relative to the expected counts (2x2 table)*/
    public double HET_CHI_SQUARED_PVALUE;

    /* The log10 of the Chi-squared pvalue for the number of HETs and HOMs */
    public double HET_LOG10_CHI_SQUARED_PVALUE;

    /* The categorical cross entropy of the counts of HETs and HOMs relative to the expected counts */
    public double HET_CROSS_ENTROPY_LOD;

    /* The Chi-squared pvalue for the number of HETs and HOMs relative to the expected counts (2x2 table)*/
    public double HOM_CHI_SQUARED_PVALUE;
    /* How many definite major allele homozygous calls were there*/
    public double HOM_LOG10_CHI_SQUARED_PVALUE;
    /* How many definite major allele homozygous calls were there*/
    /* How many definite major allele homozygous calls were there*/
    public double HOM_CROSS_ENTROPY_LOD;
    /* How many definite major allele homozygous calls were there*/
    public double DISCRIMINATORY_POWER;
    /* How many definite major allele homozygous calls were there*/
    public double LOD_SELF_CHECK;

}
