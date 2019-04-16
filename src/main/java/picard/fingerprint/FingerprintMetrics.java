/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

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
 */
public class FingerprintMetrics extends MetricBase {
    /* The Sample alias taken from RG header or #CROME line */
    public String SAMPLE_ALIAS;

    /* The originating file (if available) for this fingerprint */
    public String SOURCE;

    /* Additional information about the fingerprint */
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

    /* The Chi-squared pvalue of the observed counts vector relative to the expected counts, (3x2 table) */
    public double CHI_SQUARED_PVALUE;

    /* The log10 of the Chi-squared pvalue*/
    public double LOG10_CHI_SQUARED_PVALUE;

    /* The categorical cross entropy of the counts of genotypes relative to expected (big is bad)*/
    public double CROSS_ENTROPY_LOD;

    /* The Chi-squared pvalue for the number of HETs and HOMs relative to the expected counts (2x2 table) */
    public double HET_CHI_SQUARED_PVALUE;

    /* The log10 of the Chi-squared pvalue for the number of HETs and HOMs */
    public double HET_LOG10_CHI_SQUARED_PVALUE;

    /* The categorical cross entropy of the counts of HETs and HOMs relative to the expected counts (big is bad) */
    public double HET_CROSS_ENTROPY_LOD;

    /* The Chi-squared pvalue for the number of HOM-Allele1s and HOM_Allele2s relative to the expected counts (2x2 table) */
    public double HOM_CHI_SQUARED_PVALUE;

    /* The log10 of the Chi-squared pvalue for the number of HOM-Allele1s and HOM_Allele2s */
    public double HOM_LOG10_CHI_SQUARED_PVALUE;

    /* The categorical cross entropy of the HOM-Allele1s and HOM_Allele2s relative to the expected counts (big is bad)*/
    public double HOM_CROSS_ENTROPY_LOD;

    /* The fingerprinting LOD score this sample gets when compared to itself (big is good)*/
    public double LOD_SELF_CHECK;

    /* The difference in fingerprinting LOD between LOD_SELF_CHECK and the LOD score found when
    fingerprinting it against a random permutation of the probablity vectors (PLs) in each of the haplotypes. */
    public double DISCRIMINATORY_POWER;
}
