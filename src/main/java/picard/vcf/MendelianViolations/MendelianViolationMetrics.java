/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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
package picard.vcf.MendelianViolations;

import picard.analysis.replicates.MergeableMetricBase;
import picard.pedigree.Sex;

/**
 * Describes the type and number of mendelian violations found within a Trio.
 */
public class MendelianViolationMetrics extends MergeableMetricBase {
    public static String getExtension() {
        return "mendelian_violation_metrics";
    }

    /** The family ID assigned to the trio for which these metrics are calculated.*/
    @MergeByAssertEquals
    public String FAMILY_ID;
    /** The ID of the mother within the trio.*/
    @MergeByAssertEquals
    public String MOTHER;
    /** The ID of the father within the trio.*/
    @MergeByAssertEquals
    public String FATHER;
    /** The ID of the offspring within the trio.*/
    @MergeByAssertEquals
    public String OFFSPRING;
    /** The sex of the offspring. */
    @MergeByAssertEquals
    public Sex OFFSPRING_SEX = null;
    /** The number of biallelic, SNP sites at which all relevant samples exceeded the minimum genotype quality and depth and at least one of the samples was variant. */
    @MergeByAdding
    public long NUM_VARIANT_SITES;
    /** The number of diploid sites at which a potential de-novo mutation was observed (i.e. both parents are hom-ref, offspring is not hom-ref. */
    @MergeByAdding
    public long NUM_DIPLOID_DENOVO;
    /** The number of sites at which both parents are homozygous for a non-reference allele and the offspring is heterozygous. */
    @MergeByAdding
    public long NUM_HOMVAR_HOMVAR_HET;
    /** The number of sites at which the one parent is homozygous reference, the other homozygous variant and the offspring is homozygous. */
    @MergeByAdding
    public long NUM_HOMREF_HOMVAR_HOM;
    /** The number of sites at which one parent is homozygous, the other is heterozygous and the offspring is the alternative homozygote. */
    @MergeByAdding
    public long NUM_HOM_HET_HOM;
    /** The number of sites at which the offspring is haploid, the parent is homozygous reference and the offspring is non-reference. */
    @MergeByAdding
    public long NUM_HAPLOID_DENOVO;
    /** The number of sites at which the offspring is haploid and exhibits a reference allele that is not present in the parent. */
    @MergeByAdding
    public long NUM_HAPLOID_OTHER;
    /** The number of otherwise unclassified events. */
    @MergeByAdding
    public long NUM_OTHER;
    /** The total of all mendelian violations observed. */
    @NoMergingIsDerived
    public long TOTAL_MENDELIAN_VIOLATIONS;

    @Override
    public void calculateDerivedFields() {
        TOTAL_MENDELIAN_VIOLATIONS = NUM_DIPLOID_DENOVO + NUM_HOMVAR_HOMVAR_HET +
                NUM_HOMREF_HOMVAR_HOM + NUM_HOM_HET_HOM + NUM_HAPLOID_DENOVO +
                NUM_HAPLOID_OTHER + NUM_OTHER;
    }
}
