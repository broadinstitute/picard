/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.variant.variantcontext.Allele;
import picard.util.MathUtil;

import java.util.List;

/**
 * Represents the likelihood of the HaplotypeBlock given the GenotypeLikelihoods (GL field from a VCF, which is actually a log10-likelihood)
 * for each of the SNPs in that block. By convention the alleles stored for each SNP are in phase.
 *
 * @author Yossi Farjoun
 */
public class HaplotypeProbabilitiesFromGenotypeLikelihoods extends HaplotypeProbabilitiesUsingLogLikelihoods {

    public HaplotypeProbabilitiesFromGenotypeLikelihoods(final HaplotypeBlock haplotypeBlock) {
        super(haplotypeBlock);
    }

    /**
     * Adds a base observation with the observed quality to the evidence for this haplotype
     * based on the fact that the SNP is part of the haplotype.
     *
     * @param snp The snp in the haplotypeblock to which the likelihoods belong
     * @param alleles the (ordered) alleles to which the biallelic genotype likelihoods correspond. So that if the alleles are [A,B], the
     * @param logGenotypeLikelihoods correspond to the logLikelihoods of [AA, AB, BB]. Log is assumed to be in base 10.
     */

    public void addToLogLikelihoods(final Snp snp, final List<Allele> alleles, final double logGenotypeLikelihoods[]) {
        assertSnpPartOfHaplotype(snp);

        // only allow biallelic snps
        assert (logGenotypeLikelihoods.length == Genotype.values().length);
        assert (alleles.size() == 2);

        //make sure that alleles are comparable to SNPs
        for (int i = 0; i < 2; i++) {
            assert (alleles.get(i).getBases().length == 1);
        }

        final byte allele1 = alleles.get(0).getBases()[0];
        final byte allele2 = alleles.get(1).getBases()[0];

        // alleles as given might be swapped with alleles in haplotype block.
        // if that is the case, swap them around.

        if (snp.getAllele1() == allele1 &&
                snp.getAllele2() == allele2) {
            setLogLikelihoods(MathUtil.sum(getLogLikelihoods(),logGenotypeLikelihoods));
            return;
        }
        if (snp.getAllele2() == allele1 &&
                snp.getAllele1() == allele2) {
            final double [] ll = getLogLikelihoods();
            ll[Genotype.HOM_ALLELE1.v]  += logGenotypeLikelihoods[Genotype.HOM_ALLELE2.v];
            ll[Genotype.HET_ALLELE12.v] += logGenotypeLikelihoods[Genotype.HET_ALLELE12.v];
            ll[Genotype.HOM_ALLELE2.v]  += logGenotypeLikelihoods[Genotype.HOM_ALLELE1.v];

            setLogLikelihoods(ll);
            return;
        }

        // if we are here it means that there was a mismatch in alleles...
        assert true;
    }
}
