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

import htsjdk.samtools.util.QualityUtil;
import picard.util.MathUtil;

import static java.lang.Math.log10;

/**
 * Represents the probability of the underlying haplotype of the contaminating sample given the data.
 * By convention the alleles stored for each SNP are in phase.
 *
 * @author Yossi Farjoun
 */

public class HaplotypeProbabilitiesFromContaminatorSequence extends HaplotypeProbabilitiesFromSequence {

    public double contamination;

    // for each model (contGenotype, mainGenotype) there's a likelihood of the data. These need to be collected separately
    // and only collated once all the data is in.
    double[][] likelihoodMap = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

    public HaplotypeProbabilitiesFromContaminatorSequence(final HaplotypeBlock haplotypeBlock, final double contamination) {
        super(haplotypeBlock);

        assert (contamination <= 1.0);
        assert (contamination >= 0.0);

        this.contamination = contamination;
    }

    /**
     * Adds a base observation with the observed quality to the evidence for this haplotype
     * based on the fact that the SNP is part of the haplotype.
     */
    public void addToProbs(final Snp snp, final byte base, final byte qual) {
        assertSnpPartOfHaplotype(snp);

        // Skip bases that don't match either expected allele for this SNP
        final boolean altAllele;
        if (base == snp.getAllele1()) {
            this.obsAllele1++;
            altAllele = false;
        } else if (base == snp.getAllele2()) {
            this.obsAllele2++;
            altAllele = true;
        } else {
            this.obsAlleleOther++;
            return;
        }
        final double pErr = QualityUtil.getErrorProbabilityFromPhredScore(qual);

        // we need to keep the 9 models separate until all the reads have been seen.
        // Once we have seen all the reads, we add across the mainGeno and the three likelihoods
        // are the likelihoods of the contaminator, the main sample's genotype needs to be summed over in each case:

        for (final Genotype contGeno : Genotype.values()) {
            for (final Genotype mainGeno : Genotype.values()) {
                //theta is the expected frequency of the alternate allele
                final double theta = 0.5 * ((1 - contamination) * mainGeno.v + contamination * contGeno.v);
                likelihoodMap[contGeno.v][mainGeno.v] *= (( altAllele ? theta : (1 - theta)) * (1 - pErr) +
                                                          (!altAllele ? theta : (1 - theta)) * pErr);
            }
        }
    }

    //a function needed to update the logLikelihoods from the likelihoodMap.
    private void updateLikelihoods() {
        final double[] ll = new double[Genotype.values().length];
        for (final Genotype contGeno : Genotype.values()) {
            // p(a | g_c) = \sum_g_m { P(g_m) \prod_i P(a_i| g_m, g_c)}
            ll[contGeno.v] = log10(MathUtil.sum(MathUtil.multiply(this.getPriorProbablities(), likelihoodMap[contGeno.v])));
        }
        setLogLikelihoods(ll);
    }

    @Override
    public void merge(final HaplotypeProbabilities other) {
        super.merge(other);

        if (!this.getHaplotype().equals(other.getHaplotype())) {
            throw new IllegalArgumentException("Mismatched haplotypes in call to HaplotypeProbabilities.merge(): " +
                    getHaplotype() + ", " + other.getHaplotype());
        }

        if (!(other instanceof HaplotypeProbabilitiesFromContaminatorSequence)) {
            throw new IllegalArgumentException("Can only merge HaplotypeProbabilities of same class.");
        }

        final HaplotypeProbabilitiesFromContaminatorSequence o = (HaplotypeProbabilitiesFromContaminatorSequence) other;
        if (o.contamination != this.contamination) {
            throw new IllegalArgumentException("Can only merge HaplotypeProbabilitiesFromContaminatorSequence with the same contamination value.");
        }

        for (final Genotype contGeno : Genotype.values()) {
            this.likelihoodMap[contGeno.v] = MathUtil.multiply(this.likelihoodMap[contGeno.v], o.likelihoodMap[contGeno.v]);
        }
    }

    @Override
    public double[] getLogLikelihoods() {
        updateLikelihoods();
        return super.getLogLikelihoods();
    }
}
