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
import htsjdk.utils.ValidationUtils;
import htsjdk.variant.variantcontext.VariantContext;
import picard.util.MathUtil;

/**
 * Represents the probability of the underlying haplotype of the contaminating sample given the data.
 * By convention the alleles stored for each SNP are in phase.
 *
 * @author Yossi Farjoun
 */

public class HaplotypeProbabilitiesFromContaminatedArraysVC extends HaplotypeProbabilitiesUsingLogLikelihoods {

    public double contamination;

    // for each model (contGenotype, mainGenotype) there's a likelihood of the data. These need to be collected separately
    // and only collated once all the data is in.
    private final double[][] likelihoodMap = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};
    private boolean valuesNeedUpdating = true;

    public HaplotypeProbabilitiesFromContaminatedArraysVC(final HaplotypeBlock haplotypeBlock, final double contamination) {
        super(haplotypeBlock);

        ValidationUtils.validateArg(contamination <= 1.0, () -> "contamination must be <=1, found " + contamination);
        ValidationUtils.validateArg(contamination >= 0.0, () -> "contamination must be >=0, found " + contamination);

        this.contamination = contamination;
    }

    public HaplotypeProbabilitiesFromContaminatedArraysVC(HaplotypeProbabilitiesFromContaminatedArraysVC other) {
        super(other);

        contamination = other.contamination;
        for (final Genotype g : Genotype.values()) {
            System.arraycopy(other.likelihoodMap[g.v], 0, likelihoodMap[g.v], 0, NUM_GENOTYPES);
        }
    }

    /**
     * Adds a base observation with the observed quality to the evidence for this haplotype
     * based on the fact that the SNP is part of the haplotype.
     */
    public void addToProbs(final VariantContext vc) {
    }

    //a function needed to update the logLikelihoods from the likelihoodMap.
    private void updateLikelihoods() {

    }

    @Override
    public HaplotypeProbabilitiesFromContaminatedArraysVC deepCopy() {
        return new HaplotypeProbabilitiesFromContaminatedArraysVC(this);
    }

    @Override
    public HaplotypeProbabilitiesFromContaminatedArraysVC merge(final HaplotypeProbabilities other) {
        super.merge(other);

        if (!this.getHaplotype().equals(other.getHaplotype())) {
            throw new IllegalArgumentException("Mismatched haplotypes in call to HaplotypeProbabilities.merge(): " +
                    getHaplotype() + ", " + other.getHaplotype());
        }

        if (!(other instanceof HaplotypeProbabilitiesFromContaminatedArraysVC)) {
            throw new IllegalArgumentException("Can only merge HaplotypeProbabilities of same class.");
        }

        final HaplotypeProbabilitiesFromContaminatedArraysVC o = (HaplotypeProbabilitiesFromContaminatedArraysVC) other;
        if (o.contamination != this.contamination) {
            throw new IllegalArgumentException("Can only merge HaplotypeProbabilitiesFromContaminatorSequence with the same contamination value.");
        }

        for (final Genotype contGeno : Genotype.values()) {
            this.likelihoodMap[contGeno.v] = MathUtil.multiply(this.likelihoodMap[contGeno.v], o.likelihoodMap[contGeno.v]);
        }
        valuesNeedUpdating = true;
        return this;
    }

    @Override
    public double[] getLikelihoods() {
        updateLikelihoods();
        return super.getLikelihoods();
    }

    @Override
    public double[] getLogLikelihoods() {
        updateLikelihoods();
        return super.getLogLikelihoods();
    }
}
