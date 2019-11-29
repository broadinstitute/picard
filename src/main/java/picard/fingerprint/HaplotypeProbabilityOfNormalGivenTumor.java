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

package picard.fingerprint;

import htsjdk.samtools.util.CollectionUtil;

import java.util.Map;

/**
 * A wrapper class for any HaplotypeProbabilities instance that will assume that the given evidence is that of a tumor sample and
 * provide an hp for the normal sample that tumor came from. This models possible loss of hetrozygosity where het genotypes
 * turn into a homozygous genotype with probability pLoH.
 * <p>
 * The shortcoming of this model is that we assume that the events are all independent, but this way they are allowed.
 *
 * @author farjoun
 */

public class HaplotypeProbabilityOfNormalGivenTumor extends HaplotypeProbabilities {
    static private class TransitionMatrix {
        private final double[][] transitionMatrix;

        public TransitionMatrix(double pLoH) {
            transitionMatrix = new double[][]{
                    //This is P(g_t|g_n)
                    //tumor genotype are the columns.
                    {1, 0, 0},  //normal is hom_ref => tumor must be the same
                    {pLoH / 2, 1 - pLoH, pLoH / 2},  //normal is het     => tumor may have transitioned
                    {0, 0, 1}}; //normal is hom_var => tumor must be the same
        }

        public double[][] getTransitionMatrix() {
            return transitionMatrix;
        }
    }

    private final HaplotypeProbabilities hpOfTumor;
    static private Map<Double, TransitionMatrix> transitionMatrixMap = new CollectionUtil.DefaultingMap<>(TransitionMatrix::new, true);

    final private double[][] transitionMatrix;

    public HaplotypeProbabilityOfNormalGivenTumor(final HaplotypeProbabilities hpOfTumor, final double pLoH) {
        super(hpOfTumor.getHaplotype());
        transitionMatrix = transitionMatrixMap.get(pLoH).getTransitionMatrix();
        this.hpOfTumor = hpOfTumor;
    }

    // This function needs to be overridden since we want likelihood to mean the probability of the
    // data given a particular _normal_ genotype, however, the likelihood as given is that where the
    // genotype is of the tumor (if that's what the data was measuring)

    // P(D_t|g_n) = \sum_{g_t}  P(D_t|g_t,g_n)
    //            = \sum P(D_t|g_t, g_n) P(g_t|g_n)
    //            = \sum P(D_t|g_t) P(g_t|g_n)
    //            = hpOfTumor.getLikelihoods() * transitionMatrix
    // where the * operator is understood as linear algebra operation.

    @Override
    public double[] getLikelihoods() {
        final double[] normalHaplotypeLikelihoods = new double[3];
        final double[] tumorHaplotypeLikelihoods = hpOfTumor.getLikelihoods();
        for (final Genotype g_n : Genotype.values()) {
            normalHaplotypeLikelihoods[g_n.v] = 0D;
            for (final Genotype g_t : Genotype.values()) {
                normalHaplotypeLikelihoods[g_n.v] += tumorHaplotypeLikelihoods[g_t.v] * transitionMatrix[g_n.v][g_t.v];
            }
        }
        return normalHaplotypeLikelihoods;
    }

    @Override
    public Snp getRepresentativeSnp() {
        return hpOfTumor.getRepresentativeSnp();
    }

    @Override
    public void merge(final HaplotypeProbabilities ignored) {
        throw new IllegalArgumentException("Cannot merge HaplotypeProbabilityOfNormalGivenTumor. Merge the underlying object and create a new wrapper.");
    }

    @Override
    public int getObsAllele1() {
        return hpOfTumor.getObsAllele1();
    }

    @Override
    public int getObsAllele2() {
        return hpOfTumor.getObsAllele2();
    }

    @Override
    public int getTotalObs() {
        return hpOfTumor.getTotalObs();
    }

    @Override
    public boolean hasEvidence() {
        return hpOfTumor.hasEvidence();
    }
}
