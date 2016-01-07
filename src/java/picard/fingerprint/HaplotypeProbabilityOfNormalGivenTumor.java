/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

/**
 * A wrapper class for any HaplotypeProbabilities instance that will assume that the given evidence is that of a tumor sample and
 * provide an hp for the normal sample that tumor came from. This models possible loss of hetrozygosity where het genotypes
 * turn into a homozygous genotype with probability pLoH.
 *
 * The shortcoming of this model is that we assume that the events are all independent, but this way they are allowed.
 *
 * @author farjoun
 */

public class HaplotypeProbabilityOfNormalGivenTumor extends HaplotypeProbabilities {

    private final double[][] transitionMatrix;
    private final HaplotypeProbabilities hpOfTumor;

    public HaplotypeProbabilityOfNormalGivenTumor(final HaplotypeProbabilities hpOfTumor, final double pLoH) {
        super(hpOfTumor.getHaplotype());

        this.hpOfTumor = hpOfTumor;
        transitionMatrix = new double[][]{
                //This is P(g_t|g_n)
                //tumor genotype are the columns.
                {1,               0,        0},  //normal is hom_ref => tumor must be the same
                {pLoH / 2, 1 - pLoH, pLoH / 2},  //normal is het     => tumor might transit
                {0,               0,        1}}; //normal is hom_var => tumor must be the same
    }

    // This function needs to be overridden since we want likelihood to mean the probability of the
    // data given a particular _normal_ genotype, however, the likelihood as given is that where the
    // genotype is of the tumor (if that's what the data was measuring)

    // P(D_t|g_n) = \sum_{g_n}  P(D_t|g_t,g_n) = \sum P(D_t|g_t) P(g_t|g_n) = hpOfTumor.getLikelihoods() * transitionMatrix

    @Override
    public double[] getLikelihoods() {
        final double[] asTumorLikelihoods = new double[3];
        final double[] asNormalLikelihoods = hpOfTumor.getLikelihoods();
        for (final Genotype g_n : Genotype.values()) {
            for (final Genotype g_t : Genotype.values()) {
                asTumorLikelihoods[g_t.v] += asNormalLikelihoods[g_n.v] * transitionMatrix[g_n.v][g_t.v];
            }
        }
        return asTumorLikelihoods;
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
