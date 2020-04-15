/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import htsjdk.utils.ValidationUtils;
import picard.util.MathUtil;

import java.util.Arrays;

/**
 * Represents the probability of the underlying haplotype using logLikelihoods as the basic datum for each of the SNPs. By convention the
 * alleles stored for each SNP are in phase.
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
 */
abstract class HaplotypeProbabilitiesUsingLogLikelihoods extends HaplotypeProbabilities {

    // some derived classes might need to incorporate accumulated data before logLikelihood is usable.
    // use the getter to allow these classes to calculate the likelihood from the data.
    private final double[] loglikelihoods = new double[NUM_GENOTYPES];

    private boolean likelihoodsNeedUpdating = true;

    // stored in order to reduce computation we store these partial results.
    // they need to be recalculated if loglikelihoodNeedsUpdating
    private final double[] likelihoods = new double[NUM_GENOTYPES];
    private final double[] posteriorProbabilities = new double[NUM_GENOTYPES];

    //normalized (likelihood * prior / normalization_factor)
    private final double[] shiftedLogPosteriors = new double[NUM_GENOTYPES];

    public HaplotypeProbabilitiesUsingLogLikelihoods(final HaplotypeBlock haplotypeBlock) {
        super(haplotypeBlock);
    }

    public HaplotypeProbabilitiesUsingLogLikelihoods(final HaplotypeProbabilitiesUsingLogLikelihoods other) {
        super(other.getHaplotype());
        System.arraycopy(other.loglikelihoods, 0, loglikelihoods, 0, NUM_GENOTYPES);
        System.arraycopy(other.likelihoods, 0, likelihoods, 0, NUM_GENOTYPES);
        System.arraycopy(other.posteriorProbabilities, 0, posteriorProbabilities, 0, NUM_GENOTYPES);
        System.arraycopy(other.shiftedLogPosteriors, 0, shiftedLogPosteriors, 0, NUM_GENOTYPES);
        likelihoodsNeedUpdating = other.likelihoodsNeedUpdating;
    }

    /**
     * Simple returns the SNP from the haplotype that has the lowest genome coordinate.
     */
    @Override
    public Snp getRepresentativeSnp() {
        return getHaplotype().getFirstSnp();
    }

    @Override
    public boolean hasEvidence() {
        return Arrays.stream(getLogLikelihoods()).anyMatch(d -> d != 0);
    }

    /**
     * Merges information from another haplotype probabilities object for the same haplotype into
     * this object. Useful for when probabilities need to be merged to levels higher than the
     * read group, e.g. the sample or individual.
     *
     * @param other Another haplotype probabilities object to merge in (must of the the same class and for the same HaplotypeBlock)
     * @return
     */
    @Override
    public HaplotypeProbabilitiesUsingLogLikelihoods merge(final HaplotypeProbabilities other) {
        if (!this.getHaplotype().equals(other.getHaplotype())) {
            throw new IllegalArgumentException("Mismatched haplotypes in call to HaplotypeProbabilities.merge(): " +
                    getHaplotype() + ", " + other.getHaplotype());
        }

        if (!(other instanceof HaplotypeProbabilitiesUsingLogLikelihoods)) {
            throw new IllegalArgumentException(String.format("Can only merge HaplotypeProbabilities of same class. Found %s and %s",
                    this.getClass().toString(), other.getClass().toString()));
        }

        final HaplotypeProbabilitiesUsingLogLikelihoods o = (HaplotypeProbabilitiesUsingLogLikelihoods) other;

        setLogLikelihoods(MathUtil.sum(getLogLikelihoods(), o.getLogLikelihoods()));
        return this;
    }

    /**
     * Returns the posterior probability of the haplotypes given the evidence (uses the internal prior)
     */
    protected double[] getPosteriorProbabilities0() {
        return MathUtil.pNormalizeLogProbability(getShiftedLogPosterior0());
    }

    /**
     * getter for posteriorProbs
     */
    @Override
    public double[] getPosteriorProbabilities() {
        updateDependentValues();
        return posteriorProbabilities;
    }

    /**
     * Makes a copy of the loglikelihoods array and applies the priors.
     * returns log10( P(haplotype | evidence) ) + C where C is unknown.
     * One can recover C by normalizing, but this might be unneeded depending on the application
     * uses Bayes P(m|x)=P(x|m)*P(m)/P(x) but then doesn't divide by P(x)
     * <p>
     * uses the internal prior as P(m)
     */
    private double[] getShiftedLogPosterior0() {
        final double[] ll = this.getLogLikelihoods();
        final double[] shiftedLogPosterior = new double[NUM_GENOTYPES];
        final double[] haplotypeFrequencies = getPriorProbablities();
        for (final Genotype g : Genotype.values()) {
            shiftedLogPosterior[g.v] = ll[g.v] + Math.log10(haplotypeFrequencies[g.v]);
        }
        return shiftedLogPosterior;
    }

    /**
     * getter for shiftedLogPosterior
     */
    private double[] getShiftedLogPosterior() {
        updateDependentValues();
        return shiftedLogPosteriors;
    }

    /**
     * Converts the loglikelihoods into linear-space with normalizing.
     */
    @Override
    public double[] getLikelihoods() {
        updateDependentValues();
        return likelihoods;
    }

    public double[] getLikelihoods0() {
        return MathUtil.pNormalizeLogProbability(this.loglikelihoods);
    }

    /**
     * Since this class uses log-rawLikelihoods natively, we override and return the native variable
     */
    @Override
    public double[] getLogLikelihoods() {
        return this.loglikelihoods;
    }

    public void setLogLikelihoods(final double[] ll) {
        ValidationUtils.validateArg(ll.length == NUM_GENOTYPES,
                () -> "logLikelihood must have length 3, found " + ll.length);

//        protect from underflow
        double max = MathUtil.max(ll);
        final double[] maxRemoved = MathUtil.sum(ll, -max);

        double sum = MathUtil.sum(MathUtil.getProbabilityFromLog(maxRemoved));
        // normalize log rawLikelihoods:
        System.arraycopy(MathUtil.sum(maxRemoved, -Math.log10(sum)), 0, loglikelihoods, 0, NUM_GENOTYPES);

        likelihoodsNeedUpdating = true;
        updateDependentValues();
    }

    /**
     * Overridden to calculate the LOD from the loglikelihoods instead of the probabilities
     * because it will allow for more accurate calculation before overflowing.
     */
    @Override
    public double getLodMostProbableGenotype() {
        final double[] logs = getShiftedLogPosterior();
        double biggest = -Double.MAX_VALUE;
        double secondBiggest = biggest;

        for (double prob : logs) {
            if (prob > biggest) {
                secondBiggest = biggest;
                biggest = prob;
                continue;
            }
            if (prob > secondBiggest) {
                secondBiggest = prob;
            }
        }
        return biggest - secondBiggest;
    }

    private void updateDependentValues() {
        if (likelihoodsNeedUpdating) {
            System.arraycopy(getLikelihoods0(), 0, likelihoods, 0, NUM_GENOTYPES);
            System.arraycopy(getShiftedLogPosterior0(), 0, shiftedLogPosteriors, 0, NUM_GENOTYPES);
            System.arraycopy(getPosteriorProbabilities0(), 0, posteriorProbabilities, 0, NUM_GENOTYPES);

            likelihoodsNeedUpdating = false;
        }
    }
}
