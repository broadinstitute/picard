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

import java.util.Arrays;

import static java.lang.Math.log10;
import static picard.util.MathUtil.multiply;
import static picard.util.MathUtil.pNormalizeVector;

/**
 * Abstract class for storing and calculating various likelihoods and probabilities
 * for haplotype alleles given evidence.
 *
 * @author Tim Fennell
 */
public abstract class HaplotypeProbabilities {

    private final HaplotypeBlock haplotypeBlock;

    protected HaplotypeProbabilities(final HaplotypeBlock haplotypeBlock) {
        this.haplotypeBlock = haplotypeBlock;
    }

    /** Returns the haplotype for which the probabilities apply. */
    public HaplotypeBlock getHaplotype() {
        return this.haplotypeBlock;
    }

    public double [] getPriorProbablities(){
        return getHaplotype().getHaplotypeFrequencies();
    }

    /** Returns the probabilities, in order, of the AA, Aa and aa haplotypes.

     * Mathematically, this is P(H | D, F) where and H is the vector of possible haplotypes {AA,Aa,aa}.
     * D is the data seen by the class, and
     * F is the population frequency of each genotype.
     */
    /** Returns the posterior probabilities using the population frequency as a prior. */
    public double[] getPosteriorProbabilities() {
        return pNormalizeVector(multiply(getLikelihoods(), getPriorProbablities()));}

    /**
     * Returns the likelihoods, in order, of the AA, Aa and aa haplotypes given the evidence
     *
     * Mathematically this is P(evidence | haplotype) where haplotype={AA,Aa,aa}.
     */
    public abstract double[] getLikelihoods();

    public double[] getLogLikelihoods() {
        final double[] likelihoods = getLikelihoods();
        final double[] lLikelihoods = new double [likelihoods.length];
        for (int i = 0; i < likelihoods.length; ++i) {
            lLikelihoods[i] = Math.log10(likelihoods[i]);
        }
        return lLikelihoods;
    }

    /**
     * Returns a representative SNP for this haplotype. Different subclasses may implement this in
     * different ways, but should do so in a deterministic/repeatable fashion.
     */
    public abstract Snp getRepresentativeSnp();

    /**
     * Returns the number of observations of alleles supporting the first/major haplotype allele.
     * Strictly this doesn't make sense for all subclasses, but it's nice to have it part of the API so
     * a default implementation is provided here.
     * @return int
     */
    public int getObsAllele1() { return 0; }

    /**
     * Returns the number of observations of alleles supporting the second/minor haplotype allele.
     * Strictly this doesn't make sense for all subclasses, but it's nice to have it part of the API so
     * a default implementation is provided here.
     * @return int
     */
    public int getObsAllele2() { return 0; }

    /**
     * Returns the total number of observations of any allele.
     * Strictly this doesn't make sense for all subclasses, but it's nice to have it part of the API so
     * a default implementation is provided here.
     * @return int
     */
    public int getTotalObs() { return 0; }

    /** Returns true if evidence has been added, false if the probabilities are just the priors. */
    public boolean hasEvidence() {
        return true;
    }

	/** Merges in the likelihood information from the supplied haplotype probabilities object. */
	public abstract void merge(final HaplotypeProbabilities other);

    /**
     * Returns the index of the highest probability which can then be used to look up
     * DiploidHaplotypes or DiploidGenotypes as appropriate.
     */
    int getMostLikelyIndex() {
        final double[] probs = getPosteriorProbabilities();

        if (probs[0] > probs[1] && probs[0] > probs[2]) return 0;
        else if (probs[1] > probs[2]) return 1;
        else return 2;
    }

    /** Gets the most likely haplotype given the probabilities. */
    public DiploidHaplotype getMostLikelyHaplotype() {
        return DiploidHaplotype.values()[getMostLikelyIndex()];
    }

    /** Gets the genotype for this Snp given the most likely haplotype. */
    public DiploidGenotype getMostLikelyGenotype(final Snp snp) {
        assertSnpPartOfHaplotype(snp);
        return snp.getGenotype(getMostLikelyHaplotype());
    }

    /** Throws an exception if the passed SNP is not part of this haplotype. */
    void assertSnpPartOfHaplotype(final Snp snp) {
        if (!this.haplotypeBlock.getSnps().contains(snp)) {
            throw new IllegalArgumentException("Snp " + snp + " does not belong to haplotype " + this.haplotypeBlock);
        }
    }

    /** This function returns the scaled probability of the evidence collected
     * given a vector of priors on the haplotype using the internal likelihood, which may be
     * scaled by an unknown factor. This factor causes the result to be scaled, hence the name.
     *
     * Mathematically:
     *
     * P(Evidence| P(h_i)=F_i) = \sum_i P(Evidence | h_i) P(h_i)
     *                         = \sum_i P(Evidence | h_i) F_i
     *                         = c * \sum_i Likelihood_i * F_i
     *
     * Here, h_i are the three possible haplotypes, F_i are the given priors, and Likelihood_i
     * are the stored likelihoods which are scaled from the actually likelihoods by an unknown
     * factor, c. Note that the calculation ignores the internal haplotype probabilities (i.e. priors)
     *
     * @param genotypeFrequencies vector of (possibly scaled) probabilities of the three haplotypes
     * @return P(evidence | P_h)) / c
     */

    public double scaledEvidenceProbabilityUsingGenotypeFrequencies(final double[] genotypeFrequencies) {
        final double[] likelihoods = getLikelihoods();
        assert (genotypeFrequencies.length == likelihoods.length);

        double result = 0;
        for (int i = 0; i < likelihoods.length; ++i) {
            result += likelihoods[i] * genotypeFrequencies[i];
        }
        return result;
    }

    public double shiftedLogEvidenceProbabilityUsingGenotypeFrequencies(final double[] genotypeFrequencies) {
        return log10(scaledEvidenceProbabilityUsingGenotypeFrequencies(genotypeFrequencies));
    }

    /**
     * returns the log-probability the evidence, using as priors the posteriors of another object
     *
     * @param otherHp an additional HaplotypeProbabilities object representing the same underlying HaplotypeBlock
     * @return  log10( P(evidence| P(h_i)=P(h_i|otherHp) ) + c where c is an unknown constant
     */
    public double shiftedLogEvidenceProbabilityGivenOtherEvidence(final HaplotypeProbabilities otherHp) {
        if (!this.haplotypeBlock.equals(otherHp.getHaplotype())) {
            throw new IllegalArgumentException("Haplotypes are from different HaplotypeBlocks!");
        }
        /** Get the posterior from the other otherHp. Use this posterior as the prior to calculate probability.
         *
         *   P(hap|x,y) = P(x|hap,y) P(hap|y) / P(x|y)
         *              = P(x | hap) * P(hap | y) / P(x)
         *                likelihood * other.posterior
         *
         *              = P(x|hap) P(y|hap) P(hap)/P(x)P(y)
         *              = A P(x| hap) P(y| hap) P(hap)  # where A is an unknown scaling factor
         */
        return shiftedLogEvidenceProbabilityUsingGenotypeFrequencies(otherHp.getPosteriorProbabilities());
    }

    /**
     * Returns log (p(evidence)) + c assuming that the prior on haplotypes is given by
     * the internal haplotypeFrequencies
     */
    public double shiftedLogEvidenceProbability() {
        return shiftedLogEvidenceProbabilityUsingGenotypeFrequencies(getPriorProbablities());
    }

    /** Returns the LOD score between the most probable haplotype and the second most probable. */
    public double getLodMostProbableGenotype() {
        final double[] probs = getPosteriorProbabilities();
        final double[] logs = new double[probs.length];
        for (int i = 0; i < probs.length; ++i) {
            logs[i] = log10(probs[i]);
        }

        Arrays.sort(logs);
        return logs[2] - logs[1];
    }

    /** Log10(P(evidence| haplotype)) for the 3 different possible haplotypes
     * {aa, ab, bb}
     */

    //an enum whose only role in life is to help iterate over the 3 possible diploid genotypes
    protected enum Genotype {
        HOM_ALLELE1(0),
        HET_ALLELE12(1),
        HOM_ALLELE2(2);

        int v; //value is the number of chromosomes in the genotype that have ALLELE2.

        Genotype(final int v) {
            this.v = v;
        }
    }
}
