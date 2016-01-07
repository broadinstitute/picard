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

import picard.util.MathUtil;
import java.util.Arrays;
import static java.lang.Math.log10;

/**
 * Represents the probability of the underlying haplotype using log likelihoods as the basic datum for each of the SNPs. By convention the
 * alleles stored for each SNP are in phase.
 *
 *
 * @author Tim Fennell
 * @author Yossi Farjoun
 */
abstract class HaplotypeProbabilitiesUsingLogLikelihoods extends HaplotypeProbabilities {

    // some derived classes might need to incorporate accumulated data before logLikelihood is usable.
    // use the getter to allow these classes to calculate the likelihood from the data.
    private final double[] loglikelihoods = new double[Genotype.values().length];

    public HaplotypeProbabilitiesUsingLogLikelihoods(final HaplotypeBlock haplotypeBlock) {
        super(haplotypeBlock);
    }

    /** Simple returns the SNP from the haplotype that has the lowest genome coordinate. */
    @Override
    public Snp getRepresentativeSnp() {
        return getHaplotype().getFirstSnp();
    }

    @Override
    public boolean hasEvidence() {
        final double [] ll = this.getLogLikelihoods();
        return ll[Genotype.HOM_ALLELE1.v]  != 0 ||
               ll[Genotype.HET_ALLELE12.v] != 0 ||
               ll[Genotype.HOM_ALLELE2.v]  != 0;
    }

    /**
     * Merges information from another haplotype probabilities object for the same haplotype into
     * this object. Useful for when probabilities need to be merged to levels higher than the
     * read group, e.g. the sample or individual.
     *
     * @param other Another haplotype probabilities object to merge in (must of the the same class and for the same HaplotypeBlock)
     *
     */
    @Override
    public void merge(final HaplotypeProbabilities other) {
        if (!this.getHaplotype().equals(other.getHaplotype())) {
            throw new IllegalArgumentException("Mismatched haplotypes in call to HaplotypeProbabilities.merge(): " +
                    getHaplotype() + ", " + other.getHaplotype());
        }

        if (!(other instanceof HaplotypeProbabilitiesUsingLogLikelihoods)) {
            throw new IllegalArgumentException("Can only merge HaplotypeProbabilities of same class.");
        }

        final HaplotypeProbabilitiesUsingLogLikelihoods o = (HaplotypeProbabilitiesUsingLogLikelihoods) other;

        setLogLikelihoods(MathUtil.sum(getLogLikelihoods(), o.getLogLikelihoods()));
    }

    /**
     * Returns the posterior probability of the haplotypes given the evidence (uses the internal prior)
     *
     */
    public double[] getPosteriorProbabilities() {
        return MathUtil.pNormalizeLogProbability(getShiftedLogPosterior());
    }

    /** Makes a copy of the loglikelihoods array and applies the priors.
     * returns log10( P(haplotype | evidence) ) + C where C is unknown.
     * One can recover C by normalizing, but this might be unneeded depending on the application
     * uses Bayes P(m|x)=P(x|m)*P(m)/P(x) but then doesn't divide by P(x)
     *
     * uses the internal prior as P(m)
     * */
    private double[] getShiftedLogPosterior() {
        final double[] ll = this.getLogLikelihoods();
        final double[] shiftedLogPosterior = new double [Genotype.values().length];
        final double[] haplotypeFrequencies = getPriorProbablities();
        for (final Genotype g : Genotype.values()){
            shiftedLogPosterior[g.v] = ll[g.v] + log10(haplotypeFrequencies[g.v]);
        }
        return shiftedLogPosterior;
    }

    /**
     * Converts the loglikelihoods into linear-space.
     */
    @Override
    public double[] getLikelihoods() {
        return MathUtil.pNormalizeLogProbability(getLogLikelihoods());
    }

    /**
     * Since this class uses loglikelihoods natively, we override and return the native variable
     */
    @Override
    public double[] getLogLikelihoods() {
        return this.loglikelihoods;
    }

    public void setLogLikelihoods(final double[] ll) {
        assert (ll.length == Genotype.values().length);

        System.arraycopy(ll, 0, loglikelihoods, 0, ll.length);

    }
    /**
     * Overridden to calculate the LOD from the loglikelihoods instead of the probabilities
     * because it will allow for more accurate calculation before overflowing.
     */
    @Override
    public double getLodMostProbableGenotype() {
        final double[] logs = getShiftedLogPosterior();
        Arrays.sort(logs);
        return logs[Genotype.values().length-1] - logs[Genotype.values().length-2];
    }

}
