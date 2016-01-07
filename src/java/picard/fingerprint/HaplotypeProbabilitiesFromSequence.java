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

import htsjdk.samtools.util.QualityUtil;
import static java.lang.Math.log10;

/**
 * Represents the probability of the underlying haplotype given the data. By convention the
 * alleles stored for each SNP are in phase.
 *
 * @author Tim Fennell
 */
public class HaplotypeProbabilitiesFromSequence extends HaplotypeProbabilitiesUsingLogLikelihoods {
    protected int obsAllele1, obsAllele2, obsAlleleOther;

    public HaplotypeProbabilitiesFromSequence(final HaplotypeBlock haplotypeBlock) {
        super(haplotypeBlock);
    }

    @Override
    public boolean hasEvidence() {
        return super.hasEvidence() || obsAllele1 > 0 || obsAllele2 > 0;
    }

    /**
     * Adds a base observation with the observed quality to the evidence for this haplotype
     * based on the fact that the SNP is part of the haplotype.
     *
     * @param snp The snp in the HaplotypeBlock to which evidence is being added
     * @param base the base observed
     * @param qual the quality of the observed base
     */
    public void addToProbs(final Snp snp, final byte base, final byte qual) {
        assertSnpPartOfHaplotype(snp);
        final double [] ll = getLogLikelihoods();
        final double pError = QualityUtil.getErrorProbabilityFromPhredScore(qual);
        // Skip bases that don't match either expected allele for this SNP
        if (base == snp.getAllele1()) {
            obsAllele1++;
            for (final Genotype g:Genotype.values()){
                final double pAlt = g.v / 2d;
                ll[g.v] += log10((1d - pAlt) * (1d - pError) + pAlt * pError);
            }

        } else if (base == snp.getAllele2()) {
            obsAllele2++;
            for (final Genotype g:Genotype.values()){
                final double pAlt = 1 - g.v / 2d;
                ll[g.v] += log10((1d - pAlt) * (1d - pError) + pAlt * pError);
            }
        } else {
            obsAlleleOther++;
        }
        //technically not needed since we were changing the actual array, but good practice perhaps.
        setLogLikelihoods(ll);
    }

    /**
     * Merges information from another haplotype probabilities object for the same haplotype into
     * this object. Useful for when probabilities need to  be merged to levels higher than the
     * read group, e.g. the sample or individual.
     *
     * @param other Another haplotype probabilities object to merge in
     */
	@Override
	public void merge(final HaplotypeProbabilities other) {
        super.merge(other);

		if (!this.getHaplotype().equals(other.getHaplotype())) {
			throw new IllegalArgumentException("Mismatched haplotypes in call to HaplotypeProbabilities.merge(): " +
					getHaplotype() + ", " + other.getHaplotype());
		}

		if (! (other instanceof HaplotypeProbabilitiesFromSequence)) {
			throw new IllegalArgumentException("Can only merge() HaplotypeProbabilities of same class: Tried to merge a " +
                    this.getClass().getName() + " with a " + other.getClass().getName() +"." );
		}

		final HaplotypeProbabilitiesFromSequence o = (HaplotypeProbabilitiesFromSequence) other;
        this.obsAllele1     += o.obsAllele1;
        this.obsAllele2     += o.obsAllele2;
        this.obsAlleleOther += o.obsAlleleOther;
	}

    /** Returns the number of bases/reads that support the first allele. */
    @Override public int getObsAllele1() {
        return obsAllele1;
    }

    /** Returns the number of bases/reads that support the second allele. */
    @Override public int getObsAllele2() {
        return obsAllele2;
    }

    /** Gets the total number of observations presented at this locus. */
    @Override
    public int getTotalObs() { return obsAllele1 + obsAllele2 + obsAlleleOther; }

    /* Returns the faction of base observations that were presented that were from an allele other than the two expected ones. */
    public double getFractionUnexpectedAlleleObs() {
        return obsAlleleOther / (double) (getTotalObs());
    }
}
