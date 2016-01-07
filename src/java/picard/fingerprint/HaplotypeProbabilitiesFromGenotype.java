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

import static picard.util.MathUtil.*;

/**
 * Represents a set of HaplotypeProbabilities that were derived from a single SNP
 * genotype at a point in time.
 */
public class HaplotypeProbabilitiesFromGenotype extends HaplotypeProbabilities {
    private final Snp snp;
    private final double[] likelihoods;

    public HaplotypeProbabilitiesFromGenotype(final Snp snp, final HaplotypeBlock haplotypeBlock,
                                              final double AA, final double Aa, final double aa) {
        super(haplotypeBlock);
        this.snp = snp;
        this.likelihoods = new double[] {AA, Aa, aa};
    }

    /** Returns the SNP who's genotype was used to construct the likelihoods. */
    @Override public Snp getRepresentativeSnp() { return snp; }


    // TODO: this can't be right in general. At least one needs to divide by the prior to set things straight.
    // TODO: The only saving grace is that this is normally used for cases where the priors are large and similar to each other.


    // simply returns the _likelihoods_ that were passed into the constructor.
    public double[] getLikelihoods() {
        return likelihoods;
    }

	@Override
	public void merge(final HaplotypeProbabilities other) {
		if (!this.getHaplotype().equals(other.getHaplotype())) {
			throw new IllegalArgumentException("Mismatched haplotypes in call to HaplotypeProbabilities.merge(): " +
					getHaplotype() + ", " + other.getHaplotype());
		}

		if (! (other instanceof HaplotypeProbabilitiesFromGenotype)) {
			throw new IllegalArgumentException("Can only merge HaplotypeProbabilities of same class.");
		}

		this.likelihoods[0] = this.likelihoods[0] * other.getLikelihoods()[0];
		this.likelihoods[1] = this.likelihoods[1] * other.getLikelihoods()[1];
		this.likelihoods[2] = this.likelihoods[2] * other.getLikelihoods()[2];
	}
}
