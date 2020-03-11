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
 * Represents the results of comparing evidence for a single haplotype locus between
 * two sources of evidence (typically external genotyping data vs. sequencing data.).
 *
 * @author Tim Fennell 
 */
class LocusResult implements Comparable<LocusResult> {
    private final Snp snp;
    private final DiploidGenotype expectedGenotype;
    private final DiploidGenotype mostLikelyGenotype;
    private final int allele1Count;
    private final int allele2Count;
    private final double lodGenotype;
    private final double lodExpectedSampleTumorNormal; //LOD assuming that the first sample is from a tumor and the second is from the normal
    private final double lodExpectedSampleNormalTumor; //LOD assuming that the first sample is from the normal and the second is from a tumor

    private final double lExpectedSample;   // log probability of expected sample
    private final double lRandomSample;     // log probability of random sample

    LocusResult(final Snp snp, final DiploidGenotype expectedGenotype, final DiploidGenotype mostLikelyGenotype,
                final int allele1Count, final int allele2Count, final double lodGenotype,
                final double lExpectedSample, final double lRandomSample,
                final double lodGenotypeTumorNormal, final double lodGenotypeNormalTumor) {
        this.snp = snp;
        this.expectedGenotype = expectedGenotype;
        this.mostLikelyGenotype = mostLikelyGenotype;
        this.allele1Count = allele1Count;
        this.allele2Count = allele2Count;
        this.lodGenotype = lodGenotype;
        this.lExpectedSample = lExpectedSample;
        this.lRandomSample = lRandomSample;
        this.lodExpectedSampleTumorNormal = lodGenotypeTumorNormal;
        this.lodExpectedSampleNormalTumor = lodGenotypeNormalTumor;
    }

    public Snp getSnp() { return snp; }
    public DiploidGenotype getExpectedGenotype() { return expectedGenotype; }
    public DiploidGenotype getMostLikelyGenotype() { return mostLikelyGenotype; }
    public int getAllele1Count() { return allele1Count; }
    public int getAllele2Count() { return allele2Count; }
    public double getLodGenotype() { return lodGenotype; }
    public double getLodExpectedSampleNormalTumor() { return lodExpectedSampleNormalTumor; }
    public double getLodExpectedSampleTumorNormal() { return lodExpectedSampleTumorNormal; }
    public double lExpectedSample() { return lExpectedSample; }
    public double lRandomSample() { return lRandomSample; }

    @Override
    public int compareTo(final LocusResult that) {
        return this.snp.compareTo(that.snp);
    }
}
