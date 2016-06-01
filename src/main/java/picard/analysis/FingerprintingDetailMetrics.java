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

package picard.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * Detailed metrics about an individual SNP/Haplotype comparison within a fingerprint comparison.
 *
 * @author Tim Fennell
 */
public class FingerprintingDetailMetrics extends MetricBase {
    /** The sequencing read group from which sequence data was fingerprinted. */
    public String READ_GROUP;

    /** The name of the sample who's genotypes the sequence data was compared to. */
    public String SAMPLE;

    /**
     * The name of a representative SNP within the haplotype that was compared. Will usually be the
     * exact SNP that was genotyped externally.
     */
    public String SNP;

    /** The possible alleles for the SNP. */
    public String SNP_ALLELES;

    /** The chromosome on which the SNP resides. */
    public String CHROM;

    /** The position of the SNP on the chromosome. */
    public int    POSITION;

    /** The expected genotype of the sample at the SNP locus. */
    public String EXPECTED_GENOTYPE;

    /** The most likely genotype given the observed evidence at the SNP locus in the sequencing data. */
    public String OBSERVED_GENOTYPE;

    /** The LOD score for OBSERVED_GENOTYPE vs. the next most likely genotype in the sequencing data. */
    public double LOD;

    /** The number of observations of the first, or A,  allele of the SNP in the sequencing data. */
    public int OBS_A;

    /** The number of observations of the second, or B, allele of the SNP in the sequencing data. */
    public int OBS_B;
}
