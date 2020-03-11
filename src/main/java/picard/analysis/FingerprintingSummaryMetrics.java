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
 * Summary fingerprinting metrics and statistics about the comparison of the sequence data
 * from a single read group (lane or index within a lane) vs. a set of known genotypes for
 * the expected sample.
 *
 * @author Tim Fennell
 */
public class FingerprintingSummaryMetrics extends MetricBase {
    /** The read group from which sequence data was drawn for comparison. */
    public String READ_GROUP;

    /** The sample whose known genotypes the sequence data was compared to. */
    public String SAMPLE;

    /** The Log Likelihood of the sequence data given the expected sample's genotypes. */
    public double LL_EXPECTED_SAMPLE;

    /** The Log Likelihood of the sequence data given a random sample from the human population. */
    public double LL_RANDOM_SAMPLE;

    /**
     * The LOD for Expected Sample vs. Random Sample.  A positive LOD indicates that the sequence data
     * is more likely to come from the expected sample vs. a random sample from the population, by LOD logs.
     * I.e. a value of 6 indicates that the sequence data is 1,000,000 more likely to come from the expected
     * sample than from a random sample.  A negative LOD indicates the reverse - that the sequence data is more
     * likely to come from a random sample than from the expected sample.
     */
    public double LOD_EXPECTED_SAMPLE;

    /** The number of haplotypes that had expected genotypes to compare to. */
    public int HAPLOTYPES_WITH_GENOTYPES;

    /**
     * The subset of genotyped haplotypes for which there was sufficient sequence data to
     * confidently genotype the haplotype. Note: all haplotypes with sequence coverage contribute to the
     * LOD score, even if they cannot be "confidently checked" individually.
     * */
    public int HAPLOTYPES_CONFIDENTLY_CHECKED;

    /** The subset of confidently checked haplotypes that match the expected genotypes. */
    public int HAPLOTYPES_CONFIDENTLY_MATCHING;

    /** num of hets, observed as homs with LOD > threshold */
    public int HET_AS_HOM;

    /** num of homs, observed as hets with LOD > threshold */
    public int HOM_AS_HET;

    /** num of homs, observed as other homs with LOD > threshold */
    public int HOM_AS_OTHER_HOM;
}
