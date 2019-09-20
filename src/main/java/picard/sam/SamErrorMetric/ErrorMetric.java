/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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


package picard.sam.SamErrorMetric;

import htsjdk.samtools.util.QualityUtil;
import picard.analysis.MergeableMetricBase;

/**
 * Created by farjoun on 6/26/18.
 */
public class ErrorMetric extends MergeableMetricBase {


    protected static double PRIOR_ERROR = 0.001;
    /**
     * The value of the covariate define the bases included in this metric
     */
    @MergingIsManual
    public String COVARIATE;
    /**
     * The total number of bases included in the calculation of this metric
     */
    @MergeByAdding
    public long TOTAL_BASES;
    /**
     * Number of bases which are skipped because they overlap with a SNP variant site
     */

    @MergeByAdding
    public long SKIPPED_SNPS;

    /**
     * Number of insertions or deletions which are skipped because they overlap with an indel variant site. Note that
     * this is not the number of bases that are skipped, i.e. each insertion or deletion is counted only once.
     */
    @MergeByAdding
    public long SKIPPED_INDELS;

    public ErrorMetric(final String covariate, final long totalBases, final long skippedSNPs, final long skippedIndels) {
        this.TOTAL_BASES = totalBases;
        this.COVARIATE = covariate;
        this.SKIPPED_SNPS = skippedSNPs;
        this.SKIPPED_INDELS = skippedIndels;
    }

    // required to enable reading metric from a file.
    public ErrorMetric() {
    }

    public static void setPriorError(double priorError) {
        PRIOR_ERROR = priorError;
    }

    /**
     * compute a qscore given the number of errors and the total number of bases.
     * Uses a false count of 1 int the numerator and 1/PRIOR_ERROR in the denominator.
     */
    protected int computeQScore(final long numberOfErrors) {
        return computeQScore(numberOfErrors, TOTAL_BASES);
    }

    /**
     * compute a qscore given the number of errors and the total number of bases.
     * Uses a false count of 1 int the denominator and 1 in the numerator.
     */
    protected int computeQScore(final long numberOfErrors, final long nTotalBases) {
        return QualityUtil.getPhredScoreFromErrorProbability(
                (numberOfErrors + PRIOR_ERROR) / (nTotalBases + 1.0D));
    }
}
