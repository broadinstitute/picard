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
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.analysis.MergeableMetricBase;
import picard.util.help.HelpConstants;

/**
 * Created by farjoun on 6/26/18.
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
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

    public ErrorMetric(final String covariate, final long totalBases) {
        this.TOTAL_BASES = totalBases;
        this.COVARIATE = covariate;
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
