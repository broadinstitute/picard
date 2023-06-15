/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.util.help.HelpConstants;

/**
 * An error metric for the errors in bases.
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class BaseErrorMetric extends ErrorMetric {
    /** The number of bases that disagree with the reference */
    @MergeByAdding
    public long ERROR_BASES;

    /** The (phred) rate of bases that disagree with the reference */
    @NoMergingIsDerived
    public int Q_SCORE;

    @Override
    public void calculateDerivedFields() {
        this.Q_SCORE = computeQScore(ERROR_BASES);
    }

    public BaseErrorMetric(final String covariate, final long totalBases, final long errorBases) {
        super(covariate, totalBases);
        this.ERROR_BASES = errorBases;
    }

    // needed for reading in a metric from a file
    public BaseErrorMetric() {
    }
}
