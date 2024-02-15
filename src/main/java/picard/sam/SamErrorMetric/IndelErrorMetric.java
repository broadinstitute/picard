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
 * Metric to be used for InDel errors
 */
@DocumentedFeature(groupName = HelpConstants.DOC_CAT_METRICS, summary = HelpConstants.DOC_CAT_METRICS_SUMMARY)
public class IndelErrorMetric extends BaseErrorMetric {

    /**
     * The number of insertions. Note: This is not the number of bases that have been inserted.
     */
    @MergeByAdding
    public long NUM_INSERTIONS = 0;

    /**
     * The number of inserted bases.
     */
    @MergeByAdding
    public long NUM_INSERTED_BASES = 0;

    /**
     * The (phred) rate of insertions.
     */
    @NoMergingIsDerived
    public int INSERTIONS_Q = 0;

    /**
     * The number of deletions. Note: This is not the number of bases that have been deleted.
     */
    @MergeByAdding
    public long NUM_DELETIONS = 0L;

    /**
     * The number of deleted bases.
     */
    @MergeByAdding
    public long NUM_DELETED_BASES = 0;

    /**
     * The (phred) rate of deletions.
     */
    @NoMergingIsDerived
    public int DELETIONS_Q = 0;

    @Override
    public void calculateDerivedFields() {
        this.INSERTIONS_Q = computeQScore(NUM_INSERTIONS, TOTAL_BASES);
        this.DELETIONS_Q = computeQScore(NUM_DELETIONS, TOTAL_BASES);

    }

    public IndelErrorMetric(final String covariate,
                            final long nTotalBases,
                            final long nInserts,
                            final long nInsertedBases,
                            final long nDeletions,
                            final long nDeletedBases) {
        super(covariate, nTotalBases, nInsertedBases + nDeletedBases);

        this.NUM_INSERTIONS = nInserts;
        this.NUM_INSERTED_BASES = nInsertedBases;
        this.NUM_DELETIONS = nDeletions;
        this.NUM_DELETED_BASES = nDeletedBases;
    }

    // needed for reading in a metric from a file
    public IndelErrorMetric() {
    }
}
