/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import htsjdk.samtools.SamPairUtil.PairOrientation;
import picard.metrics.MultilevelMetrics;

/**
 * Metrics about the insert size distribution of a paired-end library, created by the
 * CollectInsertSizeMetrics program and usually written to a file with the extension
 * ".insert_size_metrics".  In addition the insert size distribution is plotted to
 * a file with the extension ".insert_size_Histogram.pdf".
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class InsertSizeMetrics extends MultilevelMetrics {

    /** The MEDIAN insert size of all paired end reads where both ends mapped to the same chromosome. */
    public double MEDIAN_INSERT_SIZE;

    /** The MODE insert size of all paired end reads where both ends mapped to the same chromosome. */
    public double MODE_INSERT_SIZE;

    /**
     * The median absolute deviation of the distribution.  If the distribution is essentially normal then
     * the standard deviation can be estimated as ~1.4826 * MAD.
     */
    public double MEDIAN_ABSOLUTE_DEVIATION;

    /** The minimum measured insert size.  This is usually 1 and not very useful as it is likely artifactual. */
    public int MIN_INSERT_SIZE;
    /**
     * The maximum measure insert size by alignment. This is usually very high representing either an artifact
     * or possibly the presence of a structural re-arrangement.
     */
    public int MAX_INSERT_SIZE;
    /**
     * The mean insert size of the "core" of the distribution. Artefactual outliers in the distribution often
     * cause calculation of nonsensical mean and stdev values.  To avoid this the distribution is first trimmed
     * to a "core" distribution of +/- N median absolute deviations around the median insert size. By default
     * N=10, but this is configurable.
     */
    public double MEAN_INSERT_SIZE;
    /** Standard deviation of insert sizes over the "core" of the distribution. */
    public double STANDARD_DEVIATION;
    /** The total number of read pairs that were examined in the entire distribution. */
    public long READ_PAIRS;
    /** The pair orientation of the reads in this data category. */
    public PairOrientation PAIR_ORIENTATION;

    /** The "width" of the bins, centered around the median, that encompass 10% of all read pairs. */
    public int WIDTH_OF_10_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 20% of all read pairs. */
    public int WIDTH_OF_20_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 30% of all read pairs. */
    public int WIDTH_OF_30_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 40% of all read pairs. */
    public int WIDTH_OF_40_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 50% of all read pairs. */
    public int WIDTH_OF_50_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 60% of all read pairs. */
    public int WIDTH_OF_60_PERCENT;
    /**
     * The "width" of the bins, centered around the median, that encompass 70% of all read pairs.
     * This metric divided by 2 should approximate the standard deviation when the insert size
     * distribution is a normal distribution.
     */
    public int WIDTH_OF_70_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 80% of all read pairs. */
    public int WIDTH_OF_80_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 90% of all read pairs. */
    public int WIDTH_OF_90_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 95% of all read pairs. */
    public int WIDTH_OF_95_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 100% of all read pairs. */
    public int WIDTH_OF_99_PERCENT;
}
