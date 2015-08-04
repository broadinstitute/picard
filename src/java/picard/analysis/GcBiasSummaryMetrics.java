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

import picard.metrics.MultilevelMetrics;

/**
 * High level metrics that capture how biased the coverage in a certain lane is.
 *
 * @author Tim Fennell
 */
public class GcBiasSummaryMetrics extends MultilevelMetrics {
    public String ACCUMULATION_LEVEL;

    /** The window size on the genome used to calculate the GC of the sequence. */
    public int WINDOW_SIZE;

    /** The total number of clusters that were seen in the gc bias calculation. */
    public int TOTAL_CLUSTERS;

    /** The total number of aligned reads used to compute the gc bias metrics. */
    public int ALIGNED_READS;

    /**
     * Illumina-style AT dropout metric.  Calculated by taking each GC bin independently and calculating
     * (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].
     */
    public double AT_DROPOUT;

    /**
     * Illumina-style GC dropout metric.  Calculated by taking each GC bin independently and calculating
     * (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].
     */
    public double GC_DROPOUT;

    /**
     * Normalized coverage over each quintile of GC content.
     */
    public double GC_NC_0_20;
    public double GC_NC_20_40;
    public double GC_NC_40_60;
    public double GC_NC_60_80;
    public double GC_NC_80_100;
}
