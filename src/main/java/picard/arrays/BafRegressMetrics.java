package picard.arrays;

/*
 * The MIT License
 *
 * Copyright (c) 2020 The Broad Institute
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

import htsjdk.samtools.metrics.MetricBase;

public class BafRegressMetrics extends MetricBase {
    /** The sample name */
    public String SAMPLE;

    /** The estimate of contamination from the model (on the 0.0-1.0 scale) */
    public double ESTIMATE;

    /** The standard error of the estimate */
    public double STDERR;

    /** The test statistic for the estimate */
    public double TVAL;

    /** The p-value of the estimate */
    public double PVAL;

    /** The log p-value of the estimate */
    public double LOG10_PVAL;

    /** The call rate of the sample (number of non-missing genotypes) */
    public double CALL_RATE;

    /** The number of homozygous genotypes used to fit the model */
    public int NHOM;
}

