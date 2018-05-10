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

package picard.analysis;

import htsjdk.samtools.metrics.MetricBase;

/**
 * TheoreticalSensitivityMetrics, are metrics calculated from TheoreticalSensitivity and parameters used in
 * the calculation.
 *
 * These metrics are intended to estimate the achievable sensitivity (for a particular log odds threshold) of
 * somatic calls for a particular depth and base quality score distribution.
 *
 * @author Mark Fleharty
 */
public class TheoreticalSensitivityMetrics extends MetricBase {
    /** The allele fraction which theoretical sensitivity is calculated for. */
    public double ALLELE_FRACTION;
    /** Estimation of sensitivity at a particular allele fraction. */
    public double THEORETICAL_SENSITIVITY;
    /** Phred-scaled value of 1-THEORETICAL_SENSITIVITY. */
    public int THEORETICAL_SENSITIVITY_Q;
    /** Log-likelihood ratio is used as a threshold to distinguish a positive site with a given allele fraction from HOM_REF. */
    public double LOG_ODDS_THRESHOLD;
    /** Number of samples drawn at each depth in the depth distribution.  Larger values allow for increased precision at the cost of compute time. */
    public int SAMPLE_SIZE;
}
