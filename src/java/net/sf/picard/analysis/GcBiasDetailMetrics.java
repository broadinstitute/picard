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

package net.sf.picard.analysis;

import net.sf.picard.metrics.MetricBase;

/**
 * Class that holds detailed metrics about reads that fall within windows of a certain
 * GC bin on the reference genome.
 *
 * @author Tim Fennell
 */
public class GcBiasDetailMetrics extends MetricBase {
    /** The G+C content of the reference sequence represented by this bin. Values are from 0% to 100% */
    public int GC;

    /** The number of windows on the reference genome that have this G+C content. */
    public int WINDOWS;

    /** The number of reads who's start position is at the start of a window of this GC. */
    public long READ_STARTS;

    /** The mean quality (determined via the error rate) of all bases of all reads that are assigned to windows of this GC. */
    public int MEAN_BASE_QUALITY;

    /**
     * The ration of "coverage" in this GC bin vs. the mean coverage of all GC bins. A number of
     * 1 represents mean coverage, a number less than one represents lower than mean coverage (e.g. 0.5
     * means half as much coverage as average) while a number greater than one represents higher than
     * mean coverage (e.g. 3.1 means this GC bin has 3.1 times more reads per window than average).
     */
    public double NORMALIZED_COVERAGE;

    /**
     * The radius of error bars in this bin based on the number of observations made. For example if
     * the normalized coverage is 0.75 and the error bar width is 0.1 then the error bars would be
     * drawn from 0.65 to 0.85. 
     */
    public double ERROR_BAR_WIDTH;

}
