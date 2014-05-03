/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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

package net.sf.picard.metrics;

/**
 *  PerRecordCollector - An interface for classes that collect data in order to generate one or more metrics.
 *      This process usually occurs in the following fashion:
 *      1. Loop through a data set (usually all records in a BAM file) and call collector.acceptRecord( data ),
 *         data in this step is usually added to metrics/histogram objects
 *      2. Call collector.finish() - perform any final calculations necessary after ALL records have been accepted
 *      3. addMetricsToFile is then used to add any metric(s) or histogram(s) to the given file
 *
 *      BEAN    - The Metric type we are generating
 *      HKEY    - The Key used in any histograms, use a Wildcard(?) type if there are no histograms
 *      ARGTYPE - Collectors are often used in groups of accumulation levels, in order to avoid recalculating
 *                any information needed by multiple collectors we allow different types of arguments that
 *                extend DefaultPerRecordCollectorArgs to accommodate any computed values
 */
public interface PerUnitMetricCollector<BEAN extends MetricBase, HKEY extends Comparable, ARGTYPE> {
    /**
     * Add a SAMRecord (with ReferenceSequence and Read Group info) to the metric(s) being calculated)
     * @param args Contains SAMRecord, SAMReadGroupRecord, ReferenceSequence of current record and any previously
     *             computed values that might be needed for this class
     */
    public void acceptRecord(final ARGTYPE args);

    /** When all records have been collected, compute any final values needed to finish constructing metrics/histogram */
    public void finish();

    /**
     * Any metrics collected will be added to the metric file provided.
     * @param file MetricsFile to which all metrics created by this collector should be added
     */
    public void addMetricsToFile(final MetricsFile<BEAN, HKEY> file);
}

