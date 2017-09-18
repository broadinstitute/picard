/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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


import htsjdk.samtools.metrics.MetricsFile;
import picard.filter.CountingFilter;
import picard.filter.CountingPairedFilter;

/**
 * Interface for processing data and generate result for CollectWgsMetrics
 * @author Mariia_Zueva@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
public interface WgsMetricsProcessor {

    /**
     * Method processes the input data and accumulates result data
     */
    void processFile();

    /**
     * Adds result metric's data to input file
     *
     * @param file               MetricsFile for result of collector's work
     * @param includeBQHistogram include base quality histogram
     * @param dupeFilter         counting filter for duplicate reads
     * @param mapqFilter         counting filter for mapping quality
     * @param pairFilter         counting filter for reads without a mapped mate pair
     */
    void addToMetricsFile(final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> file,
            final boolean includeBQHistogram,
            final CountingFilter dupeFilter,
            final CountingFilter mapqFilter,
            final CountingPairedFilter pairFilter);
}
