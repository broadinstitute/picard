/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.util;

/**
 * a collection of functions and classes for various common calculations
 */
public class MathUtils {

    /**
     * a class for calculating moving statistics - this class returns the
     * mean, variance, and std dev after accumulating any number of records.
     * taken from http://www.johndcook.com/standard_deviation.html
     */
    public static class RunningStat {
        private double oldMean, newMean, oldStdDev, newStdDev;
        private long recordCount = 0;

        /**
         * @param x the value to add to the running mean and variance
         */
        public void push(double x) {
            recordCount++;
            // See Knuth TAOCP vol 2, 3rd edition, page 232
            if (recordCount == 1) {
                oldMean = newMean = x;
                oldStdDev = 0.0;
            } else {
                newMean = oldMean + (x - oldMean) / recordCount;
                newStdDev = oldStdDev + (x - oldMean) * (x - newMean);

                // set up for next iteration
                oldMean = newMean;
                oldStdDev = newStdDev;
            }
        }

        public void clear() { recordCount = 0; }
        public final long numDataValues() { return recordCount; }
        public final double mean() { return (recordCount > 0) ? newMean : 0.0; }
        public double variance() { return ((recordCount > 1) ? newStdDev / (recordCount - 1) : 0.0); }
        public double standardDeviation() { return Math.sqrt(variance()); }
    }

}
