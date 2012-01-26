/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

package net.sf.picard.util;

import java.math.BigDecimal;

import static java.lang.Math.pow;

/**
 * General math utilities
 *
 * @author Tim Fennell
 */
public class MathUtil {
    /** Calculated the mean of an array of doubles. */
    public static double mean(final double[] in, final int start, final int stop) {
        double total = 0;
        for (int i=start; i<stop; ++i) {
            total += in[i];
        }

        return total / (stop-start);
    }

    /** Calculated the standard deviation of an array of doubles. */
    public static double stddev(final double[] in, final int start, final int length) {
        return stddev(in, start, length, mean(in, start, length));
    }

    /** Calculated the standard deviation of an array of doubles. */
    public static double stddev(final double[] in, final int start, final int stop, final double mean) {
        double total = 0;
        for (int i=start; i<stop; ++i) {
            total += (in[i] * in[i]);
        }

        return Math.sqrt((total / (stop-start)) - (mean*mean));
    }

    public static int compare(final int v1, final int v2) {
        return (v1 < v2? -1: (v1 == v2? 0: 1));
    }

    /**
     * Obtains percentage of two Longs
     * @param numerator dividend
     * @param denominator divisor
     * @return numerator/(double)denominator if both are non-null and denominator != 0, else returns null.
     */
    public static Double percentageOrNull(final Long numerator, final Long denominator) {
        if (numerator != null && denominator != null && denominator != 0) {
            return numerator.doubleValue()/denominator.doubleValue();
        } else {
            return null;
        }
    }

    /** 
     * Round off the value to the specified precision. 
     */
    public static double round(final double num, final int precision) {
        BigDecimal bd = new BigDecimal(num);
        bd = bd.setScale(precision, BigDecimal.ROUND_HALF_UP);
        return  bd.doubleValue();
    }

    /** Returns the largest value stored in the array. */
    public static double max(final double[] nums) {
        double max = nums[0];
        for (int i=1; i<nums.length; ++i) {
            if (nums[i] > max) max = nums[i];
        }

        return max;
    }

    /** Returns the smallest value stored in the array. */
    public static double min(final double[] nums) {
        double min = nums[0];
        for (int i=1; i<nums.length; ++i) {
            if (nums[i] < min) min = nums[i];
        }

        return min;
    }

    /** "Promotes" an int[] into a double array with the same values (or as close as precision allows). */
    public static double[] promote(final int[] is) {
        final double[] ds = new double[is.length];
        for (int i=0; i<is.length; ++i) ds[i] = is[i];
        return ds;
    }

    /**
     * Takes a complete set of mutually exclusive log likelihoods and converts them to probabilities
     * that sum to 1 with as much fidelity as possible.
     */
    public static double[] logLikelihoodsToProbs(final double[] likelihoods) {
        // Note: bumping all the LLs so that the biggest is 300 ensures that we have the
        // widest range possible when unlogging them before one of them underflows. 10^300 is
        // near the maximum before you hit positive infinity.

        final double maxLikelihood = max(likelihoods);
        final double bump = 300 - maxLikelihood;

        final double[] tmp = new double[likelihoods.length];
        double total = 0;
        for (int i=0; i<likelihoods.length; ++i) {
            tmp[i] = pow(10,likelihoods[i]+bump);
            total += tmp[i];
        }

        for (int i=0; i<likelihoods.length; ++i) {
            tmp[i] /= total;
        }

        return tmp;
    }
}
