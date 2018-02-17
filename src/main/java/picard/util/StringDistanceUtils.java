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
package picard.util;

import htsjdk.samtools.util.SequenceUtil;

import java.util.Arrays;

/**
 * Created by farjoun on 2/15/18.
 */
public class StringDistanceUtils {

    /**
     * Compare barcode sequence to bases from read
     *
     * @return how many bases did not match
     */
    public static int countMismatches(final byte[][] barcodeBytes, final byte[][] readSubsequence, final byte[][] qualities, final int minimumBaseQuality) {
        int numMismatches = 0;

        for (int j = 0; j < barcodeBytes.length; j++) {
            for (int i = 0; (i < barcodeBytes[j].length && readSubsequence[j].length > i); ++i) {
                if (SequenceUtil.isNoCall(readSubsequence[j][i])) {
                    continue;
                }
                if (!SequenceUtil.basesEqual(barcodeBytes[j][i], readSubsequence[j][i])) {
                    ++numMismatches;
                    continue;
                }
                if (qualities != null && qualities[j][i] < minimumBaseQuality) {
                    ++numMismatches;
                }
            }
        }

        return numMismatches;
    }

    /**
     * Compare barcode sequence to bases from read
     *
     * @return how many bases did not match
     */
    public static int countMismatchesWithIndelEvents(final byte[][] barcodeBytes, final byte[][] readSubsequence, final byte[][] qualities, final int minimumBaseQuality) {
        int numMismatches = 0;

        for (int j = 0; j < barcodeBytes.length; j++) {
            numMismatches += levenshteinDistanceForBarcodes(barcodeBytes[j], readSubsequence[j], qualities == null ? null : qualities[j], minimumBaseQuality);
        }

        return numMismatches;
    }

    public static int hammingDistance(final byte[] barcodeBases, final byte[] readBases, final byte[] readQualities, final int minimumBaseQuality) {
        int numMismatches = 0;
        for (int i = 0; i < barcodeBases.length && readBases.length > i; ++i) {
            if (SequenceUtil.isNoCall(readBases[i])) {
                continue;
            }
            if (!SequenceUtil.basesEqual(barcodeBases[i], readBases[i])) {
                ++numMismatches;
                continue;
            }
            if (readQualities != null && readQualities[i] < minimumBaseQuality) {
                ++numMismatches;
            }
        }
        return numMismatches;
    }

    static boolean anySmaller(final byte[] values, final int minValue) {
        for (int i = 0; i < values.length; i++) {
            if (values[i] < minValue) return true;
        }
        return false;
    }

    static byte[] maskSmaller(final byte[] bases, final byte[] qualities, final int minValue) {
        final byte[] maskedBases = Arrays.copyOf(bases, bases.length);

        for (int i = 0; i < qualities.length; i++) {
            if (qualities[i] < minValue) maskedBases[i] = '.';
        }
        return maskedBases;

    }

    public static int levenshteinDistanceForBarcodes(final byte[] barcodeBases, final byte[] readBases, final byte[] readQualities, final int minimumBaseQuality) {
        final byte[] maskedReadBases;
        if (readQualities == null || !anySmaller(readQualities, minimumBaseQuality)) {
            maskedReadBases = readBases;
        } else {
            maskedReadBases = maskSmaller(readBases, readQualities, minimumBaseQuality);
        }

        return levenshteinDistance(barcodeBases, readBases, 10);
    }

    /**
     * +++lifted from Commons Lang Text +++
     * <p>
     * modified to specific to comparing barcodes to reads. This means
     * ignoring insertions or deletions in the last position as it would amount to
     * double counting otherwise.
     * <p>
     * Also, in case of passing the threshold this algo has been modified to return
     * threshold + 1
     * <p>
     * This implementation only computes the distance if it's less than or
     * equal to the threshold value, returning -1 if it's greater. The
     * advantage is performance: unbounded distance is O(nm), but a bound of
     * k allows us to reduce it to O(km) time by only computing a diagonal
     * stripe of width 2k + 1 of the cost table. It is also possible to use
     * this to compute the unbounded Levenshtein distance by starting the
     * threshold at 1 and doubling each time until the distance is found;
     * this is O(dm), where d is the distance.
     * <p>
     * One subtlety comes from needing to ignore entries on the border of
     * our stripe eg. p[] = |#|#|#|* d[] = *|#|#|#| We must ignore the entry
     * to the leftRev of the leftmost member We must ignore the entry above the
     * rightmost member
     * *
     * As a concrete example, suppose s is of length 5, t is of length 7,
     * and our threshold is 1. In this case we're going to walk a stripe of
     * length 3. The matrix would look like so:
     * <p>
     * <pre>
     *    1 2 3 4 5
     * 1 |#|#| | | |
     * 2 |#|#|#| | |
     * 3 | |#|#|#| |
     * 4 | | |#|#|#|
     * 5 | | | |#|#|
     *
     *
     * This implementation decreases memory usage by using two
     * single-dimensional arrays and swapping them back and forth instead of
     * allocating an entire n by m matrix. This requires a few minor
     * changes, such as immediately returning when it's detected that the
     * stripe has run off the matrix and initially filling the arrays with
     * large values so that entries we don't compute are ignored.
     *
     * See Algorithms on Strings, Trees and Sequences by Dan Gusfield for
     * some discussion.
     */

    public static int levenshteinDistance(final byte[] barcode, final byte[] read, final int threshold) {
        int n = barcode.length; // length of left
        int m = read.length; // length of right

        if (n != m) {
            throw new IllegalArgumentException("This version of levenshteinDistance is speficially made for comparing strings " +
                    "of equal length. found " + n + " and " + m + ".");
        }
        // if one string is empty, the edit distance is necessarily the length
        // of the other
        if (n == 0) {
            return 0;
        }

        // it's easier to ignore indels in the begining than in the end...so we copy and reverse the arrays
        final byte[] barcodeRev = Arrays.copyOf(barcode, barcode.length);
        final byte[] readRev = Arrays.copyOf(read, read.length);

        SequenceUtil.reverse(barcodeRev, 0, barcodeRev.length);
        SequenceUtil.reverse(readRev, 0, readRev.length);

        int[] previousCost = new int[n + 1]; // 'previous' cost array, horizontally
        int[] cost = new int[n + 1]; // cost array, horizontally
        int[] tempD; // placeholder to assist in swapping previousCost and cost

        // fill in starting table values
        final int boundary = Math.min(n, threshold) + 1;
        for (int i = 0; i < boundary; i++) {
            previousCost[i] = 0;// indels in the beginning do not cost
        }
        // these fills ensure that the value above the rightmost entry of our
        // stripe will be ignored in following loop iterations
        Arrays.fill(previousCost, boundary, previousCost.length, Integer.MAX_VALUE);
        Arrays.fill(cost, Integer.MAX_VALUE);

        // iterates through t
        for (int j = 1; j <= m; j++) {
            final byte readJ = readRev[j - 1]; // jth character of readRev
            cost[0] = 0;

            // compute stripe indices, constrain to array size
            final int min = Math.max(1, j - threshold);
            final int max = j > Integer.MAX_VALUE - threshold ? n : Math.min(
                    n, j + threshold);

            // ignore entry barcodeRev of leftmost ?????
            if (min > 1) {
                cost[min - 1] = Integer.MAX_VALUE-10;
            }

            // iterates through [min, max] in s
            for (int i = min; i <= max; i++) {
                if (barcodeRev[i - 1] == readJ || SequenceUtil.isNoCall(readJ) ) {
                    // diagonally left and up
                    cost[i] = previousCost[i - 1];
                } else {
                    final int snpCost = previousCost[i - 1];
                    final int delCost = cost[i - 1];
                    final int insCost = previousCost[i];

                    cost[i] = 1 + Math.min(Math.min(snpCost, delCost), insCost);
                }
            }

            // copy current distance counts to 'previous row' distance counts
            tempD = previousCost;
            previousCost = cost;
            cost = tempD;
        }

        // if previousCost[n] is greater than the threshold, there's no guarantee on it
        // being the correct
        // distance
        if (previousCost[n] <= threshold) {
            return previousCost[n];
        }
        return threshold + 1;
    }
}
