/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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
 * A class for finding the distance between a single barcode and a barcode-read (with base qualities)
 */
public class SingleBarcodeDistanceMetric {

    // The bases that the (comparison) barcode has
    private final byte[] barcodeBases;

    // The bases that the read has
    private final byte[] readBases;

    // the Read qualities
    private final byte[] readQualities;

    // The minimum base quality in the read that will be considered informative.
    private final int minimumBaseQuality;

    // a distance beyond which the user doesn't care how large the distance is
    private final int maximalInterestingDistance;
    private final byte[] maskedBases;

    public SingleBarcodeDistanceMetric(final byte[] barcodeBases,
                                       final byte[] readBases,
                                       final byte[] readQualities,
                                       final int minimumBaseQuality,
                                       final int maximalInterestingDistance) {
        this.barcodeBases = barcodeBases;
        this.readQualities = readQualities;
        this.readBases = readBases;
        this.maskedBases = maskIfAnySmaller(readBases,readQualities,minimumBaseQuality);
        this.minimumBaseQuality = minimumBaseQuality;
        this.maximalInterestingDistance = maximalInterestingDistance;
    }

    public int hammingDistance() {
        int numMismatches = 0;
        for (int i = 0; i < barcodeBases.length && i < readBases.length && numMismatches <= maximalInterestingDistance; ++i) {

            if (SequenceUtil.isNoCall(readBases[i])) {
                continue;
            }

            if (!SequenceUtil.basesEqual(barcodeBases[i], readBases[i])) {
                ++numMismatches;
                continue;
            }

            // bases are equal but if quality is low we still penalize.
            // TODO: is it actually useful to penalize barcodes in this way?
            if (readQualities != null && readQualities[i] < minimumBaseQuality) {
                ++numMismatches;
            }
        }
        return numMismatches;
    }


    @Deprecated //(due to typo in original name. 1/9/2020)
    public int leniantHammingDistance() {
        return lenientHammingDistance();
    }
        /**
         * Similar to Hamming distance but this version doesn't penalize matching bases with low quality for the read.

         * @return the edit distance between the barcode(s) and the read(s)
         */
    public int lenientHammingDistance() {
        int numMismatches = 0;
        for (int i = 0; i < barcodeBases.length && i < maskedBases.length && numMismatches <= maximalInterestingDistance; ++i) {

            if (SequenceUtil.isNoCall(maskedBases[i])) {
                continue;
            }

            if (!SequenceUtil.basesEqual(barcodeBases[i], maskedBases[i])) {
                ++numMismatches;
            }
        }
        return numMismatches;
    }

    private static boolean anySmaller(final byte[] values, final int minValue) {
        if (values == null) {
            return false;
        }
        for (final byte value : values) {
            if (value < minValue) {
                return true;
            }
        }
        return false;
    }

    static private byte[] maskIfAnySmaller(final byte [] readBases, final byte[] readQualities, final int minimumBaseQuality) {

        if (!anySmaller(readQualities,minimumBaseQuality)) {
            return readBases;
        }

        final byte[] maskedBases = Arrays.copyOf(readBases, readBases.length);

        for (int i = 0; i < readQualities.length; i++) {
            if (readQualities[i] < minimumBaseQuality) {
                maskedBases[i] = '.';
            }
        }
        return maskedBases;
    }

    /**
     * +++lifted from Commons Lang Text +++
     * <p>
     * Modified to specific to comparing barcodes to reads. This means ignoring
     * insertions or deletions in the last position as it would amount to
     * double counting otherwise. Based on https://www.pnas.org/content/115/27/E6217
     * <p>
     * Also, in case of passing the threshold this algorithm has been modified to return
     * maximalInterestingDistance + 1
     * <p>
     * This implementation only computes the distance if it is less than or
     * equal to the threshold value, returning maximalInterestingDistance + 1
     * otherwise.
     * The advantage is performance: unbounded distance is O(nm), but a bound of
     * k allows us to reduce it to O(km) time by only computing a diagonal
     * stripe of width 2k + 1 of the cost table. It is also possible to use
     * this to compute the unbounded Levenshtein distance by starting the
     * threshold at 1 and doubling each time until the distance is found;
     * this is O(dm), where d is the distance.
     * <p>
     * One subtlety comes from needing to ignore entries on the border of
     * our stripe eg. p[] = |#|#|#|* d[] = *|#|#|#| We must ignore the entry
     * to the leftRev of the leftmost member We must ignore the entry above the
     * rightmost member.
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
     * </pre>
     *
     * in addition, we terminate the calculation early and return threshold +1
     * if we detect that there is no way to get to the bottom right corner without
     * going over the thershold.
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

    public  int freeDistance() {
        final int EDIT_COST = 1;
        final int n = barcodeBases.length; // length of barcodeBases

        if (n != readBases.length) {
            throw new IllegalArgumentException("This version of freeDistance is specifically made for comparing strings " +
                    "of equal length. found " + n + " and " + readBases.length + ".");
        }

        if (n == 0) {
            return 0;
        }

        // it's easier to ignore indels in the beginning than in the end...so we copy and reverse the arrays
        final byte[] barcodeRev = Arrays.copyOf(barcodeBases, barcodeBases.length);
        final byte[] readRev = Arrays.copyOf(maskedBases, maskedBases.length);

        //reverseQualities reverses without complementing....which is what I want here.
        SequenceUtil.reverseQualities(barcodeRev);
        SequenceUtil.reverseQualities(readRev);

        int[] previousCost = new int[n + 1]; // 'previous' cost array, horizontally
        int[] cost = new int[n + 1]; // cost array, horizontally
        int[] tempD; // placeholder to assist in swapping previousCost and cost

        // fill in starting table values
        final int boundary = Math.min(n, maximalInterestingDistance) + 1;
        for (int i = 0; i < boundary; i++) {
            previousCost[i] = 0; // indels in the beginning do not cost
        }
        // these fills ensure that the value above the rightmost entry of our
        // stripe will be ignored in following loop iterations
        Arrays.fill(previousCost, boundary, previousCost.length, Integer.MAX_VALUE);
        Arrays.fill(cost, Integer.MAX_VALUE);

        // iterates through t
        for (int j = 1; j <= n; j++) {
            final byte readJ = readRev[j - 1]; // jth character of readRev
            cost[0] = 0;

            // compute stripe indices, constrain to array size
            final int min = Math.max(1, j - maximalInterestingDistance);
            final int max = j > Integer.MAX_VALUE - maximalInterestingDistance ? n : Math.min(
                    n, j + maximalInterestingDistance);

            // ignore leftmost entry of barcodeRev
            if (min > 1) {
                cost[min - 1] = Integer.MAX_VALUE - maximalInterestingDistance;
            }

            // iterates through [min, max] in barcodeRev
            int minCost = Integer.MAX_VALUE;
            for (int i = min; i <= max; i++) {
                if (barcodeRev[i - 1] == readJ || SequenceUtil.isNoCall(readJ)) {
                    // diagonally left and up
                    cost[i] = previousCost[i - 1];
                } else {
                    final int snpCost = previousCost[i - 1];
                    final int delCost = cost[i - 1];
                    final int insCost = previousCost[i];

                    cost[i] = Math.min(Math.min(snpCost, delCost), insCost) + EDIT_COST;
                }
                // be adding the distnace to the "center" we can estimate the lowest cost that this entry would
                // lead to assuming there are no more SNP errors (only indels)
                final int distToCenter = Math.abs(i - j);

                minCost = Math.min(minCost, cost[i] + distToCenter * EDIT_COST);
            }

            if (minCost > maximalInterestingDistance) {
                return maximalInterestingDistance + EDIT_COST;
            }
            // copy current distance counts to 'previous row' distance counts
            tempD = previousCost;
            previousCost = cost;

            cost = tempD;
        }

        // if previousCost[n] is greater than the threshold, there's no guarantee on it
        // being the correct distance. but it's definitely bigger than
        // maximalInterestingDistance
        if (previousCost[n] > maximalInterestingDistance) {
            return maximalInterestingDistance + EDIT_COST;
        }
        return previousCost[n];
    }
}
