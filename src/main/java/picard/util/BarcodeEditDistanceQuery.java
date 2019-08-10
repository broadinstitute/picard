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

/**
 * A class for finding the distance between multiple (matched) barcodes and multiple barcode reads.
 */
public class BarcodeEditDistanceQuery {

    /** list of barcodes (e.g. for dual-indexed barcodes) */
    public final byte[][] barcodeBytes;

    /** list of barcode reads */
    public final byte[][] readSubsequence;

    /** list of quality scores for reads */
    public final byte[][] qualities;

    /** minimal base quality to condiser informative */
    public final int minimumBaseQuality;

    /** maximal edit distance between reads and barcodes to be considered */
    public final int maximalInterestingDistance;

    public BarcodeEditDistanceQuery(final byte[][] barcodeBytes,
                                    final byte[][] readSubsequence,
                                    final byte[][] qualities,
                                    final int minimumBaseQuality,
                                    final int maximalInterestingDistance) {
        this.barcodeBytes = barcodeBytes;
        this.readSubsequence = readSubsequence;
        this.qualities = qualities;
        this.minimumBaseQuality = minimumBaseQuality;
        this.maximalInterestingDistance = maximalInterestingDistance;
    }

    public SingleBarcodeDistanceMetric getSingleBarcodeDistanceQuery(final int index, final int previousMismatchCount) {
        return new SingleBarcodeDistanceMetric(barcodeBytes[index],
                readSubsequence[index],
                qualities == null ? null : qualities[index],
                minimumBaseQuality,
                maximalInterestingDistance - previousMismatchCount);
    }
}
