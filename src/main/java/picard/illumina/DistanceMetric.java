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
package picard.illumina;

import org.broadinstitute.barclay.argparser.CommandLineParser;
import picard.util.BarcodeEditDistanceQuery;
import picard.util.SingleBarcodeDistanceMetric;

public enum DistanceMetric implements CommandLineParser.ClpEnum {
    HAMMING("Hamming distance: The n-th base in the read is compared against the n-th base in the barcode. " +
            "Unequal bases and low quality bases are considered mismatches. No-call read-bases are not considered mismatches. ") {
        @Override
        protected int distance0(final SingleBarcodeDistanceMetric editDistance) {
            return editDistance.hammingDistance();
        }
    },

    LENIENT_HAMMING("Leniant Hamming distance: The n-th base in the read is compared against the n-th base in the barcode. " +
            "Unequal bases are considered mismatches. No-call read-bases, or those with low quality are not considered mismatches.") {
        @Override
        protected int distance0(final SingleBarcodeDistanceMetric editDistance) {
            return editDistance.lenientHammingDistance();
        }
    },

    FREE("FREE Metric: A Levenshtein-like metric that performs a simple Smith-Waterman with mismatch, gap open, " +
            "and gap extend costs all equal to 1. " +
            "Insertions or deletions at the ends of the read or barcode do not count toward the distance. " +
            "No-call read-bases, or those with low quality are not considered mismatches.") {
        @Override
        protected int distance0(final SingleBarcodeDistanceMetric editDistance) {
            return editDistance.freeDistance();
        }
    };

    final private String helpString;

    DistanceMetric(final String helpString) {
        this.helpString = helpString;
    }

    protected abstract int distance0(final SingleBarcodeDistanceMetric editDistance);

    public int distance(final BarcodeEditDistanceQuery editDistance) {
        int numMismatches = 0;

        for (int j = 0; j < editDistance.barcodeBytes.length; j++) {
            final SingleBarcodeDistanceMetric singleBarcodeDistanceMetric = editDistance.getSingleBarcodeDistanceQuery(j,numMismatches);
            numMismatches += distance0(singleBarcodeDistanceMetric);
            if (numMismatches > editDistance.maximalInterestingDistance) {
                return numMismatches;
            }
        }
        return numMismatches;
    }

    @Override
    public String getHelpDoc() {
        return this.helpString;
    }
}
