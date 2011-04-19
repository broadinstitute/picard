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
package net.sf.picard.illumina;

import net.sf.picard.PicardException;

/**
 * @author alecw@broadinstitute.org
 */
public class BarcodeUtil {
    /**
     * Where on read the barcode is.
     */
    public static class BarcodePosition {
        public final boolean barcodeIsInSecondRead;

        /**
         * 0-based offset into the read where the barcode starts
         */
        public final int barcodeOffset;

        public final int barcodeLength;

        /**
         * Original 1-based cycle number of barcode start
         */
        public final int barcodeCycle;

        public BarcodePosition(final int barcodeOffset, final boolean barcodeIsInSecondRead, final int barcodeLength,
                               final int barcodeCycle) {
            this.barcodeOffset = barcodeOffset;
            this.barcodeIsInSecondRead = barcodeIsInSecondRead;
            this.barcodeLength = barcodeLength;
            this.barcodeCycle = barcodeCycle;
        }
    }

    /**
     * Determine which end contains barcode, and offset in that end where the barcode starts
     * @param barcodeCycle 1-based cycle number where the barcode starts
     */
    public static BarcodePosition findBarcodeEndAndStart(final boolean isPairedEnd, final int firstEndLength, final int secondEndLength,
                                                         final int barcodeCycle, final int barcodeLength) {
        final int barcodeOffset;
        final boolean isSecondEnd;
        if (!isPairedEnd && barcodeCycle + barcodeLength - 1 > firstEndLength) {
            throw new PicardException("Barcode position it past the end of unpaired read.");
        }
        if (barcodeCycle == 1) {
            if (barcodeLength > firstEndLength) {
                throw new PicardException("First end read is not long enough for barcode.");
            }
            barcodeOffset = 0;
            isSecondEnd = false;
        } else if (barcodeCycle <= firstEndLength) {
            if (barcodeCycle + barcodeLength != firstEndLength + 1) {
                throw new PicardException("Barcode position is not at end of first end.");
            }
            barcodeOffset = barcodeCycle - 1;
            isSecondEnd = false;
        } else if (barcodeCycle == firstEndLength + 1) {
            if (barcodeLength > secondEndLength) {
                throw new PicardException("Second end read is not long enough for barcode.");
            }
            barcodeOffset = 0;
            isSecondEnd = true;
        } else {
            if (barcodeCycle + barcodeLength != firstEndLength + secondEndLength + 1) {
                throw new PicardException("Barcode position is not at end of second end.");
            }
            barcodeOffset = barcodeCycle - firstEndLength - 1;
            isSecondEnd = true;
        }
        return new BarcodePosition(barcodeOffset, isSecondEnd, barcodeLength, barcodeCycle);
    }
}
