/*
  * The MIT License
  *
  * Copyright (c) 2015 The Broad Institute
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

package picard.sam.markduplicates.util;

public class ReadEndsForMarkDuplicatesWithBarcodes extends ReadEndsForMarkDuplicates {
    public int barcode = 0; // primary barcode for this read (and pair)
    public int readOneBarcode = 0; // read one barcode, 0 if not present
    public int readTwoBarcode = 0; // read two barcode, 0 if not present or not paired

    public ReadEndsForMarkDuplicatesWithBarcodes() { }

    public ReadEndsForMarkDuplicatesWithBarcodes(final ReadEndsForMarkDuplicates read) {
        super(read);
    }

    public ReadEndsForMarkDuplicatesWithBarcodes(final ReadEndsForMarkDuplicatesWithBarcodes read) {
        super(read);
        barcode = read.barcode;
        readOneBarcode = read.readOneBarcode;
        readTwoBarcode = read.readTwoBarcode;
    }

    @Override
    public ReadEndsForMarkDuplicatesWithBarcodes clone() {
        return new ReadEndsForMarkDuplicatesWithBarcodes(this);
    }

    public static int getSizeOf() {
        return ReadEndsForMarkDuplicates.getSizeOf() + (3 * 4);
    }
}
