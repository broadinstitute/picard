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
package net.sf.picard.illumina.parser;

/**
 * Filter an IlluminaDataProvider so that only reads that match the given barcode are returned.
 *
 * @author alecw@broadinstitute.org
 */
public class BarcodeFilter implements AbstractIlluminaDataProvider {

    private final AbstractIlluminaDataProvider dataProvider;
    private final String barcode;
    private IlluminaReadData read;
    private final boolean includeUnmatched;

    /**
     * Create iterator that only returns reads that match the given barcode.
     * @param barcode Which barcode to match.
     * @param dataProvider Source of reads.
     * @param includeUnmatched Whether to include reads that match no barcode
     */
    public BarcodeFilter(final String barcode, final AbstractIlluminaDataProvider dataProvider, final boolean includeUnmatched) {
        if (barcode == null && !includeUnmatched) {
            throw new IllegalArgumentException("Barcode must be non-null when includeUnmatched is false");
        }
        this.barcode = barcode;
        this.dataProvider = dataProvider;
        this.includeUnmatched = includeUnmatched;
        advance();
    }

    /**
     * Create iterator that only returns reads that match the given barcode.
     * @param barcode Which barcode to match.
     * @param dataProvider Source of reads.
     */
    public BarcodeFilter(final String barcode, final AbstractIlluminaDataProvider dataProvider) {
        this(barcode, dataProvider, false);
    }

    private void advance() {
        while (dataProvider.hasNext()) {
            read = dataProvider.next();
            if ((barcode != null && barcode.equals(read.getMatchedBarcode())) || 
                (includeUnmatched && read.getMatchedBarcode() == null)) {
                return;
            }
        }
        read = null;
    }

    @Override
    public boolean hasNext() {
        return read != null;
    }

    @Override
    public IlluminaReadData next() {
        final IlluminaReadData ret = read;
        advance();
        return ret;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Jump so that the next record returned will be from the specified tile.
     */
    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        dataProvider.seekToTile(oneBasedTileNumber);
        advance();
    }
}
