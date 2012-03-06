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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;
import net.sf.picard.util.CollectionUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * Parse barcode file and return the matched barcode (not the actual sequence from the read).
 * If read does not match the returned BarcodeData will return null when getBarcode() is called.
 *
 * Barcode.txt file Format (consists of tab delimited columns, 1 record per row)
 * sequence_read    Matched(Y/N)    BarcodeSequenceMatched
 *
 * sequence read          - the actual bases at barcode position
 * Matched(y/n)           - Y or N indicating if there was a barcode match
 * BarcodeSequenceMatched - matched barcode sequence (empty if read did not match one of the barcodes).
 * @author jburke@broadinstitute.org
 */
class BarcodeParser implements IlluminaParser<BarcodeData> {
    private static final int Y_OR_N_COLUMN = 1;
    private static final int BARCODE_COLUMN = 2;
    private final IlluminaTextIterator textIterator;
    private final IlluminaFileMap tilesToFiles;
    
    private static final Set<IlluminaDataType> SupportedTypes = Collections.unmodifiableSet(CollectionUtil.makeSet(IlluminaDataType.Barcodes));

    public BarcodeParser(final int lane, final IlluminaFileMap tilesToFiles) {
        this.textIterator = new IlluminaTextIterator(lane, tilesToFiles, false);
        this.tilesToFiles = tilesToFiles;
    }

    @Override
    public void seekToTile(int oneBasedTileNumber) {
        textIterator.seekToTile(oneBasedTileNumber);
    }

    @Override
    public BarcodeData next() {
        final String [] fields = textIterator.next();
        final String barcode;
        if (fields[Y_OR_N_COLUMN].equals("Y")) {
            barcode = fields[BARCODE_COLUMN];
        } else {
            barcode = null;
        }

        return new BarcodeData() {
            @Override
            public String getBarcode() {
                return barcode;
            }
        };
    }

    @Override
    public boolean hasNext() {
        return textIterator.hasNext();
    }

    @Override
    public void verifyData(List<Integer> tiles, final int [] cycles) {
        if(tiles == null) {
            tiles = new ArrayList<Integer>();
            tiles.addAll(tilesToFiles.keySet());
        }
        
        for(final Integer tile : tiles) {
            final File barcodeFile = tilesToFiles.get(tile);
            if(barcodeFile == null) {
                throw new PicardException("Missing barcode file for tile(" + tile + ")");
            }

            if(!barcodeFile.exists()) {
                throw new PicardException("Barcode file (" + barcodeFile.getAbsolutePath() +" does not exist for tile(" + tile + ")");
            }
        }
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return SupportedTypes;
    }

    public void remove() {
        throw new UnsupportedOperationException("Remove is not supported by " + BarcodeParser.class.getName());
    }
}
