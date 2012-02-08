/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
import net.sf.picard.illumina.parser.readers.AbstractIlluminaPositionFileReader;
import net.sf.picard.illumina.parser.readers.ClocsFileReader;
import net.sf.picard.illumina.parser.readers.LocsFileReader;
import net.sf.picard.illumina.parser.readers.PosFileReader;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.Collections;
import java.util.Set;

import static net.sf.picard.util.CollectionUtil.makeSet;

/**
 * PosParser parses multiple files formatted as one of the three file formats that contain position information
 * only (pos, locs, and clocs).  This parser takes a map from tilesToFiles and a FileType enum value indicating
 * whether or not these are POS,LOCS, or CLOCS files.  The only client classes to this class should be IlluminaDataProvider
 * and test classes.  Check out AbstractIlluminaFileReader, PosFileReader, LocsFileReader, and ClocsFileReader for
 * more information on Position related illumina files.
 */
public class PosParser extends PerTileParser<PositionalData> {
    private static Set<IlluminaDataType> supportedTypes = Collections.unmodifiableSet(makeSet(IlluminaDataType.Position));

    /** The FileType of the files we are parsing */
    private final IlluminaFileUtil.SupportedIlluminaFormat fileType;

    public PosParser(final IlluminaFileMap tilesToFiles, final IlluminaFileUtil.SupportedIlluminaFormat fileType) {
        super(tilesToFiles);
        this.fileType = fileType;
    }

    public PosParser(final IlluminaFileMap tilesToFiles, final int startingTile, final IlluminaFileUtil.SupportedIlluminaFormat fileType) {
        super(tilesToFiles, startingTile);
        this.fileType = fileType;
    }

    /**
     * Make an CloseableIterator<PositionalData> based on the given file and fileType specified at construction.
     * This method wraps a reader in an iterator that converts it's output to the output format expected by
     * IlluminaDataProvider (PositionalData).
     * @param file A file for the current tile being parsed
     * @return An iterator over the PositionalData in that file.
     */
    @Override
    protected CloseableIterator<PositionalData> makeTileIterator(final File file) {

        final AbstractIlluminaPositionFileReader fileReader;
        switch(fileType){
            case Pos:
                fileReader = new PosFileReader(file);
                break;

            case Locs:
                fileReader = new LocsFileReader(file);
                break;

            case Clocs:
                fileReader = new ClocsFileReader(file);
                break;

            default:
                throw new PicardException("Unrecognized pos file type " + fileType.name());
        }

        return new CloseableIterator<PositionalData>() {
            private AbstractIlluminaPositionFileReader reader = fileReader;

            public void close() {
                reader.close();
            }

            public boolean hasNext() {
                return reader.hasNext();
            }

            public PositionalData next() {
                final AbstractIlluminaPositionFileReader.PositionInfo nextValue = reader.next();
                return new PositionalData() {
                    public int getXCoordinate() {
                        return nextValue.xQseqCoord;
                    }

                    public int getYCoordinate() {
                        return nextValue.yQseqCoord;
                    }

                    public int getLane() {
                        return nextValue.lane;
                    }

                    public int getTile() {
                        return nextValue.tile;
                    }
                };
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return supportedTypes;
    }

}
