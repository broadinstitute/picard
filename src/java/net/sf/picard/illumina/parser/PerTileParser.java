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
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

/** Abstract base class for Parsers that open a single tile file at a time and iterate through them. */
public abstract class PerTileParser<ILLUMINA_DATA extends IlluminaData> implements IlluminaParser<ILLUMINA_DATA>  {
    private final IlluminaFileMap tileToFiles;
    private Integer nextTile;

    private CloseableIterator<ILLUMINA_DATA> currentIterator;

    /** Factory method for the iterator of each tile */
    protected abstract CloseableIterator<ILLUMINA_DATA> makeTileIterator(final File iterator);

    public PerTileParser(final IlluminaFileMap tilesToFiles) {
        this.tileToFiles = tilesToFiles;
        this.nextTile = tilesToFiles.firstKey();
    }

    public PerTileParser(final IlluminaFileMap tilesToFiles, final int nextTile) {
        this.tileToFiles = tilesToFiles;
        this.nextTile = nextTile;

        if(!tilesToFiles.containsKey(nextTile)) {
            String tilesStr = "";
            boolean first = true;
            for(final Integer tile : tilesToFiles.keySet()) {
                if(!first) {
                    tilesStr += ", ";
                }

                tilesStr += tile;
            }

            throw new IllegalArgumentException("NextTile (" + nextTile + ") is not contained by tilesToFiles (" + tilesStr + ")");
        }
    }

    private void advanceTile() {
        if(nextTile == null){
            throw new NoSuchElementException("No more tiles to advance!");
        }

        if(currentIterator != null) {
            currentIterator.close();
        }

        currentIterator = makeTileIterator(tileToFiles.get(nextTile));
        nextTile = tileToFiles.higherKey(nextTile);
    }

    public void seekToTile(int oneBasedTileNumber) {
        nextTile = oneBasedTileNumber;

        if(!tileToFiles.containsKey(oneBasedTileNumber)) {
            String keys = "";
            int curKey = 0;
            if(tileToFiles.size() > 0) {
                curKey = tileToFiles.firstKey();
                keys += curKey;
            }

            for(int i = 1; i < tileToFiles.size(); i++) {
                curKey = tileToFiles.higherKey(curKey);
                keys += ", " + curKey;
            }

            throw new PicardException("PerTileParser does not contain key(" + oneBasedTileNumber +") keys available (" + keys + ")");
        }

        if(currentIterator != null) {
            currentIterator.close();
        }
        currentIterator = null;
    }

    public ILLUMINA_DATA next() {
        if(!hasNext()) {
            throw new NoSuchElementException();
        }

        if(currentIterator == null || !currentIterator.hasNext()) {
            advanceTile();
        }

        return currentIterator.next();
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public boolean hasNext() {
        return nextTile != null || (currentIterator != null && currentIterator.hasNext());
    }

    public void close() {
        if(currentIterator != null) {
            currentIterator.close();
        }
    }

    public void verifyData(List<Integer> tiles, final int [] cycles) {
        final List<Integer> mapTiles = new ArrayList<Integer>(this.tileToFiles.keySet());
        if(!mapTiles.containsAll(tiles)) {
            throw new PicardException("Missing tiles in PerTileParser expected(" + tilesToString(tiles) + ") but found (" + tilesToString(mapTiles) + ")");
        }

        if(!tiles.containsAll(mapTiles)) {
            throw new PicardException("Extra tiles where found in PerTileParser  expected(" + tilesToString(tiles) + ") but found (" + tilesToString(mapTiles) + ")");
        }
    }

    private static String tilesToString(final List<Integer> tiles) {
        String result = "";
        if(tiles.size() > 0) {
            result += tiles.get(0);
        }

        for(int i = 1; i < tiles.size(); i++) {
            result += ", " + tiles.get(i);
        }

        return result;
    }
}
