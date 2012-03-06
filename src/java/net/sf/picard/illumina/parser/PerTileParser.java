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
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;

/** Abstract base class for Parsers that open a single tile file at a time and iterate through them. */
public abstract class PerTileParser<ILLUMINA_DATA extends IlluminaData> implements IlluminaParser<ILLUMINA_DATA>  {
    private final IlluminaFileMap tileToFiles;
    private CloseableIterator<ILLUMINA_DATA> currentIterator;
    private Integer nextTile;
    private Integer currentTile;

    /** Factory method for the iterator of each tile */
    protected abstract CloseableIterator<ILLUMINA_DATA> makeTileIterator(final File nextTileFile);

    public PerTileParser(final IlluminaFileMap tilesToFiles) {
        this.tileToFiles = tilesToFiles;
        this.nextTile = tilesToFiles.firstKey();
        this.currentTile = null;
    }

    public PerTileParser(final IlluminaFileMap tilesToFiles, final int nextTile) {
        this.tileToFiles = tilesToFiles;
        this.currentTile = null;
        this.nextTile = nextTile;

        if(!tilesToFiles.containsKey(nextTile)) {
            throw new IllegalArgumentException("NextTile (" + nextTile + ") is not contained by tilesToFiles (" + StringUtil.join(",", new ArrayList<Integer>(tilesToFiles.keySet())));
        }
    }

    /**
     * Return the tile of the NEXT ILLUMINA_DATA object to be returned by the method next.  This might force us to advance to the
     * next file (as it will contains the data for the next) tile/ILLUMINA_DATA object.
     * @return tile number for the next ILLUMINA_DATA object to be returned
     */
    public int getTileOfNextCluster() {
        maybeAdvance();
        return currentTile;
    }

    private void advanceTile() {
        if(nextTile == null){
            throw new NoSuchElementException("No more tiles to advance!");
        }

        if(currentIterator != null) {
            currentIterator.close();
        }

        currentIterator = makeTileIterator(tileToFiles.get(nextTile));
        currentTile = nextTile;
        nextTile = tileToFiles.higherKey(nextTile);
    }

    public void seekToTile(int oneBasedTileNumber) {
        nextTile = oneBasedTileNumber;

        if(!tileToFiles.containsKey(oneBasedTileNumber)) {
            throw new PicardException("PerTileParser does not contain key(" + oneBasedTileNumber +") keys available (" + StringUtil.join(",", new ArrayList<Integer>(tileToFiles.keySet())) + ")");
        }

        if(currentIterator != null) {
            currentIterator.close();
        }
        currentIterator = null;
    }

    public void maybeAdvance() {
        if(!hasNext()) {
            throw new NoSuchElementException();
        }

        if(currentIterator == null || !currentIterator.hasNext()) {
            advanceTile();
        }
    }

    public ILLUMINA_DATA next() {
        maybeAdvance();

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
            throw new PicardException("Missing tiles in PerTileParser expected(" + StringUtil.join(",", tiles) + ") but found (" + StringUtil.join(",", mapTiles) + ")");
        }

        if(!tiles.containsAll(mapTiles)) {
            throw new PicardException("Extra tiles where found in PerTileParser  expected(" + StringUtil.join(",", tiles) + ") but found (" + StringUtil.join(",", mapTiles) + ")");
        }
    }
}
