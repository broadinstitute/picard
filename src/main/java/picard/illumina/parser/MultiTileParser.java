/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package picard.illumina.parser;

import htsjdk.samtools.util.PeekIterator;
import picard.PicardException;

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * Abstract class for files with fixed-length records for multiple tiles, e.g. .locs and .filter files.
 * @param <OUTPUT_RECORD> The kind of record to be returned (as opposed to the type of the record stored in the file).
 */
public abstract class MultiTileParser<OUTPUT_RECORD extends IlluminaData> implements IlluminaParser<OUTPUT_RECORD> {
    private final TileIndex tileIndex;
    private final Iterator<TileIndex.TileIndexRecord> tileIndexIterator;
    private final PeekIterator<Integer> requestedTilesIterator;
    private final Set<IlluminaDataType> supportedTypes;
    private int nextRecordIndex = 0;
    private int nextClusterInTile;
    private TileIndex.TileIndexRecord currentTile = null;

    /**
     * @param tileIndex Enables conversion from tile number to record number in this file.
     * @param requestedTiles Iterate over these tile numbers, which must be in ascending order.
     * @param supportedTypes The data types(s) that are provided by this file type, used to decide what file types to read.
     */
    public MultiTileParser(final TileIndex tileIndex,
                           final List<Integer> requestedTiles,
                           final Set<IlluminaDataType> supportedTypes) {
        this.tileIndex = tileIndex;
        this.tileIndexIterator = tileIndex.iterator();
        this.requestedTilesIterator = new PeekIterator<Integer>(requestedTiles.iterator());
        this.supportedTypes = supportedTypes;
    }

    @Override
    public void seekToTile(final int oneBasedTileNumber) {
        while (tileIndexIterator.hasNext()) {
            final TileIndex.TileIndexRecord next = tileIndexIterator.next();
            if (next.tile > oneBasedTileNumber) {
                throw new PicardException(
                        String.format("Cannot seek backwards: next tile %d > tile sought %d", next.tile, oneBasedTileNumber));
            } else if (next.tile == oneBasedTileNumber) {
                currentTile = next;
                break;
            }
        }
        if (nextRecordIndex > currentTile.indexOfFirstClusterInTile) {
            throw new PicardException(
                    String.format("Seem to be in wrong position %d > %d", nextRecordIndex, currentTile.indexOfFirstClusterInTile));
        }
        skipRecords(currentTile.indexOfFirstClusterInTile - nextRecordIndex);
        nextRecordIndex = currentTile.indexOfFirstClusterInTile;
        nextClusterInTile = 0;
    }

    @Override
    public OUTPUT_RECORD next() {
        if (!hasNext()) throw new NoSuchElementException();
        OUTPUT_RECORD ret = readNext();
        ++nextClusterInTile;
        ++nextRecordIndex;
        return ret;
    }

    @Override
    public boolean hasNext() {
        // Skip over any empty tiles
        while ((currentTile == null || nextClusterInTile >= currentTile.numClustersInTile) && requestedTilesIterator.hasNext()) {
            seekToTile(requestedTilesIterator.next());
        }
        return currentTile != null && nextClusterInTile < currentTile.numClustersInTile;
    }

    @Override
    public int getTileOfNextCluster() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        if (currentTile != null && nextClusterInTile < currentTile.numClustersInTile) return currentTile.tile;
        else return requestedTilesIterator.peek();
    }

    @Override
    public void verifyData(final List<Integer> tiles, final int[] cycles) {
        final List<String> tileErrors = tileIndex.verify(tiles);
        if (!tileErrors.isEmpty()) throw new PicardException(tileErrors.get(0));
        //No need to validate cycles until such time as this class is used for cycle-oriented data types
    }

    @Override
    public Set<IlluminaDataType> supportedTypes() {
        return supportedTypes;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    abstract OUTPUT_RECORD readNext();
    abstract void skipRecords(int numToSkip);
}
