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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;
import net.sf.picard.illumina.parser.readers.BclIndexReader;
import net.sf.picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import net.sf.picard.illumina.parser.readers.BclReader;
import net.sf.samtools.util.CloseableIterator;

import java.io.File;
import java.util.NoSuchElementException;

/**
 * Parse .bcl.bgzf files that contain multiple tiles in a single file.  This requires an index file that tells
 * the bgzf virtual file offset of the start of each tile in the block-compressed bcl file.
 */
public class MultiTileBclParser extends BclParser {
    private final TileIndex tileIndex;
    public MultiTileBclParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles,
                              final OutputMapping outputMapping, final boolean applyEamssFilter,
                              final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                              final TileIndex tileIndex) {
        super(directory, lane, tilesToCycleFiles, outputMapping, applyEamssFilter, bclQualityEvaluationStrategy);
        this.tileIndex = tileIndex;
        super.initialize();
    }

    /**
     * Defer initialization until after this class is fully constructed.  This is necessary because superclass
     * ctor calls makeReader, which is overridden below, and the override requires that this.tileIndex is initialized,
     * and that doesn't happen until after superclass has been constructed.
     */
    @Override
    protected void initialize() {
    }

    @Override
    protected CloseableIterator<BclReader.BclValue> makeReader(final File file, final int cycle, final int tileNumber) {
        final TileIndex.TileIndexRecord tileIndexRecord = tileIndex.findTile(tileNumber);

        final BclIndexReader bclIndexReader = new BclIndexReader(file);
        if (tileIndex.getNumTiles() != bclIndexReader.getNumTiles()) {
            throw new PicardException(String.format("%s.getNumTiles(%d) != %s.getNumTiles(%d)",
                    tileIndex.getFile().getAbsolutePath(), tileIndex.getNumTiles(), bclIndexReader.getBciFile().getAbsolutePath(), bclIndexReader.getNumTiles()));
        }

        final BclReader bclReader = BclReader.makeSeekable(file, bclQualityEvaluationStrategy);
        bclReader.seek(bclIndexReader.get(tileIndexRecord.zeroBasedTileNumber));

        return new CountLimitedIterator(bclReader, tileIndexRecord.numClustersInTile);
    }

    /**
     * An iterator wrapper that stops when it has return a pre-determined number of records even if the underlying
     * iterator still had more records.
     */
    static class CountLimitedIterator implements CloseableIterator<BclReader.BclValue> {
        private final CloseableIterator<BclReader.BclValue> underlyingIterator;
        private final int recordLimit;
        private int numRecordsRead = 0;

        CountLimitedIterator(final CloseableIterator<BclReader.BclValue> underlyingIterator, final int recordLimit) {
            this.underlyingIterator = underlyingIterator;
            this.recordLimit = recordLimit;
        }

        @Override
        public void close() {
            underlyingIterator.close();
        }

        @Override
        public boolean hasNext() {
            return numRecordsRead < recordLimit && underlyingIterator.hasNext();
        }

        @Override
        public BclReader.BclValue next() {
            if (!hasNext()) throw new NoSuchElementException();
            ++numRecordsRead;
            return underlyingIterator.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
}
