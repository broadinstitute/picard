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

import htsjdk.samtools.util.CloseableIterator;
import picard.illumina.parser.readers.BclQualityEvaluationStrategy;
import picard.illumina.parser.readers.BclReader;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Parse .bcl.bgzf files that contain multiple tiles in a single file.  This requires an index file that tells
 * the bgzf virtual file offset of the start of each tile in the block-compressed bcl file.
 */
public class MultiTileBclParser extends BclParser {
    private final TileIndex tileIndex;
    private MultiTileBclDataCycleFileParser cycleFileParser = null;
    public MultiTileBclParser(final File directory, final int lane, final CycleIlluminaFileMap tilesToCycleFiles,
                              final OutputMapping outputMapping, final boolean applyEamssFilter,
                              final BclQualityEvaluationStrategy bclQualityEvaluationStrategy,
                              final TileIndex tileIndex) {
        super(directory, lane, tilesToCycleFiles, outputMapping, applyEamssFilter, bclQualityEvaluationStrategy);
        this.tileIndex = tileIndex;
        this.initialize();
    }

    @Override
    public void initialize(){
        if(tileIndex != null){
            seekToTile(currentTile);
        }
    }

    private CountLimitedIterator makeReader(final List<File> files) {
        if(tileIndex != null) {
            final BclReader bclReader = BclReader.makeSeekable(files, bclQualityEvaluationStrategy, outputMapping.getOutputReadLengths());
            final int numClustersInTile = bclReader.seek(files, tileIndex, currentTile);
            return new CountLimitedIterator(bclReader, numClustersInTile);
        }
        else{
            return null;
        }
    }

    @Override
    protected CycleFilesParser<BclData> makeCycleFileParser(final List<File> files) {
        if (cycleFileParser == null) {
            cycleFileParser = new MultiTileBclDataCycleFileParser(files, currentTile);
        } else {
            final int numClustersInTile = cycleFileParser.getReader().seek(files, tileIndex, currentTile);
            cycleFileParser.setCurrentTile(currentTile);
            cycleFileParser.resetClusterLimit(numClustersInTile);
        }
        return cycleFileParser;
    }

    /**
     * An iterator wrapper that stops when it has return a pre-determined number of records even if the underlying
     * iterator still had more records.
     */
    static class CountLimitedIterator implements CloseableIterator<BclData> {
        public BclReader getUnderlyingIterator() {
            return underlyingIterator;
        }

        private final BclReader underlyingIterator;
        private int recordLimit;
        private int numRecordsRead = 0;

        CountLimitedIterator(final BclReader underlyingIterator, final int recordLimit) {
            this.underlyingIterator = underlyingIterator;
            this.recordLimit = recordLimit;
        }

        @Override
        public void close() {
            //underlyingIterator.close();
        }

        @Override
        public boolean hasNext() {
            return numRecordsRead < recordLimit && underlyingIterator.hasNext();
        }

        @Override
        public BclData next() {
            if (!hasNext()) throw new NoSuchElementException();
            ++numRecordsRead;
            return underlyingIterator.next();
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }


    private class MultiTileBclDataCycleFileParser implements CycleFilesParser<BclData> {
        final CountLimitedIterator reader;
        int currentTile;

        public MultiTileBclDataCycleFileParser(final List<File> files, final int currentTile) {
            this.currentTile = currentTile;
            reader = makeReader(files);
        }

        @Override
        public void close() {
            reader.close();
        }

        @Override
        public BclData next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            return reader.next();
        }

        @Override
        public boolean hasNext() {
            try {
                return reader.hasNext();
            } catch (final NullPointerException npe) {
                return false;
            }
        }

        public int getCurrentTile(){
            return currentTile;
        }

        public BclReader getReader() {
            return reader.getUnderlyingIterator();
        }

        public void resetClusterLimit(final int numClustersInTile) {
            reader.recordLimit = numClustersInTile;
            reader.numRecordsRead = 0;
        }

        public void setCurrentTile(final int currentTile) {
            this.currentTile = currentTile;
        }
    }
}
