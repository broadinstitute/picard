/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble;

import net.sf.samtools.seekablestream.SeekableStream;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broad.tribble.util.ParsingUtils;
import net.sf.samtools.seekablestream.SeekableStreamFactory;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * A reader for text feature files  (i.e. not tabix files).   This includes tribble-indexed and non-indexed files.  If
 * index both iterate() and query() methods are supported.
 * <p/>
 * Note: Non-indexed files can be gzipped, but not bgzipped.
 *
 * @author Jim Robinson
 * @since 2/11/12
 */
public class TribbleIndexedFeatureReader<T extends Feature> extends AbstractFeatureReader {
    private Index index;

    /** is the path pointing to our source data a regular file? */
    private final boolean pathIsRegularFile;

    /** a potentially reusable seekable stream for queries over regular files */
    private SeekableStream seekableStream = null;

    /**
     * @param featurePath  - path to the feature file, can be a local file path, http url, or ftp url
     * @param codec        - codec to decode the features
     * @param requireIndex - true if the reader will be queries for specific ranges.  An index (idx) file must exist
     * @throws IOException
     */
    public TribbleIndexedFeatureReader(String featurePath, FeatureCodec codec, boolean requireIndex) throws IOException {

        super(featurePath, codec);

        if (requireIndex) {
            String indexFile = Tribble.indexFile(featurePath);
            if (ParsingUtils.resourceExists(indexFile)) {
                index = IndexFactory.loadIndex(indexFile);
            } else {
                // See if the index itself is gzipped
                indexFile = indexFile + ".gz";
                if (ParsingUtils.resourceExists(indexFile)) {
                    index = IndexFactory.loadIndex(indexFile);
                } else {
                    throw new TribbleException("An index is required, but none found.");
                }
            }
        }

        // does path point to a regular file?
        this.pathIsRegularFile = SeekableStreamFactory.isFilePath(path);

        readHeader();
    }

    /**
     * Get a seekable stream appropriate to read information from the current feature path
     *
     * This function ensures that if reuseStreamInQuery returns true then this function will only
     * ever return a single unique instance of SeekableStream for all calls given this instance of
     * TribbleIndexedFeatureReader.  If reuseStreamInQuery() returns false then the returned SeekableStream
     * will be newly opened each time, and should be closed after each use.
     *
     * @return a SeekableStream
     */
    private SeekableStream getSeekableStream() throws IOException {
        final SeekableStream result;
        if ( reuseStreamInQuery() ) {
            // if the stream points to an underlying file, only create the underlying seekable stream once
            if ( seekableStream == null ) seekableStream = SeekableStreamFactory.getStreamFor(path);
            result = seekableStream;
        } else {
            // we are not reusing the stream, so make a fresh copy each time we request it
            result = SeekableStreamFactory.getStreamFor(path);
        }

        return result;
    }

    /**
     * Are we attempting to reuse the underlying stream in query() calls?
     * @return true if
     */
    private boolean reuseStreamInQuery() {
        return pathIsRegularFile;
    }

    /**
     * @param featureFile - path to the feature file, can be a local file path, http url, or ftp url
     * @param codec       - codec to decode the features
     * @param index       - a tribble Index object
     * @throws IOException
     */
    public TribbleIndexedFeatureReader(String featureFile, FeatureCodec codec, Index index) throws IOException {
        this(featureFile, codec, false); // required to read the header
        this.index = index;
    }


    public void close() throws IOException {
        // close the seekable stream if that's necessary
        if ( seekableStream != null ) seekableStream.close();
    }

    /**
     * Return the sequence (chromosome/contig) names in this file, if known.
     *
     * @return list of strings of the contig names
     */
    public List<String> getSequenceNames() {
        return index == null ? new ArrayList<String>() : new ArrayList<String>(index.getSequenceNames());
    }


    /**
     * read the header from the file
     *
     * @throws IOException throws an IOException if we can't open the file
     */
    private void readHeader() throws IOException {
        InputStream is = null;
        PositionalBufferedStream pbs = null;
        try {
            is = ParsingUtils.openInputStream(path);
            if (path.endsWith("gz")) {
                // TODO -- warning I don't think this can work, the buffered input stream screws up position
                is = new GZIPInputStream(new BufferedInputStream(is));
            }
            pbs = new PositionalBufferedStream(is);
            header = codec.readHeader(pbs);
        } catch (Exception e) {
            throw new TribbleException.MalformedFeatureFile("Unable to parse header with error: " + e.getMessage(), path, e);
        } finally {
            if ( pbs != null ) pbs.close();
            else if (is != null) is.close();
        }
    }

    /**
     * Return an iterator to iterate over features overlapping the specified interval
     *
     * Note that TribbleIndexedFeatureReader only supports issuing and manipulating a single query
     * for each reader.  That is, the behavior of the following code is undefined:
     *
     * reader = new TribbleIndexedFeatureReader()
     * Iterator it1 = reader.query("x", 10, 20)
     * Iterator it2 = reader.query("x", 1000, 1010)
     *
     * As a consequence of this, the TribbleIndexedFeatureReader are also not thread-safe.
     *
     * @param chr contig
     * @param start start position
     * @param end end position
     * @return an iterator of records in this interval
     * @throws IOException
     */
    public CloseableTribbleIterator query(String chr, int start, int end) throws IOException {

        if (index == null) {
            throw new TribbleException("Index not found for: " + path);
        }

        if (index.containsChromosome(chr)) {
            List<Block> blocks = index.getBlocks(chr, start - 1, end);
            return new QueryIterator(chr, start, end, blocks);
        } else {
            return new EmptyIterator();
        }
    }


    /**
     * @return Return an iterator to iterate over the entire file
     * @throws IOException
     */
    public CloseableTribbleIterator iterator() throws IOException {
        return new WFIterator();
    }

    /**
     * Class to iterator over an entire file.
     *
     * @param <T>
     */
    class WFIterator<T extends Feature> implements CloseableTribbleIterator {
        private T currentRecord;
        private PositionalBufferedStream stream;

        /**
         * Constructor for iterating over the entire file (seekableStream).
         *
         * @throws IOException
         */
        public WFIterator() throws IOException {
            final InputStream inputStream = ParsingUtils.openInputStream(path);

            if (path.endsWith(".gz")) {
                // Gzipped -- we need to buffer the GZIPInputStream methods as this class makes read() calls,
                // and seekableStream does not support single byte reads
                InputStream is = new GZIPInputStream(new BufferedInputStream(inputStream, 512000));
                stream = new PositionalBufferedStream(is, 1000);  // Small buffer as this is buffered already.
            } else {
                stream = new PositionalBufferedStream(inputStream, 512000);
            }

            if ( header.skipHeaderBytes() ) stream.skip(header.getHeaderEnd());
            readNextRecord();
        }


        public boolean hasNext() {
            return currentRecord != null;
        }

        public T next() {
            T ret = currentRecord;
            try {
                readNextRecord();
            } catch (IOException e) {
                throw new RuntimeException("Unable to read the next record, the last record was at " +
                        ret.getChr() + ":" + ret.getStart() + "-" + ret.getEnd(), e);
            }
            return ret;
        }


        /**
         * Advance to the next record in the query interval.
         *
         * @throws IOException
         */
        private void readNextRecord() throws IOException {
            currentRecord = null;

            while (! stream.isDone() ) {
                Feature f = null;
                try {
                    f = codec.decode(stream);

                    if (f == null) {
                        continue;
                    }

                    currentRecord = (T) f;
                    return;

                } catch (TribbleException e) {
                    e.setSource(path);
                    throw e;
                } catch (NumberFormatException e) {
                    String error = "Error parsing line at byte position: " + stream.getPosition();
                    throw new TribbleException.MalformedFeatureFile(error, path, e);
                }
            }
        }

        public void remove() {
            throw new UnsupportedOperationException("Remove is not supported in Iterators");
        }

        public void close() {
            stream.close();
        }

        public WFIterator<T> iterator() {
            return this;
        }
    }

    /**
     * Iterator for a query interval
     *
     * @param <T>
     */
    class QueryIterator<T extends Feature> implements CloseableTribbleIterator {
        private String chr;
        private String chrAlias;
        int start;
        int end;
        private T currentRecord;
        private SeekableStream mySeekableStream;
        private PositionalBufferedStream stream;
        private Iterator<Block> blockIterator;

        public QueryIterator(String chr, int start, int end, List<Block> blocks) throws IOException {
            this.chr = chr;
            this.start = start;
            this.end = end;
            mySeekableStream = getSeekableStream();
            blockIterator = blocks.iterator();
            advanceBlock();
            readNextRecord();

            // The feature chromosome might not be the query chromosome, due to alias definitions.  We assume
            // the chromosome of the first record is correct and record it here.  This is not pretty.
            chrAlias = (currentRecord == null ? chr : currentRecord.getChr());

        }


        public boolean hasNext() {
            return currentRecord != null;
        }

        public T next() {
            T ret = currentRecord;
            try {
                readNextRecord();
            } catch (IOException e) {
                throw new RuntimeException("Unable to read the next record, the last record was at " +
                        ret.getChr() + ":" + ret.getStart() + "-" + ret.getEnd(), e);
            }
            return ret;
        }


        private void advanceBlock() throws IOException {
            while (blockIterator != null && blockIterator.hasNext()) {
                Block block = blockIterator.next();
                if (block.getSize() > 0) {
                    int bufferSize = Math.min(2000000, block.getSize() > 100000000 ? 10000000 : (int)block.getSize());
                    stream = new PositionalBufferedStream(new BlockStreamWrapper(mySeekableStream, block), bufferSize);
                    // note we don't have to skip the header here as the block should never start in the header
                    return;
                }
            }

            // If we get here the blocks are exhausted, set reader to null
            if ( stream != null ) {
                stream.close();
                stream = null;
            }
        }

        /**
         * Advance to the next record in the query interval.
         *
         * @throws IOException
         */
        private void readNextRecord() throws IOException {

            if(stream == null) {
                return;  // <= no more features to read
            }

            currentRecord = null;

            while (true) {   // Loop through blocks
                while (!stream.isDone()) {  // Loop through current block
                    Feature f = null;
                    try {
                        f = codec.decode(stream);
                        if (f == null) {
                            continue;   // Skip
                        }
                        if ((chrAlias != null && !f.getChr().equals(chrAlias)) || f.getStart() > end) {
                            if (blockIterator.hasNext()) {
                                advanceBlock();
                                continue;
                            } else {
                                return;    // Done
                            }
                        }
                        if (f.getEnd() < start) {
                            continue;   // Skip
                        }

                        currentRecord = (T) f;     // Success
                        return;

                    } catch (TribbleException e) {
                        e.setSource(path);
                        throw e;
                    } catch (NumberFormatException e) {
                        String error = "Error parsing line: " + stream.getPosition();
                        throw new TribbleException.MalformedFeatureFile(error, path, e);
                    }
                }
                if (blockIterator != null && blockIterator.hasNext()) {
                    advanceBlock();   // Advance to next block
                } else {
                    return;   // No blocks left, we're done.
                }
            }
        }


        public void remove() {
            throw new UnsupportedOperationException("Remove is not supported.");
        }


        public void close() {
            // Note that this depends on BlockStreamWrapper not actually closing the underlying stream
            if ( stream != null ) stream.close();

            if ( ! reuseStreamInQuery() ) {
                // if we are going to reuse the underlying stream we don't close the underlying stream.
                try {
                    mySeekableStream.close();
                } catch (IOException e) {
                    throw new TribbleException("Couldn't close seekable stream", e);
                }
            }
        }

        public Iterator<T> iterator() {
            return this;
        }
    }


    /**
     * Wrapper around a SeekableStream that limits reading to the specified "block" of bytes.  Attempts to
     * read beyond the end of the block should return -1  (EOF).
     */
    static class BlockStreamWrapper extends InputStream {

        SeekableStream seekableStream;
        long maxPosition;

        BlockStreamWrapper(SeekableStream seekableStream, Block block) throws IOException {
            this.seekableStream = seekableStream;
            seekableStream.seek(block.getStartPosition());
            maxPosition = block.getEndPosition();
        }

        @Override
        public int read() throws IOException {
            return (seekableStream.position() > maxPosition) ? -1 : seekableStream.read();
        }

        @Override
        public int read(byte[] bytes, int off, int len) throws IOException {
            // note the careful treatment here to ensure we can continue to
            // read very long > Integer sized blocks
            long maxBytes = maxPosition - seekableStream.position();
            if (maxBytes <= 0) {
                return -1;
            }

            int bytesToRead = (int)Math.min(len, Math.min(maxBytes, Integer.MAX_VALUE));
            return seekableStream.read(bytes, off, bytesToRead);

        }
    }


}
