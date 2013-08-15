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
package org.broad.tribble.index;

import org.broad.tribble.*;
import org.broad.tribble.index.interval.IntervalIndexCreator;
import org.broad.tribble.index.interval.IntervalTreeIndex;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.index.linear.LinearIndexCreator;
import org.broad.tribble.readers.LocationAware;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.util.ParsingUtils;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.zip.GZIPInputStream;

/**
 * Factory class for creating indexes.  It is the responsibility of this class to determine and create the
 * correct index type from the input file or stream
 */
public class IndexFactory {
    /** We can optimize index-file-creation for different factors. As of this writing, those are index-file size or seeking time. */
    public enum IndexBalanceApproach {
        FOR_SIZE,
        FOR_SEEK_TIME
    }

    /**
     * an enum that contains all of the information about the index types, and how to create them
     */
    public enum IndexType {
        LINEAR(1, LinearIndexCreator.class, LinearIndex.class, LinearIndexCreator.DEFAULT_BIN_WIDTH),
        INTERVAL_TREE(2, IntervalIndexCreator.class, IntervalTreeIndex.class, IntervalIndexCreator.DEFAULT_FEATURE_COUNT),
        TABIX(3, null, null, -1);

        private final int indexValue;
        private final Class<IndexCreator> indexCreatorClass;
        private final int defaultBinSize;
        private final Class<Index> indexType;

        public int getDefaultBinSize() {
            return defaultBinSize;
        }

        public IndexCreator getIndexCreator() {
            try {
                return indexCreatorClass.newInstance();
            } catch ( InstantiationException e ) {
                throw new TribbleException("Couldn't make index creator in " + this, e);
            } catch ( IllegalAccessException e ) {
                throw new TribbleException("Couldn't make index creator in " + this, e);
            }
        }

        public boolean canCreate() {
            return indexCreatorClass != null;
        }

        IndexType(final int headerValue, final Class creator, final Class indexClass, final int defaultBinSize) {
            indexValue = headerValue;
            indexCreatorClass = creator;
            indexType = indexClass;
            this.defaultBinSize = defaultBinSize;
        }

        public int getHeaderValue() {
            return indexValue;
        }

        public Class getIndexType() {
            return indexType;
        }

        /**
         *
         * @param headerValue
         * @return The {@code IndexType} based on the {@code headerValue}
         * @throws TribbleException.UnableToCreateCorrectIndexType
         */
        public static IndexType getIndexType(final int headerValue) {
            for (final IndexType type : IndexType.values())
                if (type.indexValue == headerValue) return type;
            throw new TribbleException.UnableToCreateCorrectIndexType("Unknown index type value" + headerValue);
        }
    }


    /**
     * Load in index from the specified file.   The type of index (LinearIndex or IntervalTreeIndex) is determined
     * at run time by reading the type flag in the file.
     *
     * @param indexFile from which to load the index
     * @return the parsed index object
     */
    public static Index loadIndex(final String indexFile) {
        Index idx = null;
        InputStream is = null;
        LittleEndianInputStream dis = null;
        try {
            if (indexFile.endsWith(".gz")) {
                is = new BufferedInputStream(new GZIPInputStream(ParsingUtils.openInputStream(indexFile)), 512000);
            } else {
                is = new BufferedInputStream(ParsingUtils.openInputStream(indexFile), 512000);
            }

            dis = new LittleEndianInputStream(is);

            // Read the type and version,  then create the appropriate type
            final int magicNumber = dis.readInt();
            final int type = dis.readInt();
            final Class indexClass = IndexType.getIndexType(type).getIndexType();

            idx = (Index) indexClass.newInstance();
            idx.read(dis);
        } catch (IOException ex) {
            throw new TribbleException.UnableToReadIndexFile("Unable to read index file", indexFile, ex);
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        } finally {
            try {
                if (is != null) is.close();
                if (dis != null) dis.close();
                //log.info(String.format("Closed %s and %s", is, dis));
            } catch (IOException e) {
                //log.error("Error closing indexFile: " + indexFile, e);
            }
        }
        return idx;
    }


    /**
     * a helper method for creating a linear binned index with default bin size
     *
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @return a index
     */
    public static Index createLinearIndex(final File inputFile, final FeatureCodec codec) {
        return createIndex(inputFile, codec, IndexType.LINEAR);
    }

    /**
     * a helper method for creating a linear binned index
     *
     * @see #createIndex(File, FeatureCodec, IndexType, int)
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @param binSize   the bin size
     * @return a index
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createLinearIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec, final int binSize) {
        return createIndex(inputFile, codec, IndexType.LINEAR, binSize);
    }

    /**
     * create an interval-tree index with the default features per bin count
     *
     * @param inputFile the file containing the features
     * @param codec to decode the features
     * @return
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createIntervalIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec) {
        return createIndex(inputFile, codec, IndexType.INTERVAL_TREE);
    }


    /**
     * a helper method for creating an interval-tree index
     *
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @param binSize   the bin size
     * @return a index
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createIntervalIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec, final int binSize) {
        return createIndex(inputFile, codec, IndexType.INTERVAL_TREE, binSize);
    }

    /**
     * Create a dynamic index with the default balancing approach
     *
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @return a index
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createDynamicIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec) {
        return createDynamicIndex(inputFile, codec, IndexBalanceApproach.FOR_SEEK_TIME);
    }

    /**
     * Create a index of the specified type with default binning parameters
     *
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @param type      the type of index to create
     * @return a index
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec, final IndexType type) {
        return createIndex(inputFile, codec, type, type.getDefaultBinSize());
    }

    /**
     * Create a index of the specified type
     *
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @param type      the type of index to create
     * @param binSize   the bin size
     * @return a index
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec, final IndexType type, final int binSize) {
        if ( ! type.canCreate() )
            throw new TribbleException("Tribble can only read, not create indices of type " + type.name());

        final IndexCreator idx = type.getIndexCreator();
        idx.initialize(inputFile, binSize);
        return createIndex(inputFile, new FeatureIterator<FEATURE_TYPE, SOURCE_TYPE>(inputFile, codec), idx);
    }

    /**
     * Write the index to a file; little endian.
     * @param idx
     * @param idxFile
     * @throws IOException
     */
    public static void writeIndex(final Index idx, final File idxFile) throws IOException {
        LittleEndianOutputStream stream = null;
        try {
            stream = new LittleEndianOutputStream(new BufferedOutputStream(new FileOutputStream(idxFile)));
            idx.write(stream);
        }
        finally {
            if(stream != null) {
                stream.close();
            }
        }
    }

    /**
     * create a dynamic index, given an input file, codec, and balance approach
     *
     * @param inputFile the input file to load features from
     * @param codec     the codec to use for decoding records
     * @param iba       the index balancing approach
     * @return a index
     */
    public static <FEATURE_TYPE extends Feature, SOURCE_TYPE> Index createDynamicIndex(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE_TYPE> codec, final IndexBalanceApproach iba) {
        // get a list of index creators
        final DynamicIndexCreator indexCreator = new DynamicIndexCreator(iba);
        indexCreator.initialize(inputFile, indexCreator.defaultBinSize());
        return createIndex(inputFile, new FeatureIterator<FEATURE_TYPE, SOURCE_TYPE>(inputFile, codec), indexCreator);
    }


    private static Index createIndex(final File inputFile, final FeatureIterator iterator, final IndexCreator creator) {
        Feature lastFeature = null;
        Feature currentFeature;
        final Map<String, Feature> visitedChromos = new HashMap<String, Feature>(40);
        while (iterator.hasNext()) {
            final long position = iterator.getPosition();
            currentFeature = iterator.next();

            checkSorted(inputFile, lastFeature, currentFeature);
            //should only visit chromosomes once
            final String curChr = currentFeature.getChr();
            final String lastChr = lastFeature != null ? lastFeature.getChr() : null;
            if(!curChr.equals(lastChr)){
                if(visitedChromos.containsKey(curChr)){
                    String msg = "Input file must have contiguous chromosomes.";
                    msg += " Saw feature " + featToString(visitedChromos.get(curChr));
                    msg += " followed later by " + featToString(lastFeature);
                    msg += " and then " + featToString(currentFeature);
                    throw new TribbleException.MalformedFeatureFile(msg, inputFile.getAbsolutePath());
                }else{
                    visitedChromos.put(curChr, currentFeature);
                }
            }

            creator.addFeature(currentFeature, position);

            lastFeature = currentFeature;
        }

        iterator.close();
        return creator.finalizeIndex(iterator.getPosition());
    }

    private static String featToString(final Feature feature){
        return feature.getChr() + ":" + feature.getStart() + "-" + feature.getEnd();
    }

    private static void checkSorted(final File inputFile, final Feature lastFeature, final Feature currentFeature){
        // if the last currentFeature is after the current currentFeature, exception out
        if (lastFeature != null && currentFeature.getStart() < lastFeature.getStart() && lastFeature.getChr().equals(currentFeature.getChr()))
            throw new TribbleException.MalformedFeatureFile("Input file is not sorted by start position. \n" +
                    "We saw a record with a start of " + currentFeature.getChr() + ":" + currentFeature.getStart() +
                    " after a record with a start of " + lastFeature.getChr() + ":" + lastFeature.getStart(), inputFile.getAbsolutePath());
    }


    /**
     * Iterator for reading features from a file, given a {@code FeatureCodec}.
     */
    static class FeatureIterator<FEATURE_TYPE extends Feature, SOURCE> implements CloseableTribbleIterator<Feature> {
        // the stream we use to get features
        private SOURCE source;
        // the next feature
        private Feature nextFeature;
        // our codec
        private final FeatureCodec<FEATURE_TYPE, SOURCE> codec;
        private final File inputFile;

        // we also need cache our position
        private long cachedPosition;

        /**
         *
         * @param inputFile The file from which to read. Stream for reading is opened on construction.
         * @param codec
         */
        public FeatureIterator(final File inputFile, final FeatureCodec<FEATURE_TYPE, SOURCE> codec) {
            this.codec = codec;
            this.inputFile = inputFile;
            final FeatureCodecHeader header = readHeader();
            source = (SOURCE) codec.makeIndexableSourceFromStream(initStream(inputFile, header.getHeaderEnd()));
            readNextFeature();
        }

        /**
         * Some codecs,  e.g. VCF files,  need the header to decode features.  This is a rather poor design,
         * the internal header is set as a side-affect of reading it, but we have to live with it for now.
         */
        private FeatureCodecHeader readHeader() {
            try {
                final SOURCE source = this.codec.makeSourceFromStream(initStream(inputFile, 0));
                final FeatureCodecHeader header = this.codec.readHeader(source);
                codec.close(source);
                return header;
            } catch (IOException e) {
                throw new TribbleException.InvalidHeader("Error reading header " + e.getMessage());
            }
        }

        private PositionalBufferedStream initStream(final File inputFile, final long skip) {
            try {
                final FileInputStream is = new FileInputStream(inputFile);
                final PositionalBufferedStream pbs = new PositionalBufferedStream(is);
                if ( skip > 0 ) pbs.skip(skip);
                return pbs;
            } catch (FileNotFoundException e) {
                throw new TribbleException.FeatureFileDoesntExist("Unable to open the input file, most likely the file doesn't exist.", inputFile.getAbsolutePath());
            } catch (IOException e) {
                throw new TribbleException.MalformedFeatureFile("Error initializing stream", inputFile.getAbsolutePath(), e);
            }
        }

        public boolean hasNext() {
            return nextFeature != null;
        }

        public Feature next() {
            final Feature ret = nextFeature;
            readNextFeature();
            return ret;
        }

        /**
         * @throws UnsupportedOperationException
         */
        public void remove() {
            throw new UnsupportedOperationException("We cannot remove");
        }


        /**
         * @return the file position from the underlying reader
         */
        public long getPosition() {
            return (hasNext()) ? cachedPosition : ((LocationAware) source).getPosition();
        }

        @Override
        public Iterator<Feature> iterator() {
            return this;
        }

        @Override
        public void close() {
            codec.close(source);
        }

        /**
         * Read the next feature from the stream
         * @throws TribbleException.MalformedFeatureFile
         */
        private void readNextFeature() {
            cachedPosition = ((LocationAware) source).getPosition();
            try {
                nextFeature = null;
                while (nextFeature == null && !codec.isDone(source)) {
                    nextFeature = codec.decodeLoc(source);
                }
            } catch (IOException e) {
                throw new TribbleException.MalformedFeatureFile("Unable to read a line from the file", inputFile.getAbsolutePath(), e);
            }
        }
    }
}
