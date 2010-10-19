/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.FileAppendStreamLRUCache;
import net.sf.samtools.util.CloserUtil;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Holds info about a mate pair for use when processing a coordinate sorted file.  When one read of a pair is encountered,
 * the caller should add a record to this map.  When the other read of a pair is encountered, the record should be removed.
 * This class assumes that reads will be processed in order of reference sequence index.  When the map is queried for
 * a record for a given reference sequence index, all the records for that sequence are loaded from temp file into RAM, so there
 * must be sufficient RAM to hold all the records for one reference sequence.  If the records are not processed in
 * reference sequence order, loading and unloading of records will cause performance to be terrible.
 * @param <KEY> KEY + reference sequence index are used to identify the record being stored or retrieved.
 * @param <REC> The type of record being retrieved.
 */
public class CoordinateSortedPairInfoMap<KEY, REC> {
    /**
     * directory where files will go
     */
    private final File workDir = IoUtil.createTempDir("CSPI.", null);
    private int sequenceIndexOfMapInRam = -1;
    private Map<KEY, REC> mapInRam = null;
    private final FileAppendStreamLRUCache outputStreams;
    private final List<Integer> sizeOfMapOnDisk = new ArrayList<Integer>();
    private final Codec<KEY, REC> elementCodec;
    private final List<File> mapsOnDisk = new ArrayList<File>();

    CoordinateSortedPairInfoMap(int maxOpenFiles, Codec<KEY, REC> elementCodec) {
        this.elementCodec = elementCodec;
        workDir.deleteOnExit();
        outputStreams = new FileAppendStreamLRUCache(maxOpenFiles);
    }

    /**
     *
     * @param sequenceIndex
     * @param key
     * @return The record corresponding to the given sequenceIndex and key, or null if it is not present.
     */
    public REC remove(final int sequenceIndex, final String key) {
        ensureSequenceLoaded(sequenceIndex);
        return mapInRam.remove(key);
    }

    private void ensureSequenceLoaded(final int sequenceIndex) {
        try {
            if (sequenceIndexOfMapInRam == sequenceIndex) {
                return;
            }

            // Spill map in RAM to disk
            if (mapInRam != null) {
                final OutputStream os = getOutputStreamForSequence(sequenceIndexOfMapInRam);
                elementCodec.setOutputStream(os);
                for (final Map.Entry<KEY, REC> entry : mapInRam.entrySet()) {
                    elementCodec.encode(entry.getKey(), entry.getValue());
                }
                sizeOfMapOnDisk.set(sequenceIndexOfMapInRam, mapInRam.size());
                mapsOnDisk.set(sequenceIndexOfMapInRam, null);
                mapInRam.clear();
            } else {
                mapInRam = new HashMap<KEY, REC>();
            }

            sequenceIndexOfMapInRam = sequenceIndex;

            // Load map from disk if it existed
            if (mapsOnDisk.size() > sequenceIndex) {
                File mapOnDisk = mapsOnDisk.get(sequenceIndex);
                if (outputStreams.containsKey(mapOnDisk)) {
                    outputStreams.remove(mapOnDisk).close();
                }
                final int numRecords = sizeOfMapOnDisk.get(sequenceIndex);
                sizeOfMapOnDisk.set(sequenceIndex, 0);
                final File file = makeFileForSequence(sequenceIndex);
                if (file.exists()) {
                    FileInputStream is = null;
                    try {
                        is = new FileInputStream(file);
                        elementCodec.setInputStream(is);
                        for (int i = 0; i < numRecords; ++i) {
                            final KeyAndRecord<KEY, REC> keyAndRecord = elementCodec.decode();
                            if (mapInRam.containsKey(keyAndRecord.getKey()))
                                throw new PicardException("Value was put into PairInfoMap more than once.  " +
                                        sequenceIndex + ": " + keyAndRecord.getKey());
                            mapInRam.put(keyAndRecord.getKey(), keyAndRecord.getRecord());
                        }
                    } finally {
                        CloserUtil.close(is);
                    }
                    net.sf.samtools.util.IOUtil.deleteFiles(file);
                }
            }
        } catch (IOException e) {
            throw new PicardException("Error loading new map from disk.", e);
        }
    }

    /**
     * Store the record with the given sequence index and key.  It is assumed that value did not previously exist
     * in the map, and an exception is thrown (possibly at a later time) if that is not the case.
     * @param sequenceIndex
     * @param key
     * @param record
     */
    public void put(final int sequenceIndex, final KEY key, final REC record) {
        if (sequenceIndex == sequenceIndexOfMapInRam) {
            // Store in RAM map
            if (mapInRam.containsKey(key))
                throw new IllegalArgumentException("Putting value into PairInfoMap that already existed. " +
                        sequenceIndex + ": " + key);
            mapInRam.put(key, record);
        } else {
            // Append to file
            final OutputStream os = getOutputStreamForSequence(sequenceIndex);
            elementCodec.setOutputStream(os);
            elementCodec.encode(key, record);
            sizeOfMapOnDisk.set(sequenceIndex, sizeOfMapOnDisk.get(sequenceIndex) + 1);
        }
    }

    private File makeFileForSequence(final int index) {
        final File file = new File(workDir, index + ".tmp");
        file.deleteOnExit();
        return file;
    }

    private OutputStream getOutputStreamForSequence(final int mateSequenceIndex) {
        while (mateSequenceIndex >= mapsOnDisk.size()) {
            mapsOnDisk.add(makeFileForSequence(mapsOnDisk.size()));
            sizeOfMapOnDisk.add(0);
        }
        return outputStreams.get(mapsOnDisk.get(mateSequenceIndex));
    }

    public int size() {
        int total = sizeInRam();
        for (final Integer mapSize : sizeOfMapOnDisk) {
            if (mapSize != null) {
                total += mapSize;
            }
        }
        return total;
    }

    /**
     * @return number of elements stored in RAM.  Always <= size()
     */
    public int sizeInRam() {
        return mapInRam != null? mapInRam.size(): 0;
    }

    private static class MapEntry {
        private String key;
        private ReadEnds readEnds;

        public String getKey() {
            return key;
        }

        public void setKey(final String key) {
            this.key = key;
        }

        public ReadEnds getReadEnds() {
            return readEnds;
        }

        public void setReadEnds(final ReadEnds readEnds) {
            this.readEnds = readEnds;
        }
    }

    public interface KeyAndRecord<KEY, REC> {
        KEY getKey();
        REC getRecord();
        
    }
    /**
     * Client must implement this class, which defines the way in which records are written to and
     * read from file.
     */
    public interface Codec<KEY, REC> {
        /**
         * Where to write encoded output
         * @param os
         */
        void setOutputStream(OutputStream os);

        /**
         * Where to read encoded input from
         * @param is
         */
        void setInputStream(InputStream is);

        /**
         * Write object to output stream.  If the key is part of the record, then there is no need to write
         * it separately.
         */
        void encode(KEY key, REC record);

        /**
         * Read the next key and record from the input stream and convert into a java object.
         * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
         * a record.
         */
        KeyAndRecord<KEY, REC> decode();

    }
}
