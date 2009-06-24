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
package net.sf.picard.sam;

import net.sf.picard.io.IoUtil;
import net.sf.picard.PicardException;
import net.sf.samtools.util.CloserUtil;

import java.io.*;
import java.util.Map;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Disk-based implementation of ReadEndsMap.  A subdirectory of the system tmpdir is created to store
 * files, one for each reference sequence.  The reference sequence that is currently being queried (i.e. the
 * sequence for which remove() has been most recently called) is stored in RAM.  ReadEnds for all other sequences
 * are stored on disk.
 *
 * When put() is called for a sequence that is the current one in RAM, the ReadEnds object is merely put into the
 * in-memory map.  If put() is called for a sequence ID that is not the current RAM one, the ReadEnds object is
 * appended to the file for that sequence, creating the file if necessary.
 *
 * When remove() is called for a sequence that is the current one in RAM, remove() is called on the in-memory map.
 * If remove() is called for a sequence other than the current RAM sequence, then the current RAM sequence is written
 * to disk, the new sequence is read from disk into RAM map, and the file for the new sequence is deleted.
 *
 * If things work properly, and reads are processed in genomic order, records will be written for mates that are in
 * a later sequence.  When the mate is reached in the input SAM file, the file that was written will be deleted.
 * This should result in all temporary files being deleted by the time all the reads are processed.  The temp
 * directory is marked to be deleted on exit so everything should get cleaned up.
 * 
 * @author alecw@broadinstitute.org
 */
class DiskReadEndsMap implements ReadEndsMap {

    /**
     * directory where files will go
     */
    private final File workDir = IoUtil.createTempDir("DREM.", null);
    private int sequenceIndexOfMapInRam = -1;
    private Map<String, ReadEnds> mapInRam = null;

    /**
     * When a ReadEnds is put for a later sequence, an element is created in this list if not already
     * present, and the ReadEnds is appended to the file.  When remove() is called for a particular sequence,
     * the file is closed, read into RAM, deleted, and the element in this list corresponding to the sequence is
     * set to null.
     */
    private final List<FileOutputStream> mapsOnDisk = new ArrayList<FileOutputStream>();

    /**
     * As ReadEnds objects are written to disk, the number of objects in each file on disk is tracked, so
     * that size() method works properly.
     */
    private final List<Integer> sizeOfMapOnDisk = new ArrayList<Integer>();
    final private ElementCodec elementCodec = new ElementCodec();

    DiskReadEndsMap() {
        workDir.deleteOnExit();
    }

    public ReadEnds remove(final int mateSequenceIndex, final String key) {
        ensureSequenceLoaded(mateSequenceIndex);
        return mapInRam.remove(key);
    }

    private void ensureSequenceLoaded(final int mateSequenceIndex) {
        try {
            if (sequenceIndexOfMapInRam == mateSequenceIndex) {
                return;
            }

            // Spill map in RAM to disk
            if (mapInRam != null) {
                if (mapsOnDisk.size() > sequenceIndexOfMapInRam && mapsOnDisk.get(sequenceIndexOfMapInRam) != null) {
                    // Should always be null for the sequence that is currently in RAM.
                    throw new IllegalStateException("Internal error");
                }
                final OutputStream os = getOutputStreamForSequence(sequenceIndexOfMapInRam);
                elementCodec.setOutputStream(os);
                for (final Map.Entry<String, ReadEnds> entry : mapInRam.entrySet()) {
                    elementCodec.encode(entry.getKey(), entry.getValue());
                }
                sizeOfMapOnDisk.set(sequenceIndexOfMapInRam, mapInRam.size());
                mapsOnDisk.set(sequenceIndexOfMapInRam, null);
                mapInRam.clear();
            } else {
                mapInRam = new HashMap<String, ReadEnds>();
            }

            sequenceIndexOfMapInRam = mateSequenceIndex;

            // Load map from disk if it existed
            if (mapsOnDisk.size() > mateSequenceIndex && mapsOnDisk.get(mateSequenceIndex) != null) {
                mapsOnDisk.get(mateSequenceIndex).close();
                mapsOnDisk.set(mateSequenceIndex, null);
                final int numRecords = sizeOfMapOnDisk.get(mateSequenceIndex);
                sizeOfMapOnDisk.set(mateSequenceIndex, 0);
                final File file = makeFileForSequence(mateSequenceIndex);
                FileInputStream is = null;
                try {
                    is = new FileInputStream(file);
                    elementCodec.setInputStream(is);
                    for (int i = 0; i < numRecords; ++i) {
                        final MapEntry entry = elementCodec.decode();
                        mapInRam.put(entry.key, entry.readEnds);
                    }
                } finally {
                    CloserUtil.close(is);
                }
                net.sf.samtools.util.IOUtil.deleteFiles(file);
            }
        } catch (IOException e) {
            throw new PicardException("Error loading new map from disk.", e);
        }
    }

    public void put(final int mateSequenceIndex, final String key, final ReadEnds readEnds) {
        if (mateSequenceIndex == sequenceIndexOfMapInRam) {
            // Store in RAM map
            mapInRam.put(key, readEnds);
        } else {
            // Append to file
            final OutputStream os = getOutputStreamForSequence(mateSequenceIndex);
            elementCodec.setOutputStream(os);
            elementCodec.encode(key, readEnds);
            sizeOfMapOnDisk.set(mateSequenceIndex, sizeOfMapOnDisk.get(mateSequenceIndex) + 1);
        }
    }

    private File makeFileForSequence(final int index) {
        final File file = new File(workDir, index + ".read_ends");
        file.deleteOnExit();
        return file;
    }

    private OutputStream getOutputStreamForSequence(final int mateSequenceIndex) {
        try {
            while (mateSequenceIndex >= mapsOnDisk.size()) {
                mapsOnDisk.add(null);
                sizeOfMapOnDisk.add(0);
            }
            FileOutputStream ret = mapsOnDisk.get(mateSequenceIndex);
            if (ret == null) {
                // Create file if it isn't already open.  File should not exist on disk,
                // but if it does due to a bug it is truncated, which is the right thing.
                ret = new FileOutputStream(makeFileForSequence(mateSequenceIndex));
                mapsOnDisk.set(mateSequenceIndex, ret);
            }
            return ret;
        } catch (FileNotFoundException e) {
            throw new PicardException("Error creating temporary ReadEnds file", e);
        }
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

    private static class ElementCodec {
        private final ReadEndsCodec readEndsCodec = new ReadEndsCodec();

        public void setInputStream(final InputStream is) {
            readEndsCodec.setInputStream(is);
        }

        public void setOutputStream(final OutputStream os) {
            readEndsCodec.setOutputStream(os);
        }

        public MapEntry decode() {
            try {
                final MapEntry ret = new MapEntry();
                ret.setKey(readEndsCodec.getInputStream().readUTF());
                ret.setReadEnds(readEndsCodec.decode());
                return ret;
            } catch (IOException e) {
                throw new PicardException("Error loading ReadEndsMap from disk", e);
            }
        }

        public void encode(final String key, final ReadEnds readEnds) {
            try {
                readEndsCodec.getOutputStream().writeUTF(key);
                readEndsCodec.encode(readEnds);
            } catch (IOException e) {
                throw new PicardException("Error spilling ReadEndsMap to disk.", e);
            }
        }
    }
    
}
