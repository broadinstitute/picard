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
    private final CoordinateSortedPairInfoMap<String, ReadEnds> pairInfoMap;
    DiskReadEndsMap(int maxOpenFiles) {
        pairInfoMap = new CoordinateSortedPairInfoMap<String, ReadEnds>(maxOpenFiles, new Codec());
    }

    public ReadEnds remove(int mateSequenceIndex, String key) {
        return pairInfoMap.remove(mateSequenceIndex, key);
    }

    public void put(int mateSequenceIndex, String key, ReadEnds readEnds) {
        pairInfoMap.put(mateSequenceIndex, key, readEnds);
    }

    public int size() {
        return pairInfoMap.size();
    }

    public int sizeInRam() {
        return pairInfoMap.sizeInRam();
    }

    private static class KeyAndRecord implements CoordinateSortedPairInfoMap.KeyAndRecord<String, ReadEnds> {
        private final String key;
        private final ReadEnds record;

        private KeyAndRecord(String key, ReadEnds record) {
            this.key = key;
            this.record = record;
        }

        public String getKey() {
            return key;
        }

        public ReadEnds getRecord() {
            return record;
        }
    }

    private static class Codec implements CoordinateSortedPairInfoMap.Codec<String, ReadEnds> {
        private final ReadEndsCodec readEndsCodec = new ReadEndsCodec();

        public void setInputStream(final InputStream is) {
            readEndsCodec.setInputStream(is);
        }

        public void setOutputStream(final OutputStream os) {
            readEndsCodec.setOutputStream(os);
        }

        public CoordinateSortedPairInfoMap.KeyAndRecord<String, ReadEnds> decode() {
            try {
                final String key = readEndsCodec.getInputStream().readUTF();
                final ReadEnds record = readEndsCodec.decode();
                return new KeyAndRecord(key, record);
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
