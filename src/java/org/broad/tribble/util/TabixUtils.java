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
package org.broad.tribble.util;


import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.BlockCompressedInputStream;
import org.broad.tribble.TribbleException;
import org.broad.tribble.readers.TabixReader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * classes that have anything to do with tabix
 */
public class TabixUtils {

    public static final String STANDARD_INDEX_EXTENSION = ".tbi";

    public static class TPair64 implements Comparable<TPair64> {
        public long u, v;

        public TPair64(final long _u, final long _v) {
            u = _u;
            v = _v;
        }

        public TPair64(final TPair64 p) {
            u = p.u;
            v = p.v;
        }

        public int compareTo(final TPair64 p) {
            return u == p.u ? 0 : ((u < p.u) ^ (u < 0) ^ (p.u < 0)) ? -1 : 1; // unsigned 64-bit comparison
        }
    }

    public static class TIndex {
        public HashMap<Integer, TPair64[]> b; // binning index
        public long[] l; // linear index
    }


    public static class TIntv {
        public int tid, beg, end;
    }


    public static boolean less64(final long u, final long v) { // unsigned 64-bit comparison
        return (u < v) ^ (u < 0) ^ (v < 0);
    }

    /**
     * Generates the SAMSequenceDictionary from the given tabix index file
     *
     * @param tabixIndex the tabix index file
     * @return non-null sequence dictionary
     */
    public static SAMSequenceDictionary getSequenceDictionary(final File tabixIndex) {
        if (tabixIndex == null) throw new IllegalArgumentException();

        try {
            final BlockCompressedInputStream is = new BlockCompressedInputStream(tabixIndex);

            // read preliminary bytes
            byte[] buf = new byte[32];
            is.read(buf, 0, 32);

            // read sequence dictionary
            int i, j, len = TabixReader.readInt(is);
            buf = new byte[len];
            is.read(buf);

            final List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>();
            for (i = j = 0; i < buf.length; ++i) {
                if (buf[i] == 0) {
                    byte[] b = new byte[i - j];
                    System.arraycopy(buf, j, b, 0, b.length);
                    sequences.add(new SAMSequenceRecord(new String(b)));
                    j = i + 1;
                }
            }
            is.close();

            return new SAMSequenceDictionary(sequences);
        } catch (Exception e) {
            throw new TribbleException("Unable to read tabix index: " + e.getMessage());
        }
    }
}
