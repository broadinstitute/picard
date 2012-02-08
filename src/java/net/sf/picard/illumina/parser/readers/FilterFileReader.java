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
package net.sf.picard.illumina.parser.readers;

import net.sf.picard.PicardException;
import net.sf.picard.util.UnsignedTypeUtil;

import java.io.File;
import java.nio.ByteBuffer;
import java.util.Iterator;

/**
 * Illumina uses an algorithm described in "Theory of RTA" that determines whether or not a cluster passes filter("PF") or not.
 * These values are written as sequential bytes to Filter Files.  The structure of a filter file is as follows:
 * Bytes 0-3  : 0
 * Bytes 4-7  : unsigned int version
 * Bytes 8-11 : unsigned int numClusters
 */
public class FilterFileReader implements Iterator<Boolean> {
    /** Number of bytes in the files header that will be skipped by the iterator*/
    private static final int HEADER_SIZE  = 12;

    /** Expected Version */
    public final int EXPECTED_VERSION = 3;

    /** Iterator over each cluster in the FilterFile */
    private final BinaryFileIterator<Byte> bbIterator;

    /** Version number found in the FilterFile, this should equal 3 */
    public final int version;

    /** The number of cluster's pf values stored in this file */
    public final long numClusters;

    /** Byte representing a cluster failing filter(not a PF read), we test this exactly at
     * the moment but technically the standard  may be to check only lowest significant bit */
    private final static byte FailedFilter = 0x00;

    /** Byte representing a cluster passing filter(a PF read), we test this exactly at
     * the moment but technically the standard  may be to check only lowest significant bit */
    private final static byte PassedFilter = 0x01;

    /** The index of the current cluster within the file*/
    private int currentCluster;

    public FilterFileReader(final File file) {
        bbIterator = MMapBackedIteratorFactory.getByteIterator(HEADER_SIZE, file);
        final ByteBuffer headerBuf = bbIterator.getHeaderBytes();

        for(int i = 0; i < 4; i++) {
            final byte b = headerBuf.get();
            if(b != 0) {
                throw new PicardException("The first four bytes of a Filter File should be 0 but byte " + i + " was " + b + " in file " + file.getAbsolutePath());
            }
        }

        version = headerBuf.getInt();
        if(version != EXPECTED_VERSION) {
            throw new PicardException("Expected version is " + EXPECTED_VERSION + " but version found was "  + version + " in file " + file.getAbsolutePath());
        }

        numClusters = UnsignedTypeUtil.uIntToLong(headerBuf.getInt());
        bbIterator.assertTotalElementsEqual(numClusters);

        currentCluster = 0;
    }

    public boolean hasNext() {
        return currentCluster < numClusters;
    }

    public Boolean next() {
        final byte value = bbIterator.next();
        currentCluster += 1;
        if(value == PassedFilter) {
            return true;
        } else if(value == FailedFilter) {
            return false;
        } else {
            String hexVal = Integer.toHexString(value);
            hexVal = (hexVal.length() < 2 ? "0x0" : "0x") + hexVal;
            throw new PicardException("Didn't recognized PF Byte (" + hexVal + ")" + " for element (" + currentCluster + ") in file(" + bbIterator.getFile().getAbsolutePath() + ")");
        }
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}
