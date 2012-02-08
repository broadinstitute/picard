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
 * BCL Files are base call and quality score binary files containing a (base,quality) pair for successive clusters.
 * The file is structured as followed:
 *  Bytes 1-4 : unsigned int numClusters
 *  Bytes 5-numClusters + 5 : 1 byte base/quality score
 *
 *  The base/quality scores are organized as follows (with one exception, SEE BELOW):
 *  The right 2 most bits (these are the LEAST significant bits) indicate the base, where
 *  A=00(0x00), C=01(0x01), G=10(0x02), and T=11(0x03)
 *
 *  The remaining bytes compose the quality score which is an unsigned int.
 *
 *  EXCEPTION: If a byte is entirely 0 (e.g. byteRead == 0) then it is a no call, the base
 *  becomes '.' and the Quality becomes 2, the default illumina masking value
 *
 *  (E.g. if we get a value in binary of 10001011 it gets transformed as follows:
 *
 *  Value read: 10001011(0x8B)
 *
 *  Quality     Base
 *
 *  100010      11
 *  00100010    0x03
 *  0x22        T
 *  34          T
 *
 *  So the output base/quality will be a (T/34)
 */
public class BclReader implements Iterator<BclReader.BclValue> {

    /** The byte iteterator created by fileReader that reads successive byte values (after the header) */
    private final BinaryFileIterator<Byte> bbIterator;

    /** The size of the opening header (consisting solely of numClusters*/
    private static final int HEADER_SIZE = 4;

    /** The number of clusters provided in this BCL */
    public final long numClusters;

    /* This strips off the leading bits in order to compare the read byte against the static values below */
    private final static byte BASE_MASK = 0x0003;
    private final static byte A_VAL = 0x00;
    private final static byte C_VAL = 0x01;
    private final static byte G_VAL = 0x02;
    private final static byte T_VAL = 0x03;

    /** The index to the next cluster that will be returned by this reader */
    private long nextCluster;

    public class BclValue {
        public final byte base;
        public final byte quality;

        public BclValue(final byte base, final byte quality) {
            this.base = base;
            this.quality = quality;
        }
    }

    public BclReader(final File file) {
        bbIterator = MMapBackedIteratorFactory.getByteIterator(HEADER_SIZE, file);
        final ByteBuffer headerBuf = bbIterator.getHeaderBytes();
        numClusters = UnsignedTypeUtil.uIntToLong(headerBuf.getInt());
        bbIterator.assertTotalElementsEqual(numClusters);
        nextCluster = 0;
    }

    public boolean hasNext() {
        return nextCluster < numClusters;
    }

    public BclValue next() {
        final byte base;
        final byte quality;

        final byte element = bbIterator.next();

        if(element == 0) { //NO CALL, don't confuse with an A call
            base = '.';
            quality = 2;
        } else {
            switch(element & BASE_MASK) {
                case A_VAL:
                    base = 'A';
                    break;

                case C_VAL:
                    base = 'C';
                    break;

                case G_VAL:
                    base = 'G';
                    break;

                case T_VAL:
                    base = 'T';
                    break;

                default:
                    throw new PicardException("Impossible case! BCL Base value neither A, C, G, nor T! Value(" + (element & BASE_MASK) + ") + in file(" + bbIterator.getFile().getAbsolutePath() + ")");
            }

            quality = (byte)(UnsignedTypeUtil.uByteToInt(element) >>> 2);
            if(quality == 0 || quality == 1) {
                throw new PicardException("If base is NOT a NO CALL then it should have a quality of 2 or greater!  Quality Found(" + quality + ")  Cluster(" + nextCluster + ")");
            }
        }

        ++nextCluster;
        return new BclValue(base, quality);
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }
}

