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
package net.sf.samtools;

import java.nio.ByteBuffer;

/**
 * Converter between disk and in-memory (object, not String) CIGAR representation.
 */
class BinaryCigarCodec {
    private static final BinaryCigarCodec singleton = new BinaryCigarCodec();

    /**
     * It is not necssary to get the singleton but it is preferrable to use the same one
     * over and over vs. creating a new object for each BAMRecord.  This class has no state
     * so this is thread-safe.
     */
    static BinaryCigarCodec getSingleton() {
        return singleton;
    }

    /**
     * Convert CIGAR from object representation to disk representation.
     * @return Array of unsigned ints, one for each element of CIGAR.
     */
    int[] encode(final Cigar cigar) {
        if (cigar.numCigarElements() == 0) {
            return new int[0];
        }

        // Binary rep can be no longer than 1/2 of text rep
        // Although this is documented as uint, I think lengths will never get that long,
        // and it's a pain in Java.
        final int[] binaryCigar = new int[cigar.numCigarElements()];
        int binaryCigarLength = 0;
        for (int i = 0; i < cigar.numCigarElements(); ++i) {
            final CigarElement cigarElement = cigar.getCigarElement(i);
            final int op = CigarOperator.enumToBinary(cigarElement.getOperator());
            binaryCigar[binaryCigarLength++] = cigarElement.getLength() << 4 | op;
        }
        return binaryCigar;
    }

    /**
     * Convert CIGAR from disk representation to object.
     * @param binaryCigar ByteArray that is assumed to have byte order set appropriately for extracting ints.
     */
    Cigar decode(final ByteBuffer binaryCigar) {
        final Cigar ret = new Cigar();
        while (binaryCigar.hasRemaining()) {
            final int cigarette = binaryCigar.getInt();
            ret.add(binaryCigarToCigarElement(cigarette));
        }
        return ret;
    }

    /**
     * Convert CIGAR from disk representation to object.
     * @param binaryCigar Array of unsigned ints, one for each CIGAR element.
     */
    Cigar decode(final int[] binaryCigar) {
        final Cigar ret = new Cigar();
        for (final int cigarette : binaryCigar) {
            ret.add(binaryCigarToCigarElement(cigarette));
        }
        return ret;
    }

    /**
     * @param cigarette CIGAR element (operator + length) encoded as an unsigned int.
     * @return Object representation of the CIGAR element.
     */
    private static CigarElement binaryCigarToCigarElement(final int cigarette) {
        final int binaryOp = cigarette & 0xf;
        final int length = cigarette >> 4;
        return new CigarElement(length, CigarOperator.binaryToEnum(binaryOp));
    }
}
