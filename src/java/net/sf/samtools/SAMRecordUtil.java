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

import net.sf.samtools.util.StringUtil;

/**
 * @author alecw@broadinstitute.org
 */
public class SAMRecordUtil {
    /** Byte typed variables for all normal bases. */
    private static final byte a='a', c='c', g='g', t='t', A='A', C='C', G='G', T='T';

    /** Returns the complement of a single byte. */
    public static byte complement(final byte b) {
        switch (b) {
            case a: return t;
            case c: return g;
            case g: return c;
            case t: return a;
            case A: return T;
            case C: return G;
            case G: return C;
            case T: return A;
            default: return b;
        }
    }

    /** Reverses and complements the bases in place. */
    public static void reverseComplement(final byte[] bases) {
        final int lastIndex = bases.length - 1;

        int i, j;
        for (i=0, j=lastIndex; i<j; ++i, --j) {
            final byte tmp = complement(bases[i]);
            bases[i] = complement(bases[j]);
            bases[j] = tmp;
        }
        if (bases.length % 2 == 1) {
            bases[i] = complement(bases[i]);
        }
    }

    /**
     * Reverse-complement all known sequence and base quality attributes of the SAMRecord.
     */
    public static void reverseComplement(final SAMRecord rec) {
        final byte[] readBases = rec.getReadBases();
        reverseComplement(readBases);
        rec.setReadBases(readBases);
        final byte qualities[] = rec.getBaseQualities();
        reverseArray(qualities);
        rec.setBaseQualities(qualities);
        final byte[] sqTagValue = (byte[])rec.getAttribute(SAMTagUtil.getSingleton().SQ);
        if (sqTagValue != null) {
            SQTagUtil.reverseComplementSqArray(sqTagValue);
            rec.setAttribute(SAMTagUtil.getSingleton().SQ, sqTagValue);
        }
        final String e2TagValue = (String)rec.getAttribute(SAMTagUtil.getSingleton().E2);
        if (e2TagValue != null) {
            final byte[] secondaryBases = StringUtil.stringToBytes(e2TagValue);
            reverseComplement(secondaryBases);
            rec.setAttribute(SAMTagUtil.getSingleton().E2, StringUtil.bytesToString(secondaryBases));
        }
        final String u2TagValue = (String)rec.getAttribute(SAMTagUtil.getSingleton().U2);
        if (u2TagValue != null) {
            final byte[] secondaryQualities = StringUtil.stringToBytes(u2TagValue);
            reverseArray(secondaryQualities);
            rec.setAttribute(SAMTagUtil.getSingleton().U2, StringUtil.bytesToString(secondaryQualities));
        }
    }

    /**
     * Reverse the given array in place.
     */
    public static void reverseArray(final byte[] array) {
        final int lastIndex = array.length - 1;
        int i, j;
        for (i=0, j=lastIndex; i<j; ++i, --j) {
            final byte tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }
    }
}
