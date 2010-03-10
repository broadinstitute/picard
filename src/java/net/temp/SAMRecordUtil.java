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
import net.sf.samtools.util.SequenceUtil;

/**
 * @author alecw@broadinstitute.org
 */
public class SAMRecordUtil {

    /** List of String tags that must be reversed if present when a SAMRecord is reverseComplemented */
    private static final short[] STRING_TAGS_TO_REVERSE = {
            SAMTagUtil.getSingleton().U2,
            SAMTagUtil.getSingleton().OQ
    };

    /**
     * Reverse-complement all known sequence and base quality attributes of the SAMRecord.
     */
    public static void reverseComplement(final SAMRecord rec) {
        final byte[] readBases = rec.getReadBases();
        SequenceUtil.reverseComplement(readBases);
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
            SequenceUtil.reverseComplement(secondaryBases);
            rec.setAttribute(SAMTagUtil.getSingleton().E2, StringUtil.bytesToString(secondaryBases));
        }
        for (final short stringTag : STRING_TAGS_TO_REVERSE) {
            final String value = (String)rec.getAttribute(stringTag);
            if (value != null) {
                rec.setAttribute(stringTag, StringUtil.reverseString(value));
            }
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
