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

import net.sf.samtools.util.CoordMath;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.util.SequenceUtil;

import java.util.List;

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

    /** Mishali
     *
     * @return Cigar object for the mate read, or null if there is none.
     */
    public static Cigar getMateCigarMishali(final SAMRecord rec) {
        if (!rec.getMateUnmappedFlag()) {
            final String cigarString = (String)rec.getAttribute(SAMTagUtil.getSingleton().MC);
            if (cigarString == null)
                throw new SAMException("Tag MC for mate cigar not found.");
            else {
                final Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
                if (rec.getValidationStringency() != SAMFileReader.ValidationStringency.SILENT && !rec.getMateUnmappedFlag()) {
                    //validateCigar function for mate needed
                    //		    SAMUtils.processValidationErrors(validateCigar(-1L), -1L, getValidationStringency());
                }
                return cigar;
            }
        } else
            throw new RuntimeException("getMateCigar called on an unmapped mate.");
    }

    /** ggrant
     * Replacing Mishali's method with one that doesn't throw Exceptions if MC not found..
     * TODO - verify that this doesn't break her methods
     *
     * @return Cigar object for the mate read, or null if there is none.
     */
    public static Cigar getMateCigar(final SAMRecord rec) {
        return rec.getMateCigar();
    }


    /** Mishali
     * @return 1-based inclusive rightmost position of the clipped mate sequence, or 0 read if unmapped.
     */
    public static int getMateAlignmentEnd(final SAMRecord rec) {
        if (rec.getMateUnmappedFlag()) {
            throw new RuntimeException("getMateAlignmentEnd called on an unmapped mate.");
        }
        //int mateAlignmentEnd = rec.getMateAlignmentStart();
        //	mateAlignmentEnd = mateAlignmentEnd + getMateCigar(rec).getReferenceLength() - 1;
        return CoordMath.getEnd(rec.getMateAlignmentStart(), getMateCigar(rec).getReferenceLength());
    }

    /** Mishali
     * @return the mate alignment start (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment start of 100 but the first 4 bases were clipped (hard or soft clipped)
     * then this method will return 96.
     *
     * Invalid to call on an unmapped read.
     */
    public static int getMateUnclippedStart(final SAMRecord rec) {
        if (rec.getMateUnmappedFlag())
            throw new RuntimeException("getMateUnclippedStart called on an unmapped mate.");
        int pos = rec.getMateAlignmentStart();

        for (final CigarElement cig : getMateCigar(rec).getCigarElements()) {
            final CigarOperator op = cig.getOperator();
            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                pos -= cig.getLength();
            }
            else {
                break;
            }
        }

        return pos;
    }

    /** Mishali
     * @return the alignment end (1-based, inclusive) adjusted for clipped bases.  For example if the read
     * has an alignment end of 100 but the last 7 bases were clipped (hard or soft clipped)
     * then this method will return 107.
     *
     * Invalid to call on an unmapped read.
     */
    public static int getMateUnclippedEnd(final SAMRecord rec) {
        if (rec.getMateUnmappedFlag())
            throw new RuntimeException("getMateUnclippedEnd called on an unmapped mate.");

        int pos = getMateAlignmentEnd(rec);
        final List<CigarElement> cigs = getMateCigar(rec).getCigarElements();
        for (int i=cigs.size() - 1; i>=0; --i) {
            final CigarElement cig = cigs.get(i);
            final CigarOperator op = cig.getOperator();

            if (op == CigarOperator.SOFT_CLIP || op == CigarOperator.HARD_CLIP) {
                pos += cig.getLength();
            }
            else {
                break;
            }
        }

        return pos;
    }

}
