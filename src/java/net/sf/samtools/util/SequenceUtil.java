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
package net.sf.samtools.util;

import net.sf.samtools.*;

import java.util.List;

public class SequenceUtil {
    /** Byte typed variables for all normal bases. */
    public static final byte a='a', c='c', g='g', t='t', A='A', C='C', G='G', T='T';

    /**
     * Calculate the reverse complement of the specified sequence
     * (Stolen from Reseq)
     *
     * @param sequenceData
     * @return reverse complement
     */
    public static String reverseComplement(final String sequenceData) {
        final byte[] bases = net.sf.samtools.util.StringUtil.stringToBytes(sequenceData);
        reverseComplement(bases);
        return net.sf.samtools.util.StringUtil.bytesToString(bases);
    }

    /** Attempts to efficiently compare two bases stored as bytes for equality. */
    public static boolean basesEqual(byte lhs, byte rhs) {
        if (lhs == rhs) return true;
        else {
            if (lhs > 90) lhs -= 32;
            if (rhs > 90) rhs -= 32;
        }

        return lhs == rhs;
    }
    
    /**
     * returns true if the value of base represents a no call
     */
    public static boolean isNoCall(final byte base) {
        return base == 'N' || base == 'n' || base == '.';
    }

    /**
     * Throws an exception if both parameters are not null sequenceListsEqual returns false
     * @param s1 a list of sequence headers
     * @param s2 a second list of sequence headers
     */
    public static void assertSequenceListsEqual(final List<SAMSequenceRecord> s1, final List<SAMSequenceRecord> s2) {
        if (s1 == null || s2 == null) return;
        
        if (s1.size() != s2.size()) {
            throw new SequenceListsDifferException("Sequence dictionaries are not the same size (" +
                    s1.size() + ", " + s2.size() + ")");
        }
        for (int i = 0; i < s1.size(); ++i) {
            if (!s1.get(i).isSameSequence(s2.get(i))) {
                String s1Attrs = "";
                for (final java.util.Map.Entry<String, Object> entry : s1.get(i).getAttributes()) {
                    s1Attrs += "/" + entry.getKey() + "=" + entry.getValue();
                }
                String s2Attrs = "";
                for (final java.util.Map.Entry<String, Object> entry : s2.get(i).getAttributes()) {
                    s2Attrs += "/" + entry.getKey() + "=" + entry.getValue();
                }
                throw new SequenceListsDifferException("Sequences at index " + i + " don't match: " +
                    s1.get(i).getSequenceIndex() + "/" + s1.get(i).getSequenceLength() + "/" + s1.get(i).getSequenceName() + s1Attrs +
                    " " + s2.get(i).getSequenceIndex() + "/" + s2.get(i).getSequenceLength() + "/" + s2.get(i).getSequenceName() + s2Attrs);
            }
        }
    }

    public static class SequenceListsDifferException extends SAMException {
        public SequenceListsDifferException() {
        }

        public SequenceListsDifferException(final String s) {
            super(s);
        }

        public SequenceListsDifferException(final String s, final Throwable throwable) {
            super(s, throwable);
        }

        public SequenceListsDifferException(final Throwable throwable) {
            super(throwable);
        }
    }

    /**
     * Throws an exception if both parameters are non-null and unequal.
     */
    public static void assertSequenceDictionariesEqual(final SAMSequenceDictionary s1, final SAMSequenceDictionary s2) {
        if (s1 == null || s2 == null) return;
        assertSequenceListsEqual(s1.getSequences(), s2.getSequences());
    }

    /**
     * Create a simple ungapped cigar string, which might have soft clipping at either end
     * @param alignmentStart raw aligment start, which may result in read hanging off beginning or end of read
     * @return cigar string that may have S operator at beginning or end, and has M operator for the rest of the read
     */
    public static String makeCigarStringWithPossibleClipping(final int alignmentStart, final int readLength, final int referenceSequenceLength) {
        int start = alignmentStart;
        int leftSoftClip = 0;
        if (start < 1) {
            leftSoftClip = 1 - start;
            start = 1;
        }
        int rightSoftClip = 0;
        if (alignmentStart + readLength > referenceSequenceLength + 1) {
            rightSoftClip = alignmentStart + readLength - referenceSequenceLength - 1;
        }
        // CIGAR is trivial because there are no indels or clipping in Gerald
        final int matchLength = readLength - leftSoftClip - rightSoftClip;
        if (matchLength < 1) {
            throw new SAMException("Unexpected cigar string with no M op for read.");
        }
        return makeSoftClipCigar(leftSoftClip) + Integer.toString(matchLength) + "M" + makeSoftClipCigar(rightSoftClip);
    }

    /**
     * Create a cigar string for a gapped alignment, which may have soft clipping at either end
     * @param alignmentStart raw alignment start, which may result in read hanging off beginning or end of read
     * @param readLength
     * @param referenceSequenceLength
     * @param indelPosition number of matching bases before indel.  Must be > 0
     * @param indelLength length of indel.  Positive for insertion, negative for deletion.
     * @return cigar string that may have S operator at beginning or end, has one or two M operators, and an I or a D.
     */
    public static String makeCigarStringWithIndelPossibleClipping(final int alignmentStart,
                                                                  final int readLength,
                                                                  final int referenceSequenceLength,
                                                                  final int indelPosition,
                                                                  final int indelLength) {
        int start = alignmentStart;
        int leftSoftClip = 0;
        if (start < 1) {
            leftSoftClip = 1 - start;
            start = 1;
        }
        int rightSoftClip = 0;
        final int alignmentEnd = alignmentStart + readLength - indelLength;
        if (alignmentEnd > referenceSequenceLength + 1) {
            rightSoftClip = alignmentEnd - referenceSequenceLength - 1;
        }
        if (leftSoftClip >= indelPosition) {
            throw new IllegalStateException("Soft clipping entire pre-indel match. leftSoftClip: " + leftSoftClip +
            "; indelPosition: " + indelPosition);
        }
        // CIGAR is trivial because there are no indels or clipping in Gerald
        final int firstMatchLength = indelPosition - leftSoftClip;
        final int secondMatchLength = readLength - indelPosition - (indelLength > 0? indelLength: 0) - rightSoftClip;
        if (secondMatchLength < 1) {
            throw new SAMException("Unexpected cigar string with no M op for read.");
        }
        return makeSoftClipCigar(leftSoftClip) + Integer.toString(firstMatchLength) + "M" +
                Math.abs(indelLength) + (indelLength > 0? "I": "D") +
                Integer.toString(secondMatchLength) + "M" +
                makeSoftClipCigar(rightSoftClip);
    }

    public static String makeSoftClipCigar(final int clipLength) {
        if (clipLength == 0) {
            return "";
        }
        return Integer.toString(clipLength) + "S";
    }

    /** Calculates the number of mismatches between the read and the reference sequence provided. */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases) {
        int mismatches = 0;

        final byte[] readBases = read.getReadBases();

        for (final AlignmentBlock block : read.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - 1;
            final int length = block.getLength();

            for (int i=0; i<length; ++i) {
                if (!basesEqual(readBases[readBlockStart+i], referenceBases[referenceBlockStart+i])) {
                    ++mismatches;
                }
            }
        }
        return mismatches;
    }

    /**
     * Calculates the sum of qualities for mismatched bases in the read.
     * @param referenceBases Array of ASCII bytes in which the 0th position in the array corresponds
     * to the first element of the reference sequence to which read is aligned. 
     */
    public static int sumQualitiesOfMismatches(final SAMRecord read, final byte[] referenceBases) {
        return sumQualitiesOfMismatches(read, referenceBases, 0);
    }

    /**
     * Calculates the sum of qualities for mismatched bases in the read.
     * @param referenceBases Array of ASCII bytes that covers at least the the portion of the reference sequence
     * to which read is aligned from getReferenceStart to getReferenceEnd.
     * @param referenceOffset 0-based offset of the first element of referenceBases relative to the start
     * of that reference sequence. 
     */
    public static int sumQualitiesOfMismatches(final SAMRecord read, final byte[] referenceBases,
                                               final int referenceOffset) {
        int qualities = 0;

        final byte[] readBases = read.getReadBases();
        final byte[] readQualities = read.getBaseQualities();

        if (read.getAlignmentStart() <= referenceOffset) {
            throw new IllegalArgumentException("read.getAlignmentStart(" + read.getAlignmentStart() +
                    ") <= referenceOffset(" + referenceOffset + ")");
        }

        for (final AlignmentBlock block : read.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - 1 - referenceOffset;
            final int length = block.getLength();

            for (int i=0; i<length; ++i) {
                if (!basesEqual(readBases[readBlockStart+i], referenceBases[referenceBlockStart+i])) {
                    qualities += readQualities[readBlockStart+i];
                }
            }
        }

        return qualities;
    }

    /**
     * Calculates the for the predefined NM tag from the SAM spec. To the result of
     * countMismatches() it adds 1 for each indel.
     */
    public static int calculateSamNmTag(final SAMRecord read, final byte[] referenceBases) {
        int samNm = countMismatches(read, referenceBases);
        for (final CigarElement el : read.getCigar().getCigarElements()) {
            if (el.getOperator() == CigarOperator.INSERTION || el.getOperator() == CigarOperator.DELETION) {
                samNm += el.getLength();
            }
        }
        return samNm;
    }

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
}