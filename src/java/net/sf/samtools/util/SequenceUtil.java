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

import java.io.File;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class SequenceUtil {
    /** Byte typed variables for all normal bases. */
    public static final byte a='a', c='c', g='g', t='t', n='n', A='A', C='C', G='G', T='T', N='N';

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

    /** Returns true if the byte is in [acgtACGT]. */
    public static boolean isValidBase(final byte b){
        return b == a || b == A ||
               b == c || b == C ||
               b == g || b == G ||
               b == t || b == T;
    }

    /** Calculates the fraction of bases that are G/C in the sequence. */
    public static double calculateGc(final byte[] bases) {
        int gcs = 0;
        for (int i=0; i<bases.length; ++i) {
            final byte b = bases[i];
            if (b == 'C' || b == 'G' || b == 'c' || b == 'g') ++gcs;
        }

        return gcs / (double) bases.length;
    }

    /**
     * Throws an exception only if both parameters are not null
     * @param s1 a list of sequence headers
     * @param s2 a second list of sequence headers
     */
    public static void assertSequenceListsEqual(final List<SAMSequenceRecord> s1, final List<SAMSequenceRecord> s2) {
        if (s1 != null && s2 != null) {

            if (s1.size() != s2.size()) {
                throw new SequenceListsDifferException(
                    "Sequence dictionaries are not the same size (" + s1.size() + ", " + s2.size() +
                        ")");
            }

            for (int i = 0; i < s1.size(); ++i) {
                if (!s1.get(i).isSameSequence(s2.get(i))) {
                    String s1Attrs = "";
                    for (final java.util.Map.Entry<String, String> entry : s1.get(i)
                        .getAttributes()) {
                        s1Attrs += "/" + entry.getKey() + "=" + entry.getValue();
                    }
                    String s2Attrs = "";
                    for (final java.util.Map.Entry<String, String> entry : s2.get(i)
                        .getAttributes()) {
                        s2Attrs += "/" + entry.getKey() + "=" + entry.getValue();
                    }
                    throw new SequenceListsDifferException(
                        "Sequences at index " + i + " don't match: " +
                            s1.get(i).getSequenceIndex() + "/" + s1.get(i).getSequenceLength() +
                            "/" + s1.get(i).getSequenceName() + s1Attrs + " " +
                            s2.get(i).getSequenceIndex() + "/" + s2.get(i).getSequenceLength() +
                            "/" + s2.get(i).getSequenceName() + s2Attrs);
                }
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
     * Returns true if both parameters are null or equal, otherwise returns false
     */
    public static boolean areSequenceDictionariesEqual(final SAMSequenceDictionary s1, final SAMSequenceDictionary s2) {
        if (s1 == null && s2 == null) return true;
        if (s1 == null || s2 == null) return false;

        try {
            assertSequenceListsEqual(s1.getSequences(), s2.getSequences());
            return true;
        } catch (SequenceListsDifferException e) {
            return false;
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
     * Throws an exception if both parameters are non-null and unequal, including the filenames.
     */
    public static void assertSequenceDictionariesEqual(final SAMSequenceDictionary s1, final SAMSequenceDictionary s2,
                                                       final File f1, final File f2) {
        try {
            assertSequenceDictionariesEqual(s1, s2);
        } catch (SequenceListsDifferException e) {
            throw new SequenceListsDifferException("In files " + f1.getAbsolutePath() + " and " + f2.getAbsolutePath(), e);
        }
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
        return countMismatches(read, referenceBases, 0, false);
    }

    /** Calculates the number of mismatches between the read and the reference sequence provided. */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases, final int referenceOffset) {
        return countMismatches(read, referenceBases, referenceOffset, false);
    }

    /**
     * Calculates the number of mismatches between the read and the reference sequence provided.
     *
     * @param referenceBases Array of ASCII bytes that covers at least the the portion of the reference sequence
     * to which read is aligned from getReferenceStart to getReferenceEnd.
     * @param referenceOffset 0-based offset of the first element of referenceBases relative to the start
     * of that reference sequence.
     * @param bisulfiteSequence If this is true, it is assumed that the reads were bisulfite treated
     *      and C->T on the positive strand and G->A on the negative strand will not be counted
     *      as mismatches.
     */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases, final int referenceOffset,
                                      final boolean bisulfiteSequence) {
        try {
            int mismatches = 0;

            final byte[] readBases = read.getReadBases();

            for (final AlignmentBlock block : read.getAlignmentBlocks()) {
                final int readBlockStart = block.getReadStart() - 1;
                final int referenceBlockStart = block.getReferenceStart() - 1 - referenceOffset;
                final int length = block.getLength();

                for (int i=0; i<length; ++i) {
                    if (!bisulfiteSequence) {
                        if (!basesEqual(readBases[readBlockStart+i], referenceBases[referenceBlockStart+i])) {
                            ++mismatches;
                        }
                    }
                    else {
                        if (!bisulfiteBasesEqual(read.getReadNegativeStrandFlag(), readBases[readBlockStart+i],
                                referenceBases[referenceBlockStart+i])) {
                            ++mismatches;
                        }
                    }
                }
            }
            return mismatches;
        } catch (Exception e) {
            throw new SAMException("Exception counting mismatches for read " + read, e);
        }
    }

    /**
     * Calculates the number of mismatches between the read and the reference sequence provided.
     *
     * @param referenceBases Array of ASCII bytes that covers at least the the portion of the reference sequence
     * to which read is aligned from getReferenceStart to getReferenceEnd.
     * @param bisulfiteSequence If this is true, it is assumed that the reads were bisulfite treated
     *      and C->T on the positive strand and G->A on the negative strand will not be counted
     *      as mismatches.
     */
    public static int countMismatches(final SAMRecord read, final byte[] referenceBases, final boolean bisulfiteSequence) {
        return countMismatches(read, referenceBases, 0, bisulfiteSequence);
    }

    /**
     * Sadly, this is a duplicate of the method above, except that it takes char[] for referenceBases rather
     * than byte[].  This is because GATK needs it this way.
     *
     * TODO: Remove this method when GATK map method is changed to take refseq as byte[].
     */
    private static int countMismatches(final SAMRecord read, final char[] referenceBases, final int referenceOffset) {
        int mismatches = 0;

        final byte[] readBases = read.getReadBases();

        for (final AlignmentBlock block : read.getAlignmentBlocks()) {
            final int readBlockStart = block.getReadStart() - 1;
            final int referenceBlockStart = block.getReferenceStart() - 1 - referenceOffset;
            final int length = block.getLength();

            for (int i=0; i<length; ++i) {
                if (!basesEqual(readBases[readBlockStart+i], StringUtil.charToByte(referenceBases[referenceBlockStart+i]))) {
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
        return sumQualitiesOfMismatches(read, referenceBases, 0, false);
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
        return sumQualitiesOfMismatches(read, referenceBases, referenceOffset, false);
    }

    /**
     * Calculates the sum of qualities for mismatched bases in the read.
     * @param referenceBases Array of ASCII bytes that covers at least the the portion of the reference sequence
     * to which read is aligned from getReferenceStart to getReferenceEnd.
     * @param referenceOffset 0-based offset of the first element of referenceBases relative to the start
     * of that reference sequence. 
     * @param bisulfiteSequence If this is true, it is assumed that the reads were bisulfite treated
     *      and C->T on the positive strand and G->A on the negative strand will not be counted
     *      as mismatches.
     */
    public static int sumQualitiesOfMismatches(final SAMRecord read, final byte[] referenceBases,
                                               final int referenceOffset, final boolean bisulfiteSequence) {
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
                if (!bisulfiteSequence) {
                    if (!basesEqual(readBases[readBlockStart+i], referenceBases[referenceBlockStart+i])) {
                        qualities += readQualities[readBlockStart+i];
                    }

                }
                else {
                    if (!bisulfiteBasesEqual(read.getReadNegativeStrandFlag(), readBases[readBlockStart+i],
                            referenceBases[referenceBlockStart+i])) {
                        qualities += readQualities[readBlockStart+i];
                    }
                }
            }
        }

        return qualities;
    }

    /**
     * Sadly, this is a duplicate of the method above, except that it takes char[] for referenceBases rather
     * than byte[].  This is because GATK needs it this way.
     *
     * TODO: Remove this method when GATK map method is changed to take refseq as byte[].
     */
    public static int sumQualitiesOfMismatches(final SAMRecord read, final char[] referenceBases,
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
                if (!basesEqual(readBases[readBlockStart+i], StringUtil.charToByte(referenceBases[referenceBlockStart+i]))) {
                    qualities += readQualities[readBlockStart+i];
                }
            }
        }

        return qualities;
    }

    public static int countInsertedBases(final Cigar cigar) {
        int ret = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.INSERTION) ret += element.getLength();
        }
        return ret;
    }

    public static int countDeletedBases(final Cigar cigar) {
        int ret = 0;
        for (final CigarElement element : cigar.getCigarElements()) {
            if (element.getOperator() == CigarOperator.DELETION) ret += element.getLength();
        }
        return ret;
    }

    public static int countInsertedBases(final SAMRecord read) {
        return countInsertedBases(read.getCigar());
    }

    public static int countDeletedBases(final SAMRecord read) {
        return countDeletedBases(read.getCigar());
    }

    /**
     * Calculates the for the predefined NM tag from the SAM spec. To the result of
     * countMismatches() it adds 1 for each indel.
     */
    public static int calculateSamNmTag(final SAMRecord read, final byte[] referenceBases) {
        return calculateSamNmTag(read, referenceBases, 0, false);
    }

    /**
     * Calculates the for the predefined NM tag from the SAM spec. To the result of
     * countMismatches() it adds 1 for each indel.

     * @param referenceOffset 0-based offset of the first element of referenceBases relative to the start
     * of that reference sequence.
     */
    public static int calculateSamNmTag(final SAMRecord read, final byte[] referenceBases,
                                        final int referenceOffset) {
        return calculateSamNmTag(read, referenceBases, referenceOffset, false);
    }

    /**
     * Calculates the for the predefined NM tag from the SAM spec. To the result of
     * countMismatches() it adds 1 for each indel.

     * @param referenceOffset 0-based offset of the first element of referenceBases relative to the start
     * of that reference sequence.
     * @param bisulfiteSequence If this is true, it is assumed that the reads were bisulfite treated
     *      and C->T on the positive strand and G->A on the negative strand will not be counted
     *      as mismatches.
     */
    public static int calculateSamNmTag(final SAMRecord read, final byte[] referenceBases,
                                        final int referenceOffset, final boolean bisulfiteSequence) {
        int samNm = countMismatches(read, referenceBases, referenceOffset, bisulfiteSequence);
        for (final CigarElement el : read.getCigar().getCigarElements()) {
            if (el.getOperator() == CigarOperator.INSERTION || el.getOperator() == CigarOperator.DELETION) {
                samNm += el.getLength();
            }
        }
        return samNm;
    }

    /**
     * Sadly, this is a duplicate of the method above, except that it takes char[] for referenceBases rather
     * than byte[].  This is because GATK needs it this way.
     *
     * TODO: Remove this method when GATK map method is changed to take refseq as byte[].
     */
    public static int calculateSamNmTag(final SAMRecord read, final char[] referenceBases,
                                               final int referenceOffset) {
        int samNm = countMismatches(read, referenceBases, referenceOffset);
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

    
    /** Reverses the quals in place. */
    public static void reverseQualities(final byte[] quals) {
        final int lastIndex = quals.length - 1;

        int i, j;
        for (i=0, j=lastIndex; i<j; ++i, --j) {
            final byte tmp = quals[i];
            quals[i] = quals[j];
            quals[j] = tmp;
        }
    }
    
    /** Returns true if the bases are equal OR if the mismatch cannot be accounted for by
     * bisfulite treatment.  C->T on the positive strand and G->A on the negative strand
     * do not count as mismatches */
    public static boolean bisulfiteBasesEqual(final boolean negativeStrand, final byte read, final byte reference) {

        if (basesEqual(read, reference)) {
            return true;
        }

        if (negativeStrand) {
            if (basesEqual(reference, (byte)'G') && basesEqual(read, (byte)'A')) {
                return true;
            }
        }
        else {
            if (basesEqual(reference, (byte)'C') && basesEqual(read, (byte)'T')) {
                return true;
            }
        }
        return false;
    }


    /*
     * Regexp for MD string.
     *
     * \G = end of previous match.
     * (?:[0-9]+) non-capturing (why non-capturing?) group of digits.  For this number of bases read matches reference.
     *  - or -
     * Single reference base for case in which reference differs from read.
     *  - or -
     * ^one or more reference bases that are deleted in read.
     *
     */
    static final Pattern mdPat = Pattern.compile("\\G(?:([0-9]+)|([ACTGNactgn])|(\\^[ACTGNactgn]+))");

    /**
     * Produce reference bases from an aligned SAMRecord with MD string and Cigar.
     * @param rec Must contain non-empty CIGAR and MD attribute.
     * @param includeReferenceBasesForDeletions If true, include reference bases that are deleted in the read.
     * This will make the returned array not line up with the read if there are deletions.
     * @return References bases corresponding to the read.  If there is an insertion in the read, reference contains
     * '-'.  If the read is soft-clipped, reference contains '0'.  If there is a skipped region and
     * includeReferenceBasesForDeletions==true, reference will have Ns for the skipped region.
     */
    public static byte[] makeReferenceFromAlignment(final SAMRecord rec, final boolean includeReferenceBasesForDeletions) {
        final String md = rec.getStringAttribute(SAMTag.MD.name());
        if (md == null) {
            throw new SAMException("Cannot create reference from SAMRecord with no MD tag, read: " + rec.getReadName());
        }
        // Not sure how long output will be, but it will be no longer than this.
        int maxOutputLength = 0;
        final Cigar cigar = rec.getCigar();
        if (cigar == null) {
            throw new SAMException("Cannot create reference from SAMRecord with no CIGAR, read: " + rec.getReadName());
        }
        for (final CigarElement cigarElement : cigar.getCigarElements()) {
            maxOutputLength += cigarElement.getLength();
        }
        final byte[] ret = new byte[maxOutputLength];
        int outIndex = 0;

        final Matcher match = mdPat.matcher(md);
        int curSeqPos = 0;

        int savedBases = 0;
        final byte[] seq = rec.getReadBases();
        for (final CigarElement cigEl : cigar.getCigarElements())
        {
            final int cigElLen = cigEl.getLength();
            final CigarOperator cigElOp = cigEl.getOperator();


            if (cigElOp == CigarOperator.SKIPPED_REGION) {
                // We've decided that MD tag will not contain bases for skipped regions, as they
                // could be megabases long, so just put N in there if caller wants reference bases,
                // otherwise ignore skipped regions.
                if (includeReferenceBasesForDeletions) {
                    for (int i = 0; i < cigElLen; ++i) {
                    ret[outIndex++] = N;
                    }
                }
            }
            // If it consumes reference bases, it's either a match or a deletion in the sequence
            // read.  Either way, we're going to need to parse through the MD.
            else if (cigElOp.consumesReferenceBases()) {
                // We have a match region, go through the MD
                int basesMatched = 0;

                // Do we have any saved matched bases?
                while ((savedBases>0) && (basesMatched < cigElLen))
                {
                    ret[outIndex++] = seq[curSeqPos++];
                    savedBases--;
                    basesMatched++;
                }

                while (basesMatched < cigElLen)
                {
                    boolean matched = match.find();
                    if (matched)
                    {
                        String mg;
                        if ( ((mg = match.group(1)) !=null) && (mg.length() > 0) )
                        {
                            // It's a number , meaning a series of matches
                            final int num = Integer.parseInt(mg);
                            for (int i = 0; i < num; i++)
                            {
                                if (basesMatched<cigElLen)
                                {
                                    ret[outIndex++] = seq[curSeqPos++];
                                }
                                else
                                {
                                    savedBases++;
                                }
                                basesMatched++;
                            }
                        }

                        else if ( ((mg = match.group(2)) !=null) && (mg.length() > 0) )
                        {
                            // It's a single nucleotide, meaning a mismatch
                            if (basesMatched<cigElLen)
                            {
                                ret[outIndex++] = StringUtil.charToByte(mg.charAt(0));
                                curSeqPos++;
                            }
                            else
                            {
                                throw new IllegalStateException("Should never happen.");
                            }
                            basesMatched++;
                        }
                        else if ( ((mg = match.group(3)) !=null) && (mg.length() > 0) )
                        {
                            // It's a deletion, starting with a caret
                            // don't include caret
                            if (includeReferenceBasesForDeletions) {
                                final byte[] deletedBases = StringUtil.stringToBytes(mg);
                                System.arraycopy(deletedBases, 1, ret, outIndex, deletedBases.length - 1);
                                outIndex += deletedBases.length - 1;
                            }
                            basesMatched += mg.length() - 1;

                            // Check just to make sure.
                            if (basesMatched != cigElLen)
                            {
                                throw new SAMException("Got a deletion in CIGAR (" + cigar + ", deletion " + cigElLen +
                                        " length) with an unequal ref insertion in MD (" + md + ", md " + basesMatched + " length");
                            }
                            if (cigElOp != CigarOperator.DELETION)
                            {
                                throw new SAMException ("Got an insertion in MD ("+md+") without a corresponding deletion in cigar ("+cigar+")");
                            }

                        }
                        else
                        {
                            matched = false;
                        }
                    }

                    if (!matched)
                    {
                        throw new SAMException("Illegal MD pattern: " + md + " for read " + rec.getReadName() +
                                " with CIGAR " + rec.getCigarString());
                    }
                }

            }
            else if (cigElOp.consumesReadBases())
            {
                // We have an insertion in read
                for (int i = 0; i < cigElLen; i++)
                {
                    final char c = (cigElOp == CigarOperator.SOFT_CLIP) ? '0' : '-';
                    ret[outIndex++] =  StringUtil.charToByte(c);
                    curSeqPos++;
                }
            }
            else
            {
                // It's an op that consumes neither read nor reference bases.  Do we just ignore??
            }

        }
        if (outIndex < ret.length) {
            final byte[] shorter = new byte[outIndex];
            System.arraycopy(ret, 0, shorter, 0, outIndex);
            return shorter;
        }
        return ret;
    }
}
