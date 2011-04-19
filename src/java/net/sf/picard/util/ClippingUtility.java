/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

package net.sf.picard.util;

import net.sf.picard.util.IlluminaUtil;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StringUtil;
import net.sf.samtools.SAMRecord;
import net.sf.picard.sam.ReservedTagConstants;

/**
 * Utilities to clip the adapater sequence from an illumina read
 *
 * @author Tim Fennell (adapted by mborkan@broadinstitute.org)
 */
public class ClippingUtility {

    /**
     * The default value used for the minimum number of contiguous bases to match against.
     */
    public static final int MIN_MATCH_BASES = 12;
    /**
     * The default value used for the minimum number of contiguous bases to match against in a paired end read
     */
    public static final int MIN_MATCH_PE_BASES = 6;

    /**
     * The default value used for the maximum error rate when matching read bases to clippable sequence.
     */
    public static final double MAX_ERROR_RATE = 0.10;
    /**
     * The default value used for the maximum error rate when matching paired end read bases to clippable sequence.
     */
    public static final double MAX_PE_ERROR_RATE = 0.10;

    /**
     * The value returned by methods returning int when no match is found.
     */
    public static final int NO_MATCH = -1;

    /**
     * Invokes adapterTrimIlluminRead with default parameters for a single read.
     * If the read is a negative strand, its bases will be reverse complemented
     * Simpler, more common of two overloads.
     *
     * @param read    SAM/BAM read to trim
     * @param adapter which adapters to use (indexed, paired_end, or single_end)
     */
    public static void adapterTrimIlluminaSingleRead(final SAMRecord read, final IlluminaUtil.AdapterPair adapter) {
        adapterTrimIlluminaSingleRead(read, adapter, MIN_MATCH_BASES, MAX_ERROR_RATE);
    }

    /**
     * Invokes adapterTrimIlluminRead with explicit matching thresholds for a single read.
     * If the read is a negative strand, a copy of its bases will be reverse complemented
     * More general from of two overloads.
     *
     * @param read          SAM/BAM read to trim
     * @param adapter       which adapters to use (indexed, paired_end, or single_end)
     * @param minMatchBases minimum number of contiguous bases to match against in a read
     * @param maxErrorRate  maximum error rate when matching read bases
     */
    public static void adapterTrimIlluminaSingleRead(final SAMRecord read, final IlluminaUtil.AdapterPair adapter,
        final int minMatchBases, final double maxErrorRate) {
        final byte[] bases = read.getReadBases();
        int indexOfAdapterSequence = findIndex(bases, read.getReadNegativeStrandFlag(),
                adapter.get3PrimeAdapterBytes(), minMatchBases, maxErrorRate);
        if (indexOfAdapterSequence != NO_MATCH) {
            // Convert to a one-based index for storage on the record.
            read.setAttribute(ReservedTagConstants.XT, indexOfAdapterSequence + 1);
        }
    }

    /**
     * Invokes adapterTrimIlluminaRead with default less stringent parameters for a pair of reads.
     * If the read is a negative strand, its bases will be reverse complemented
     * Simpler, more common of two overloads.
     * Returns a warning string when the trim positions found differed for each read.
     *
     * @param read1    first read of the pair
     * @param read2    second read of the pair
     * @param adapters which adapters to use (indexed, paired_end, or single_end)
     * @return warning information about the trimming for logging
     */
    public static String adapterTrimIlluminaPairedReads(final SAMRecord read1, final SAMRecord read2, final IlluminaUtil.AdapterPair adapters) {
        return adapterTrimIlluminaPairedReads(read1, read2, adapters, MIN_MATCH_PE_BASES, MAX_PE_ERROR_RATE);
    }

    /**
     * Invokes adapterTrimIlluminaRead with explicit parameters for a pair of reads.
     * More general form of two overloads.
     * Returns a warning string when the trim positions found differed for each read.
     *
     * @param read1         first read of the pair.
     * If read1 is a negative strand, a copy of its bases will be reverse complemented.
     * @param read2         second read of the pair.
     * If read2 is a negative strand, a copy of its bases will be reverse complemented
     * @param adapters      which adapters to use (indexed, paired_end, or single_end)
     * @param minMatchBases minimum number of contiguous bases to match against in a read
     * @param maxErrorRate  maximum error rate when matching read bases
     * @return warning      information about the trimming for logging
     */
    public static String adapterTrimIlluminaPairedReads(final SAMRecord read1, final SAMRecord read2,
        final IlluminaUtil.AdapterPair adapters, final int minMatchBases, final double maxErrorRate) {

        String warnString = null;

        final byte[] bases1 = read1.getReadBases();
        final int index1 = findIndex(bases1, read1.getReadNegativeStrandFlag(),
                adapters.get3PrimeAdapterBytes(), minMatchBases, maxErrorRate);
        final byte[] bases2 = read2.getReadBases();
        final int index2 = findIndex(bases2, read2.getReadNegativeStrandFlag(),
                adapters.get5PrimeAdapterBytesInReadOrder(), minMatchBases, maxErrorRate);

        if (index1 == index2) {
            if (index1 == NO_MATCH) {
                // neither matched, nothing to do
            } else {
                // both match at same place, trim there
                read1.setAttribute(ReservedTagConstants.XT, index1 + 1);
                read2.setAttribute(ReservedTagConstants.XT, index2 + 1);
            }
        } else if (index1 == NO_MATCH || index2 == NO_MATCH) {
            // one of them matched, but the other didn't.
            // Trim them both at the matching point (if that makes sense)
            // and if one side still matches at default/stricter settings
            int stricterMinMatchBases = 2 * minMatchBases;
            warnString = oneSidedMatch(read1, read2, bases1, bases2,
                index1, index2, adapters, stricterMinMatchBases, maxErrorRate);

        } else {
            // both matched at different positions. Do nothing
            warnString = "Adapters mismatch at position " + index1 + " " + StringUtil.bytesToString(bases1) +
                " " + read1 + " and reverse " + index2 + " " + StringUtil.bytesToString(bases2);
        }
        return warnString;
    }

    private static int findIndex(byte[] bases, boolean isNegativeStrand,
                                 byte[] adapterBytes, int minMatchBases, double maxErrorRate) {
        byte[] copiedBases = null;
        if (isNegativeStrand){
            copiedBases = new byte[bases.length];
            System.arraycopy(bases, 0, copiedBases, 0, bases.length);
            SequenceUtil.reverseComplement(copiedBases);
        } else {
            copiedBases = bases;
        }
        return findIndexOfClipSequence(copiedBases, adapterBytes, minMatchBases, maxErrorRate);

    }

    // If the adapter still matches using normal (stricter) thresholds, trims both at the same position
    // pre-condition - the adapter matches only one set of bases using the relaxed thresholds.
    private static String oneSidedMatch(SAMRecord read1, SAMRecord read2,
                                        byte[] bases1, byte[] bases2,
                                        int index1, int index2,
                                        IlluminaUtil.AdapterPair adapters,
                                        final int stricterMinMatchBases,
                                        final double maxErrorRate) {
        String warnString;
        String willTrim = "won't"; // 'won't' implies nonsense lengths
        // check if still matches at normal/stricter thresholds

        // if the error_rate is the same for the strict and the less-strict,
        // then we don't need to match again.
        // Just check if the readLength - index of successful match > stricterMinMatchBases
        boolean needsRematch = true;
        boolean stillMatches = false;
        if (MAX_ERROR_RATE == maxErrorRate){
            needsRematch = false;  // no need to rematch
            int successIndex = index1 == NO_MATCH ? index2 : index1;
            SAMRecord successfulRead = index1 == NO_MATCH ? read2 : read1;
            if (successfulRead.getReadLength() - successIndex >= stricterMinMatchBases){
                stillMatches = true;
            }
        }

        if ( (!needsRematch && !stillMatches) ||
             needsRematch && indexOfStrictMatchingBase(index1, bases1, bases2, adapters, stricterMinMatchBases) == NO_MATCH) {
            willTrim = "will not"; // 'will not' implies rematch failed
        } else {
            // trim both sides
            int trimIndex = index1 == NO_MATCH ? index2 : index1;
            if (bases1.length > trimIndex) {
                read1.setAttribute(ReservedTagConstants.XT, trimIndex + 1);
                willTrim = "will only trim at position " + trimIndex;
                // 'will only' implies just 1st read
            }
            if (bases2.length > trimIndex) {
                read2.setAttribute(ReservedTagConstants.XT, trimIndex + 1);
                willTrim = "will trim at position " + trimIndex;
                // will implies just 2nd read or both
            }
        }
        warnString = "No adapter match " + willTrim +
                " trim in paired read of length " + bases1.length + " and length " + bases2.length + " " +
                (index1 == NO_MATCH ? read2 + " " + StringUtil.bytesToString(bases2)
                        : " reverse " + read1 + " " + StringUtil.bytesToString(bases1)) +
                " after strict check using minMaxBases=" + stricterMinMatchBases;
        return warnString;
    }

    // does the adapter still match even at normal thresholds
    // preCondition - the adapter matches only one set of bases at the relaxed thresholds
    private static int indexOfStrictMatchingBase(int index1, byte[] bases1, byte[] bases2,
                                                 IlluminaUtil.AdapterPair adapters,
                                                 int stricterMinMatchBases) {
        byte[] basesToCheck, adapterToCheck;
        if (index1 == NO_MATCH) {
            basesToCheck = bases2;
            adapterToCheck = adapters.get5PrimeAdapterBytesInReadOrder();
        } else {
            basesToCheck = bases1;
            adapterToCheck = adapters.get3PrimeAdapterBytesInReadOrder();
        }
        int newIndex = findIndexOfClipSequence(basesToCheck, adapterToCheck,
                stricterMinMatchBases, MAX_ERROR_RATE);
        return newIndex;
    }

    /**
     * Finds the first index of the "toClip" sequence in the "read" sequence requiring at least minMatch
     * bases of pairwise alignment with a maximum number of errors dictacted by maxErrorRate.
     */
   public static int findIndexOfClipSequence(final byte[] read, final byte[] toClip, final int minMatch, final double maxErrorRate) {
        // If the read's too short we can't possibly match it
        if (read == null || read.length < minMatch) return NO_MATCH;
        final int minClipPosition = Math.max (0, read.length - toClip.length);

        // Walk backwards down the read looking for the sequence
        READ_LOOP:
        for (int start = read.length - minMatch; start > minClipPosition -1; --start) {
            final int length = Math.min(read.length - start, toClip.length);
            final int mismatchesAllowed = (int) (length * maxErrorRate);
            int mismatches = 0;

            for (int i = 0; i < length; ++i) {
                if (!SequenceUtil.isNoCall(toClip[i]) && !SequenceUtil.basesEqual(toClip[i], read[start + i])) {
                    if (++mismatches > mismatchesAllowed) continue READ_LOOP;
                }
            }

            // If we got this far without breaking out, then it matches
            return start;
        }

        return NO_MATCH;
    }
}