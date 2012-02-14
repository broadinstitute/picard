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

package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.util.CoordMath;

import java.util.Iterator;
import java.util.List;

/**
 * Utility methods for pairs of SAMRecords
 */
public class SamPairUtil {

    /**
    * The possible orientations of paired reads.
    *
    * F = mapped to forward strand
    * R = mapped to reverse strand
    *
    * FR means the read that's mapped to the forward strand comes before the
    * read mapped to the reverse strand when their 5'-end coordinates are
    * compared.
    */
    public static enum PairOrientation
    {
        FR,     // ( 5' --F-->       <--R-- 5'  )  - aka. innie
        RF,     // (   <--R-- 5'   5' --F-->    )  - aka. outie
        TANDEM; // ( 5' --F-->   5' --F-->  or  (  <--R-- 5'   <--R-- 5'  )

    };


    /**
    * Computes the pair orientation of the given SAMRecord.
    * @param r
    * @return PairOrientation of the given SAMRecord.
    * @throws IllegalArgumentException If the record is not a paired read, or
    * one or both reads are unmapped.
    */
    public static PairOrientation getPairOrientation(SAMRecord r)
    {
        final boolean readIsOnReverseStrand = r.getReadNegativeStrandFlag();

        if(r.getReadUnmappedFlag() || !r.getReadPairedFlag() || r.getMateUnmappedFlag()) {
            throw new IllegalArgumentException("Invalid SAMRecord: " + r.getReadName() + ". This method only works for SAMRecords " +
                    "that are paired reads with both reads aligned.");
        }

        if(readIsOnReverseStrand == r.getMateNegativeStrandFlag() )  {
            return PairOrientation.TANDEM;
        }

        final long positiveStrandFivePrimePos = ( readIsOnReverseStrand
            ?  r.getMateAlignmentStart()  //mate's 5' position  ( x---> )
            :  r.getAlignmentStart() );   //read's 5' position  ( x---> )

        final long negativeStrandFivePrimePos = ( readIsOnReverseStrand
            ?  r.getAlignmentEnd()                                   //read's 5' position  ( <---x )
            :  r.getAlignmentStart() + r.getInferredInsertSize() );  //mate's 5' position  ( <---x )

        return ( positiveStrandFivePrimePos < negativeStrandFivePrimePos
            ? PairOrientation.FR
            : PairOrientation.RF );
    }



    // TODO: KT and TF say this is more complicated than what I have here
    public static boolean isProperPair(final SAMRecord firstEnd, final SAMRecord secondEnd,
                                       final List<PairOrientation> expectedOrientations) {
        // are both records mapped?
        if (firstEnd.getReadUnmappedFlag() || secondEnd.getReadUnmappedFlag()) {
            return false;
        }
        if (firstEnd.getReferenceName().equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
            return false;
        }
        // AND are they both mapped to the same chromosome

        if (!firstEnd.getReferenceName().equals(secondEnd.getReferenceName())) {
            return false;
        }

        // AND is the pair orientation in the set of expected orientations
        final PairOrientation actual = getPairOrientation(firstEnd);
        return expectedOrientations.contains(actual);
    }

    public static void assertMate(final SAMRecord firstOfPair, final SAMRecord secondOfPair) {
        // Validate paired reads arrive as first of pair, then second of pair

        if (firstOfPair == null) {
            throw new PicardException(
                "First record does not exist - cannot perform mate assertion!");
        } else if (secondOfPair == null) {
            throw new PicardException(
                firstOfPair.toString() + " is missing its mate");
        } else if (!firstOfPair.getReadPairedFlag()) {
            throw new PicardException(
                "First record is not marked as paired: " + firstOfPair.toString());
        } else if (!secondOfPair.getReadPairedFlag()) {
            throw new PicardException(
                "Second record is not marked as paired: " + secondOfPair.toString());
        } else if (!firstOfPair.getFirstOfPairFlag()) {
            throw new PicardException(
                "First record is not marked as first of pair: " + firstOfPair.toString());
        } else if (!secondOfPair.getSecondOfPairFlag()) {
            throw new PicardException(
                "Second record is not marked as second of pair: " + secondOfPair.toString());
        } else if (!firstOfPair.getReadName().equals(secondOfPair.getReadName())) {
            throw new PicardException(
                "First [" + firstOfPair.getReadName() + "] and Second [" +
                    secondOfPair.getReadName() + "] readnames do not match!");
        }
    }

    /**
     * Obtain the secondOfPair mate belonging to the firstOfPair SAMRecord
     * (assumed to be in the next element of the specified samRecordIterator)
     * @param samRecordIterator the iterator assumed to contain the secondOfPair SAMRecord in the
     * next element in the iteration
     * @param firstOfPair the firstOfPair SAMRecord
     * @return the secondOfPair SAMRecord
     * @throws PicardException when the secondOfPair mate cannot be obtained due to assertion failures
     */
    public static SAMRecord obtainAssertedMate(final Iterator<SAMRecord> samRecordIterator,
                                         final SAMRecord firstOfPair) {
        if (samRecordIterator.hasNext()) {
            final SAMRecord secondOfPair = samRecordIterator.next();
            assertMate(firstOfPair, secondOfPair);
            return secondOfPair;
        } else {
            throw new PicardException(
                "Second record does not exist: " + firstOfPair.getReadName());
        }
    }

    /**
     * Compute SAMRecord insert size
     * @param firstEnd
     * @param secondEnd
     * @return note that when storing insert size on the secondEnd, the return value must be negated.
     */
    public static int computeInsertSize(final SAMRecord firstEnd, final SAMRecord secondEnd) {
        if (firstEnd.getReadUnmappedFlag() || secondEnd.getReadUnmappedFlag()) {
            return 0;
        }
        if (!firstEnd.getReferenceName().equals(secondEnd.getReferenceName())) {
            return 0;
        }
        int firstEnd5PrimePosition = firstEnd.getReadNegativeStrandFlag()? firstEnd.getAlignmentEnd(): firstEnd.getAlignmentStart();
        int secondEnd5PrimePosition = secondEnd.getReadNegativeStrandFlag()? secondEnd.getAlignmentEnd(): secondEnd.getAlignmentStart();
        return CoordMath.getLength(firstEnd5PrimePosition, secondEnd5PrimePosition);
    }

    /**
     * Write the mate info for two SAMRecords
     */
    public static void setMateInfo(final SAMRecord rec1, final SAMRecord rec2, final SAMFileHeader header) {
        // If neither read is unmapped just set their mate info
        if (!rec1.getReadUnmappedFlag() && !rec2.getReadUnmappedFlag()) {

            rec1.setMateReferenceIndex(rec2.getReferenceIndex());
            rec1.setMateAlignmentStart(rec2.getAlignmentStart());
            rec1.setMateNegativeStrandFlag(rec2.getReadNegativeStrandFlag());
            rec1.setMateUnmappedFlag(false);
            rec1.setAttribute(SAMTag.MQ.name(), rec2.getMappingQuality());

            rec2.setMateReferenceIndex(rec1.getReferenceIndex());
            rec2.setMateAlignmentStart(rec1.getAlignmentStart());
            rec2.setMateNegativeStrandFlag(rec1.getReadNegativeStrandFlag());
            rec2.setMateUnmappedFlag(false);
            rec2.setAttribute(SAMTag.MQ.name(), rec1.getMappingQuality());
        }
        // Else if they're both unmapped set that straight
        else if (rec1.getReadUnmappedFlag() && rec2.getReadUnmappedFlag()) {
            rec1.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec1.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            rec1.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec1.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            rec1.setMateNegativeStrandFlag(rec2.getReadNegativeStrandFlag());
            rec1.setMateUnmappedFlag(true);
            rec1.setAttribute(SAMTag.MQ.name(), null);
            rec1.setInferredInsertSize(0);

            rec2.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec2.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            rec2.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec2.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            rec2.setMateNegativeStrandFlag(rec1.getReadNegativeStrandFlag());
            rec2.setMateUnmappedFlag(true);
            rec2.setAttribute(SAMTag.MQ.name(), null);
            rec2.setInferredInsertSize(0);
        }
        // And if only one is mapped copy it's coordinate information to the mate
        else {
            final SAMRecord mapped   = rec1.getReadUnmappedFlag() ? rec2 : rec1;
            final SAMRecord unmapped = rec1.getReadUnmappedFlag() ? rec1 : rec2;
            unmapped.setReferenceIndex(mapped.getReferenceIndex());
            unmapped.setAlignmentStart(mapped.getAlignmentStart());

            mapped.setMateReferenceIndex(unmapped.getReferenceIndex());
            mapped.setMateAlignmentStart(unmapped.getAlignmentStart());
            mapped.setMateNegativeStrandFlag(unmapped.getReadNegativeStrandFlag());
            mapped.setMateUnmappedFlag(true);
            mapped.setInferredInsertSize(0);

            unmapped.setMateReferenceIndex(mapped.getReferenceIndex());
            unmapped.setMateAlignmentStart(mapped.getAlignmentStart());
            unmapped.setMateNegativeStrandFlag(mapped.getReadNegativeStrandFlag());
            unmapped.setMateUnmappedFlag(false);
            unmapped.setInferredInsertSize(0);
        }

        final int insertSize = SamPairUtil.computeInsertSize(rec1, rec2);
        rec1.setInferredInsertSize(insertSize);
        rec2.setInferredInsertSize(-insertSize);
    }

    public static void setProperPairAndMateInfo(final SAMRecord rec1, final SAMRecord rec2,
                                                final SAMFileHeader header,
                                                final List<PairOrientation> exepectedOrientations) {
        setMateInfo(rec1, rec2, header);
        boolean properPair =  (!rec1.getReadUnmappedFlag() && !rec2.getReadUnmappedFlag())
            ? isProperPair(rec1, rec2, exepectedOrientations)
            : false;
        rec1.setProperPairFlag(properPair);
        rec2.setProperPairFlag(properPair);
    }
}
