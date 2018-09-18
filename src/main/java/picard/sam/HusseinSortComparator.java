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
package picard.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SamReader;

import java.io.Serializable;

import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX;
import static htsjdk.samtools.SAMRecord.NO_ALIGNMENT_START;

/**
 * Comparator for sorting SAMRecords by chr + location + bases for unmappeds.
 *
 */
public class HusseinSortComparator implements SAMRecordComparator, Serializable {
    protected SamReader queryReader;

    public HusseinSortComparator(SamReader queryReader) {
        this.queryReader = queryReader;
    }

    public String getMate1Bases(SAMRecord rec) {
        if(rec.getFirstOfPairFlag()) {
            return rec.getReadString();
        } else {
            synchronized(this) {
                SAMRecord mate = this.queryReader.queryMate(rec);
                return mate.getReadString();
            }
        }

    }

    @Override
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {

        if (null == samRecord1.getHeader() || null == samRecord2.getHeader()) {
            throw new IllegalArgumentException("Records must have non-null SAMFileHeaders to be compared");
        }

        int refIndex1 = samRecord1.getReferenceIndex() == NO_ALIGNMENT_REFERENCE_INDEX ? samRecord1.getMateReferenceIndex() : samRecord1.getReferenceIndex();
        int refIndex2 = samRecord2.getReferenceIndex() == NO_ALIGNMENT_REFERENCE_INDEX ? samRecord2.getMateReferenceIndex() : samRecord2.getReferenceIndex();

        int pos1 = samRecord1.getAlignmentStart() == NO_ALIGNMENT_START ? samRecord1.getMateAlignmentStart() : samRecord1.getAlignmentStart();
        int pos2 = samRecord2.getAlignmentStart() == NO_ALIGNMENT_START ? samRecord2.getMateAlignmentStart() : samRecord2.getAlignmentStart();

        if(refIndex1 != refIndex2) {
            return refIndex1 - refIndex2;
        } else if( pos1 != pos2 ) {
            return pos1 - pos2;
        }

        //we're now either dealing with two unmapped reads, or a pair where one is mapped and one isn't
        String bases1 = getMate1Bases(samRecord1);
        String bases2 = getMate1Bases(samRecord2);

        if(!bases1.equals(bases2)) {
            return bases1.compareTo(bases2);
        } else if(samRecord1.getFirstOfPairFlag() != samRecord2.getFirstOfPairFlag()) {
            //by this point we're either an unmapped pair, or a pair mapped to the same place.
            //we'll sort by first-ness, then by secondary-ness
            return samRecord1.getFirstOfPairFlag() ? -1 : 1;
        } else if(samRecord1.isSecondaryOrSupplementary() != samRecord2.isSecondaryOrSupplementary()){
            return samRecord1.isSecondaryOrSupplementary() ? 1 : -1;
        } else {
            return samRecord1.getFlags() - samRecord2.getFlags();
        }
    }

    @Override
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return compare(samRecord1, samRecord2);
    }
}
