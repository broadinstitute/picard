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

/**
 * Comparator for sorting SAMRecords by coordinate.  Note that the header is required because
 * the order of sequences in the header defines the major sort order.
 *
 * Ideally this method would only return 0 for completely equal SAMRecords, so that sort is
 * completely deterministic.  This implementation does not achieve this completely, but it
 * comes pretty close, while avoiding decoding the variable length fields, except for read name,
 * which is decoded if coordinate and strand are equal.
 *
 * Extreme care must be taken to ensure the following:
 * if A == B, then B == A
 * if A < B, then B > A
 * if A < B && B < C, then A < C
 *
 */
public class SAMRecordCoordinateComparator implements SAMRecordComparator {
    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = fileOrderCompare(samRecord1, samRecord2);
        if (cmp != 0) {
            return cmp;
        }
        // Test of negative strand flag is not really necessary, because it is tested
        // with cmp if getFlags, but it is left here because that is the way it was done
        // in the past.
        if (samRecord1.getReadNegativeStrandFlag() == samRecord2.getReadNegativeStrandFlag()) {
            cmp = samRecord1.getReadName().compareTo(samRecord2.getReadName());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getFlags(), samRecord2.getFlags());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getMappingQuality(), samRecord2.getMappingQuality());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getMateReferenceIndex(), samRecord2.getMateReferenceIndex());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getMateAlignmentStart(), samRecord2.getMateAlignmentStart());
            if (cmp != 0) return cmp;
            cmp = compareInts(samRecord1.getInferredInsertSize(), samRecord2.getInferredInsertSize());
            return cmp;

        }
        else return (samRecord1.getReadNegativeStrandFlag()? 1: -1);
    }

    private int compareInts(int i1, int i2) {
        if (i1 < i2) return -1;
        else if (i1 > i2) return 1;
        else return 0;
    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.  If read is paired and unmapped, use the mate mapping to sort.
     *
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        final int refIndex1 = samRecord1.getReferenceIndex();
        final int refIndex2 = samRecord2.getReferenceIndex();
        if (refIndex1 == -1) {
            return (refIndex2 == -1? 0: 1);
        } else if (refIndex2 == -1) {
            return -1;
        }
        final int cmp = refIndex1 - refIndex2;
        if (cmp != 0) {
            return cmp;
        }
        return samRecord1.getAlignmentStart() - samRecord2.getAlignmentStart();
    }
}
