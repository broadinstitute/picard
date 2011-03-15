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
 * Comparator for "queryname" ordering of SAMRecords.
 */
public class SAMRecordQueryNameComparator implements SAMRecordComparator {

    public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        int cmp = fileOrderCompare(samRecord1, samRecord2);
        if (cmp != 0) {
            return cmp;
        }

        final boolean r1Paired = samRecord1.getReadPairedFlag();
        final boolean r2Paired = samRecord2.getReadPairedFlag();

        if (r1Paired || r2Paired) {
            if (!r1Paired) return 1;
            else if (!r2Paired) return -1;
            else if (samRecord1.getFirstOfPairFlag()  && samRecord2.getSecondOfPairFlag()) return -1;
            else if (samRecord1.getSecondOfPairFlag() && samRecord2.getFirstOfPairFlag()) return 1;
        }

        if (samRecord1.getReadNegativeStrandFlag() != samRecord2.getReadNegativeStrandFlag()) {
            return (samRecord1.getReadNegativeStrandFlag()? 1: -1);
        }
        if (samRecord1.getNotPrimaryAlignmentFlag() != samRecord2.getNotPrimaryAlignmentFlag()) {
            return samRecord2.getNotPrimaryAlignmentFlag()? -1: 1;
        }
        final Integer hitIndex1 = samRecord1.getIntegerAttribute(SAMTag.HI.name());
        final Integer hitIndex2 = samRecord2.getIntegerAttribute(SAMTag.HI.name());
        if (hitIndex1 != null) {
            if (hitIndex2 == null) return 1;
            else {
                cmp = hitIndex1.compareTo(hitIndex2);
                if (cmp != 0) return cmp;
            }
        } else if (hitIndex2 != null) return -1;
        return 0;
    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.
     *
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return compareReadNames(samRecord1.getReadName(), samRecord2.getReadName());
    }

    /**
     * Encapsulate algorithm for comparing read names in queryname-sorted file, since there have been
     * conversations about changing the behavior.
     */
    public static int compareReadNames(final String readName1, final String readName2) {
        return readName1.compareTo(readName2);
    }
}
