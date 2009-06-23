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
        final int cmp = fileOrderCompare(samRecord1, samRecord2);
        if (cmp != 0) {
            return cmp;
        }
        if (samRecord1.getReadPairedFlag()) {
            if (samRecord1.getFirstOfPairFlag() && samRecord2.getSecondOfPairFlag()) {
                return -1;
            }
            else if (samRecord2.getFirstOfPairFlag() && samRecord1.getSecondOfPairFlag()) {
                return 1;
            }
        }

        if (samRecord1.getReadNegativeStrandFlag() == samRecord2.getReadNegativeStrandFlag()) {
            return 0;
        }
        return (samRecord1.getReadNegativeStrandFlag()? 1: -1);
    }

    /**
     * Less stringent compare method than the regular compare.  If the two records
     * are equal enough that their ordering in a sorted SAM file would be arbitrary,
     * this method returns 0.
     *
     * @return negative if samRecord1 < samRecord2,  0 if equal, else positive
     */
    public int fileOrderCompare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
        return samRecord1.getReadName().compareTo(samRecord2.getReadName());
    }
}
