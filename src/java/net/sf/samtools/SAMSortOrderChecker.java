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
 * Encapsulates simple check for SAMRecord order.
 * @author alecw@broadinstitute.org
 */
public class SAMSortOrderChecker {
    private final SAMFileHeader.SortOrder sortOrder;
    private SAMRecord prev;
    private final SAMRecordComparator comparator;

    public SAMSortOrderChecker(final SAMFileHeader.SortOrder sortOrder) {
        this.sortOrder = sortOrder;
        switch (sortOrder) {
            case coordinate:
                comparator = new SAMRecordCoordinateComparator();
                break;
            case queryname:
                comparator = new SAMRecordQueryNameComparator();
                break;
            case unsorted:
            default:
                comparator = null;
                break;
        }
    }

    /**
     * Check if given SAMRecord violates sort order relative to previous SAMRecord.
     * @return True if sort order is unsorted, if this is the first record, or if previous <= rec.
     */
    public boolean isSorted(final SAMRecord rec) {
        if (comparator == null) {
            return true;
        }
        boolean ret = true;
        if (prev != null) {
            ret = comparator.fileOrderCompare(prev, rec) <= 0;
        }
        prev = rec;
        return ret;
    }

    public SAMRecord getPreviousRecord() {
        return prev;
    }

    /**
     * Return the sort key used for the given sort order.  Useful in error messages.
     */
    public String getSortKey(final SAMRecord rec) {
        switch (sortOrder) {

            case coordinate:
                return rec.getReferenceName() + ":" + rec.getAlignmentStart();
            case queryname:
                return rec.getReadName();
            case unsorted:
            default:
                return null;
        }
    }
}
