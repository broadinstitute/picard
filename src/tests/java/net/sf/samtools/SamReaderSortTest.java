package net.sf.samtools;

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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Tests for the implementation of SAMRecordIterator in SAMFileReader
 *
 * @author ktibbett@broadinstitute.org
 */
public class SamReaderSortTest {

    public static final String COORDINATE_SORTED_FILE = "testdata/net/sf/samtools/coordinate_sorted.sam";
    public static final String QUERYNAME_SORTED_FILE = "testdata/net/sf/samtools/queryname_sorted.sam";
    public static final String QUERYNAME_SORTED_NO_HEADER_SORT = "testdata/net/sf/samtools/unsorted.sam";

    @Test(expectedExceptions = IllegalStateException.class)
    public void testSortsDisagree() throws Exception {
        SAMRecordIterator it = new SAMFileReader(new File(COORDINATE_SORTED_FILE)).iterator();
        try {
            it.assertSorted(SAMFileHeader.SortOrder.queryname);
            while(it.hasNext()) {
                it.next();
            }
            Assert.fail("Queryname assertion should have failed on coordinate sorted file but didn't");
        }
        finally {
            it.close();
        }
    }

    @Test(dataProvider="validSorts")
    public void testSortAssertionValid(String file, SAMFileHeader.SortOrder order) {
        SAMRecordIterator it = new SAMFileReader(new File(file)).iterator();
        try {
            it.assertSorted(order);
            while (it.hasNext()) {
                it.next();
            }
        }
        finally {
            it.close();
        }
    }

    @DataProvider(name="validSorts")
    public Object[][] getValidSorts() {
        return new Object[][] {
            {COORDINATE_SORTED_FILE, SAMFileHeader.SortOrder.coordinate},
            {QUERYNAME_SORTED_FILE, SAMFileHeader.SortOrder.queryname},
            {QUERYNAME_SORTED_NO_HEADER_SORT, SAMFileHeader.SortOrder.queryname},
            {COORDINATE_SORTED_FILE, SAMFileHeader.SortOrder.unsorted}
        };
    }


    @Test(dataProvider="invalidSorts", expectedExceptions = IllegalStateException.class)
    public void testSortAssertionFails(String file, SAMFileHeader.SortOrder order) throws Exception {
        SAMRecordIterator it = new SAMFileReader(new File(file)).iterator();
        try {
            it.assertSorted(order);
            while (it.hasNext()) {
                it.next();
            }
            Assert.fail("Iterated successfully over " + file + " with invalid sort assertion: " + order.name());
        }
        finally {
            it.close();
        }
    }

    @DataProvider(name="invalidSorts")
    public Object[][] getInvalidSorts() {
        return new Object[][] {
            {QUERYNAME_SORTED_NO_HEADER_SORT, SAMFileHeader.SortOrder.coordinate}
        };
    }
}
