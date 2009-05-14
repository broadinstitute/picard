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
package net.sf.picard.filter;

import net.sf.picard.sam.ReservedTagConstants;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

import java.util.List;
import java.util.Arrays;

import net.sf.samtools.SAMRecordSetBuilder;
import net.sf.samtools.SAMRecord;

/**
 * Tests for the TagFilter class
 */
public class TagFilterTest {
    private final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();


    /**
     * Basic positive and negative tests for the TagFilter
     *
     * @param tag               The tag to be tested
     * @param validValues       The values the filter should test for
     * @param testValue         The value to test for in the record
     * @param expectedResult    The expected result (true is the sequence should match the filter, otherwise false)
     */
    @Test(dataProvider="data")
    public void testTagFilter(final String testName, final String tag, final List<Object> validValues,
                              final Object testValue, final boolean expectedResult) {
        final TagFilter filter = new TagFilter(tag, validValues);
        builder.addUnmappedFragment("testfrag");
        final SAMRecord record = builder.iterator().next();
        if (testValue != null) {
            record.setAttribute(tag, testValue);
        }
        Assert.assertEquals(filter.filterOut(record), expectedResult, testName);
    }


    /**
     * Data for various sequences which may or may not match the filter.
     */
    @DataProvider(name = "data")
    private Object[][] getTagFilterTestData()
    {
        return new Object[][]{
            {"Basic positive test", ReservedTagConstants.XN, Arrays.asList(1), 1, true},
            {"Multi-value positive test", ReservedTagConstants.XN, Arrays.asList(1,2,3), 1, true},
            {"Incorrect value negative test", ReservedTagConstants.XN, Arrays.asList(1), 2, false},
            {"Null value negative test", ReservedTagConstants.XN, Arrays.asList(1), null, false} 
        };
    }
}