package net.sf.samtools;

/**
 * The MIT License
 * <p/>
 * Copyright (c) 2014 The Broad Institute
 * <p/>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * <p/>
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * <p/>
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

public class SamHeaderRecordComparatorTest {

	@DataProvider(name="UsualSuspects")
	public Object[][] createData() {
		final SAMReadGroupRecord left = new SAMReadGroupRecord("left");
		left.setPlatformUnit("left.1");
		left.setLibrary("library");

		final SAMReadGroupRecord right = new SAMReadGroupRecord("right");
		right.setPlatformUnit("right.1");
		right.setLibrary("library");
		right.setDescription("description");

		return new Object[][] {{ left, right }};
	}

	@Test(dataProvider="UsualSuspects")
	public void testEqualRecords(final SAMReadGroupRecord left, final SAMReadGroupRecord right) {
		final SAMHeaderRecordComparator<SAMReadGroupRecord> comparator = new SAMHeaderRecordComparator<SAMReadGroupRecord>(SAMReadGroupRecord.PLATFORM_UNIT_TAG);
		Assert.assertEquals(0, comparator.compare(left, left)); // see what I did there?
	}

	@Test(dataProvider="UsualSuspects")
	public void testUnequalRecords(final SAMReadGroupRecord left, final SAMReadGroupRecord right) {
		final SAMHeaderRecordComparator<SAMReadGroupRecord> comparator = new SAMHeaderRecordComparator<SAMReadGroupRecord>(SAMReadGroupRecord.PLATFORM_UNIT_TAG);
		Assert.assertTrue(comparator.compare(left, right) < 0);
		Assert.assertTrue(comparator.compare(right, left) > 0);
	}

	@Test(dataProvider="UsualSuspects")
	public void testNullAttributes(final SAMReadGroupRecord left, final SAMReadGroupRecord right) {
		final SAMHeaderRecordComparator<SAMReadGroupRecord> comparator = new SAMHeaderRecordComparator<SAMReadGroupRecord>(SAMReadGroupRecord.FLOW_ORDER_TAG);
		Assert.assertEquals(0, comparator.compare(left, right)); // neither record has this attribute
	}

	@Test(dataProvider="UsualSuspects")
	public void testOneNullAttribute(final SAMReadGroupRecord left, final SAMReadGroupRecord right) {
		final SAMHeaderRecordComparator<SAMReadGroupRecord> comparator = new SAMHeaderRecordComparator<SAMReadGroupRecord>(SAMReadGroupRecord.DESCRIPTION_TAG);
		Assert.assertTrue(comparator.compare(left, right) < 0);
		Assert.assertTrue(comparator.compare(right, left) > 0);
	}
}
