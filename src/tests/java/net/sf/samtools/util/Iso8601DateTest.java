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
package net.sf.samtools.util;

import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.Date;

/**
 * @author alecw@broadinstitute.org
 */
public class Iso8601DateTest {
    @Test
    public void testBasic() {
        final String dateStr = "2008-12-15";
        Iso8601Date first = new Iso8601Date(dateStr);
        String firstFormatted = first.toString();
        Iso8601Date second = new Iso8601Date(firstFormatted);
        Assert.assertEquals(first, second);
        String secondFormatted = second.toString();
        Assert.assertEquals(firstFormatted, secondFormatted);
    }

    @Test
    public void testMillisecondTruncation() {
        // Create a Date with milliseconds
        final Date now = new Date();
        if (now.getTime() % 1000 == 0) {
            now.setTime(now.getTime() + 3);
        }
        Iso8601Date isoDate = new Iso8601Date(now);
        Assert.assertEquals(isoDate.getTime() % 1000, 0);
        Assert.assertEquals(isoDate.getTime() / 1000, now.getTime() / 1000);
    }
}
