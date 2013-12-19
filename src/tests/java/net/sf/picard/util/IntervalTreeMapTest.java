/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package net.sf.picard.util;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Iterator;
import java.util.Set;

public class IntervalTreeMapTest {
    @Test
    public void testBasic() {
        IntervalTreeMap<Interval> m=new IntervalTreeMap<Interval>();

        Interval chr1Interval = new Interval("chr1", 1,100);
        m.put(chr1Interval, chr1Interval);
        Interval chr2Interval = new Interval("chr2", 1,200);
        m.put(chr2Interval, chr2Interval);

        final Iterator<Interval> iterator = m.keySet().iterator();
        Assert.assertEquals(iterator.next(), chr1Interval);
        Assert.assertEquals(iterator.next(), chr2Interval);
        Assert.assertFalse(iterator.hasNext());
    }
}
