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

import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * @author alecw@broadinstitute.org
 */
public class CigarTest {
    private final TextCigarCodec codec = TextCigarCodec.getSingleton();

    @Test
    public void testPositive() {
        Assert.assertNotNull(codec);
        Assert.assertNull(codec.decode("").isValid(null, -1));
        Assert.assertNull(codec.decode("10M5N1I12M").isValid(null, -1));
        Assert.assertNull(codec.decode("10M1I5N1I12M").isValid(null, -1));
        Assert.assertNull(codec.decode("9M1D5N1I12M").isValid(null, -1));
    }
}
