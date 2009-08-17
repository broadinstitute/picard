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
import org.testng.annotations.DataProvider;
import org.testng.Assert;

/**
 * @author alecw@broadinstitute.org
 */
public class StringUtilTest {
    @Test(dataProvider = "provider")
    public void testSplit(final String input, final String[] expectedResult, final boolean concatenateExcess) {
        String[] ret = new String[expectedResult.length];
        int tokensExpected;
        for (tokensExpected = 0; tokensExpected < expectedResult.length && expectedResult[tokensExpected] != null;
             ++tokensExpected) {
        }
        final int tokensFound;
        if (concatenateExcess) {
            tokensFound = StringUtil.splitConcatenateExcessTokens(input, ret, ':');
        } else {
           tokensFound = StringUtil.split(input, ret, ':');
        }
        Assert.assertEquals(tokensFound, tokensExpected);
        Assert.assertEquals(ret, expectedResult);
    }

    @DataProvider(name="provider")
    public Object[][] splitScenarios() {
        return new Object[][] {
                {"A:BB:C", new String[]{"A", "BB", "C"}, false},
                {"A:BB:C", new String[]{"A", "BB", "C"}, true},
                {"A:BB", new String[]{"A", "BB", null}, false},
                {"A:BB", new String[]{"A", "BB", null}, true},
                {"A:BB:", new String[]{"A", "BB", null}, false},
                {"A:BB:", new String[]{"A", "BB", null}, true},
                {"A:BB:C:DDD", new String[]{"A", "BB", "C"}, false},
                {"A:BB:C:DDD", new String[]{"A", "BB", "C:DDD"}, true},
                {"A:", new String[]{"A", null, null}, false},
                {"A:", new String[]{"A", null, null}, true},
                {"A", new String[]{"A", null, null}, false},
                {"A", new String[]{"A", null, null}, true},
                {"A:BB:C", new String[]{"A", "BB", "C"}, false},
                {"A:BB:C:", new String[]{"A", "BB", "C:"}, true}, 
        };
    }
}
