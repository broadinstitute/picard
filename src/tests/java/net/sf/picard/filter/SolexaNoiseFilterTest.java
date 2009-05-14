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

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;
import net.sf.samtools.SAMRecordSetBuilder;
import net.sf.samtools.SAMRecord;

/**
 * Basic test for the SolexaNoiseFilter
 */
public class SolexaNoiseFilterTest {

    private final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
    private final SolexaNoiseFilter filter = new SolexaNoiseFilter();

    /**
     * Basic positive and negative tests for the PolyANoiseFilter
     *
     * @param sequence          The sequence to be tested
     * @param expectedResult    The expected result (true is the sequence should match the filter, otherwise false)
     */
    @Test(dataProvider="data")
    public void testSolexaNoiseFilter(final String testName, final String sequence, final boolean expectedResult) {
        builder.addUnmappedFragment("testfrag");
        final SAMRecord record = builder.iterator().next();
        record.setReadString(sequence);
        Assert.assertEquals(filter.filterOut(record), expectedResult, testName);
    }


    /**
     * Data for various sequences which may or may not match the filter.
     */
    @DataProvider(name = "data")
    private Object[][] getSolexaNoiseTestData()
    {
        return new Object[][]{
            {"36-base read all a's filter out", "AAAAAaaaaaAAAAAAAAAAAAAAAAAAAAaaaaaa", true},
            {"36-base read with n, filter out", "AAAAAaaaaaAAAAAAAAAAAAAAAAAAAAaaaaan", true}, 
            {"51-base read, final base mismatch", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAT", false},
            {"51-base read, middle base mismatch", "aaaaaaaaaaaaaaaaaaaaaaaaaaTaaaaaaaaaaaaaaaaaaaaaaaa", false},
            {"76-base read, a's and n's, filter out",
                    "aaaaaaaaaaaaaaaaaNNaaaaaaaaaaaaaaaaaaaaaanaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", true},
            {"76-base doesn't match",
                    "NNNATAAAnnnnnnnnnnTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", false},
        };
    }

}
