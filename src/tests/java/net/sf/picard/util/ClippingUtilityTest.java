/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

import net.sf.picard.util.IlluminaUtil.IlluminaAdapterPair;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;
import net.sf.samtools.util.StringUtil;

/**
 *
 */
public class ClippingUtilityTest {

   @Test(dataProvider="clipTestData")
   public void testBasicClip(final String testName, final String read, final String clip, final int minMatch, final double errRate, final int expected) {
       final byte[] r = (read == null) ? null : StringUtil.stringToBytes(read);
       final byte[] c = (clip == null) ? null : StringUtil.stringToBytes(clip);

       final int result = ClippingUtility.findIndexOfClipSequence(r, c, minMatch, errRate);
       Assert.assertEquals(result, expected, testName);
   }

   @Test(dataProvider="clipPairedTestData")
   public void testPairedEndClip(final String testName, final String read1, final String read2, final AdapterPair expected) {

       final SAMRecord rec1 = new SAMRecord(new SAMFileHeader());
       rec1.setReadString(read1);
       rec1.setFirstOfPairFlag(true);
       final SAMRecord rec2 = new SAMRecord(new SAMFileHeader());
       rec2.setReadString(read2);
       rec2.setSecondOfPairFlag(true);

       final AdapterPair result = ClippingUtility.adapterTrimIlluminaPairedReads(rec1, rec2,
               new AdapterPair[] { IlluminaAdapterPair.INDEXED, IlluminaAdapterPair.PAIRED_END });
       if (result != null) {
           Assert.assertEquals(result.getName(), expected.getName(), testName);
       }
       else {
           Assert.assertEquals(result, expected, testName);
       }
   }

    @DataProvider(name="clipTestData")
    public Object[][] getClipTestData() {
        final String FORWARD = IlluminaAdapterPair.PAIRED_END.get3PrimeAdapter();
        final String SE_FORWARD = IlluminaAdapterPair.SINGLE_END.get3PrimeAdapter();
        final String REVERSE = IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterInReadOrder();
        return new Object[][] { 
            new Object[] {"Simple test 1", "AAAAACCCCCAGATCGGAAGAGCG", "AGATCGGAAGAGCG", 6, 0.15, 10},
            new Object[] {"Simple test 2", "AAAAACCCCCGGGGGAGATCGGAAGAGCG", "AGATCGGAAGAGCG", 6, 0.15, 15},
            new Object[] {"No adapter", "AAAAACCCCCGGGGGTTTTT", "AGATCGGAAGAGCG", 6, 0.15, -1},
            new Object[] {"Partial adapter", "AAAAACCCCCAGATCGGAA", "AGATCGGAAGAGCG", 6, 0.15, 10},
// no longer support clips in middle of read
//          new Object[] {"Adapter+Primer", "AAAAACCCCCAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG", "AGATCGGAAGAGCG", 6, 0.15, 10},
            new Object[] {"No sequence", null, "AGATCGGAAGAGCG", 6, 0.15, -1},
            new Object[] {"Read too short", "AGATCGGAAG", "AGATCGGAAGAGCG", 11, 0.15, -1},
            new Object[] {"All No Calls", "AAACCCNNNNNNNNNNNNNNNNNNNN", "AGATCGGAAGAGCG", 6, 0.15, -1},
            new Object[] {"From Test Data1", "CGGCATTCCTGCTGAACCGAGATCGGAAGAGCGTCGTGTAGGGAAAGGGGGTGGATCTCGGTGGGCGGCGTGTTGT", REVERSE, 5, 0.15, 19},
            new Object[] {"From Test Data2a", "CGGAAGAGCGGTTCAGCAGGAATGCCGAGATCGGAA", REVERSE, 5, 0.14, 27},  // XT:i:28
            new Object[] {"From Test Data2b", "CGGAAGAGCGGTTCAGCAGGAATGCCGAGATCGGAA", REVERSE, 10, 0.14, -1},  // only matches 9
            new Object[] {"From PE Test Data1", "CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGGT", FORWARD, 5, 0.14, 19},
            new Object[] {"From PE Test Data2", "CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGGT", REVERSE, 5, 0.14, -1},
            new Object[] {"From Test 8-clipped", "TGGGGTGGTTATTGTTGATTTTTGTTTGTGTGTTAGGTTGTTTGTGTTAGTTTTTTATTTTATTTTCGAGATCGGA", FORWARD, 8, 0.14, 68},
            new Object[] {"50% can be bad", "AAAAACCCCCAGATCGGAAGAGCG", "AGATCGGAAGAGCG", 5, 0.5, 18},        // 18?

            new Object[] {"From 30E54AAXX.5.a", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTG", "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", 10, 0.14, 0},  // 0!? From KT test case
            new Object[] {"From 30E54AAXX.5.b", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTG", SE_FORWARD, 10, 0.14, 0},  // 1?? From KT test case
        };
    }

    @DataProvider(name="clipPairedTestData")
    public Object[][] getClipPairedTestData() {
        return new Object[][] {
            new Object[] {"Basic positive paired test matching",    "CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG", "GCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCGG", IlluminaAdapterPair.INDEXED},
            // Tbis matches on one side and matches at higher stringency on that one side
            new Object[] {"Basic positive one-sided test matching", "CTACTGGCGCTGAAAAGATCGGAAGAGCGGTTCAGC", "AAAAAATTTTTTCCCCCCGGGGGGAAAAAATTTTTT", IlluminaAdapterPair.PAIRED_END},
            // Tbis matches on one side and does not match at higher stringency on that one side
            new Object[] {"Basic negative one-sided test matching", "GGCGCTGAAACTACTGGCGCTGAAAAGATCGGAAGA", "AAAAAATTTTTTCCCCCCGGGGGGAAAAAATTTTTT", null},
            // These match but at different positions.  No clip should be done.
            new Object[] {"Mis-match paired test matching",    "CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG", "AGCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCG", null},
        };
    }
}
