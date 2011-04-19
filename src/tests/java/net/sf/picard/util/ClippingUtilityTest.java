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

import net.sf.picard.util.IlluminaUtil.AdapterPair;
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
   public void testBasicClip(String testName, String read, String clip, int minMatch, double errRate, int expected) {
       byte[] r = (read == null) ? null : StringUtil.stringToBytes(read);
       byte[] c = (clip == null) ? null : StringUtil.stringToBytes(clip);

       int result = ClippingUtility.findIndexOfClipSequence(r, c, minMatch, errRate);
       Assert.assertEquals(result, expected, testName);
   }

   @Test(dataProvider="clipPairedTestData")
   public void testPairedEndClip(String testName, String read1, String read2, String expected) {

       SAMRecord rec1 = new SAMRecord(new SAMFileHeader());
       rec1.setReadString(read1);
       rec1.setFirstOfPairFlag(true);
       SAMRecord rec2 = new SAMRecord(new SAMFileHeader());
       rec2.setReadString(read2);
       rec2.setSecondOfPairFlag(true);

       String result = ClippingUtility.adapterTrimIlluminaPairedReads(rec1, rec2, AdapterPair.PAIRED_END);
       Assert.assertEquals(result, expected, testName);
   }

    @DataProvider(name="clipTestData")
    public Object[][] getClipTestData() {
        final String FORWARD = AdapterPair.PAIRED_END.get3PrimeAdapter();
        final String SE_FORWARD = AdapterPair.SINGLE_END.get3PrimeAdapter();
        final String REVERSE = AdapterPair.PAIRED_END.get5PrimeAdapterInReadOrder();
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

            new Object[] {"From 30E54AAXX.5 should be 1", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTG", "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", 10, 0.14, 0},  // 0!? From KT test case
            new Object[] {"From 30E54AAXX.5 should be 1", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTGA", "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG", 10, 0.14, -1}, // adapter not at end
            new Object[] {"2From 30E54AAXX.5 should be 1", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTG", SE_FORWARD, 10, 0.14, 0},  // 1?? From KT test case
           new Object[] {"From 30E54AAXX.5", "ATATCTGAAGATCTCGTATGCCGTCTTCTGCTTGAAACAAAAAAATTGAGTCGGTTCATNTTTTCTTTTCTTCCAT", SE_FORWARD, 10, 0.14, -1},  // 1?? From KT test case
        };
    }

    @DataProvider(name="clipPairedTestData")
    public Object[][] getClipPairedTestData() {
        return new Object[][] {
                // todo - test a one-sided match.  test a mismatched-pair.
                // todo distinguish no-match return and matched retur
         //   new Object[] {"Paired test 1", "CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG", "GCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCGGA",
         //           "Adapters mismatch at position 29 CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG and reverse 30 GCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCGGA"},
         new Object[] {"Paired test match", "CTACTGGCGCTGAAACTGAGCAGCCAAGCAGATCGG", "GCTTGGCTGCTCAGTTTCAGCGCCAGTAGAGATCGG", null},
        // new Object[] {"Paired test - one-sided", "C.CCG.......G.GGTGCATGGGCTCCAACGTGGTGTCCTGTGGAGCTGTTGGGCCTGGGCAGGCGGCACAGATC", "TGCCAGTAGTTTTGGGTCAAGCCCTCACCTGATTCCACGCTTCATAGCTTCAGCCGTTCCCATCATACTACTAGCT",
        //         "No adapter match will not trim in paired read of length 76 and length 76  reverse null 76b aligned read. CNCCGNNNNNNNGNGGTGCATGGGCTCCAACGTGGTGTCCTGTGGAGCTGTTGGGCCTGGGCAGGCGGCACAGATC after strict check using minMaxBases=10"},
        };
    }
}
