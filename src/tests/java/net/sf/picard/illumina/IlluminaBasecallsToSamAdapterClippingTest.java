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
package net.sf.picard.illumina;

import net.sf.samtools.*;
import net.sf.picard.sam.ReservedTagConstants;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.File;


/**
 * Run IlluminaBasecallsToSam on a sample tests, then sanity-check the generated SAM file
 */
public class IlluminaBasecallsToSamAdapterClippingTest {

    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina");
    //private static final File SEQUENCE_DICTIONARY = new File(TEST_DATA_DIR, "Homo_sapiens_assembly18.seqdict.sam");
    private static final String ALIAS = "myalias";
    private static final String RUN_BARCODE = "305PJAAXX080716";   // todo - not sure if this is right
    private static final String READ_GROUP_NAME = "0";

    /**
     * Run IlluminaBasecallsToSam on a few test cases, and verify that results agree with hand-checked expectation.
     */
    @Test(dataProvider="data")
    public void testBasic(String LANE) throws Exception {
        // Create the SAM file from Gerald output
        final File samFile = File.createTempFile("." + LANE + ".illuminaBasecallsToSam", ".sam");
        samFile.deleteOnExit();
        final String[] illuminaArgv = {
                "BASECALLS_DIR=" + TEST_DATA_DIR + "/IlluminaTests/BasecallsDir",
                "LANE=" + LANE,
                "RUN_BARCODE=" + RUN_BARCODE,
                "OUTPUT=" + samFile,
                "ALIAS=" + ALIAS,
        };
        Assert.assertEquals(new IlluminaBasecallsToSam().instanceMain(illuminaArgv), 0);

        System.out.println ("Ouput Sam file is in " + samFile.getAbsolutePath());

        // Read the file and confirm it contains what is expected
        final SAMFileReader samReader = new SAMFileReader(samFile);

        // look for clipped adaptor attribute in lane 3 PE (2) and in lane 6 (1) non-PE
        int count = 0;   int matchCount = 0;
        for (SAMRecord record : samReader) {
            if (record.getIntegerAttribute(ReservedTagConstants.XT) != null) {
                count ++;
                if ((count == 1 || count ==2) && LANE.equals("3")){
                    Assert.assertEquals (65, (int)record.getIntegerAttribute(ReservedTagConstants.XT));
                    matchCount++;
                } else if (count == 1 && LANE.equals("6")){
                    Assert.assertEquals (141, (int)record.getIntegerAttribute(ReservedTagConstants.XT));
                    matchCount++;
                } else {
                    Assert.assertTrue(false, "Lane:" + LANE + " " + count +
                        " Unexpected adapter trim " + record.getIntegerAttribute(ReservedTagConstants.XT) +
                        " count=" + count + " matchCount=" + matchCount +
                        " at record " + record);
                }
            }
        }
    }

    @DataProvider(name="data")
    private Object[][] getIlluminaBasecallsToSamTestData(){
        return new Object[][] {
                {"1"},
                {"2"},
                {"4"},
                {"5"},
        };
    }


}