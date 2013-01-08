package org.broad.tribble.readers;


import org.broad.tribble.TestUtils;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.testng.AssertJUnit.assertTrue;


/**
 * Created by IntelliJ IDEA.
 * User: jrobinso
 * Date: Jul 6, 2010
 * Time: 8:57:40 PM
 * To change this template use File | Settings | File Templates.
 */
public class TabixReaderTest {

    static String tabixFile = TestUtils.DATA_DIR + "tabix/trioDup.vcf.gz";
    static TabixReader tabixReader;
    static List<String> sequenceNames;

    @BeforeClass
    public void setup() throws IOException {
        tabixReader = new TabixReader(tabixFile);
        sequenceNames = new ArrayList<String>(tabixReader.mChr2tid.keySet());
    }

    @AfterClass
    public void teardown() throws Exception {
        // close?  TabixReader doesn't have a close method!
    }

    @Test
    public void testSequenceNames() {


        String[] expectedSeqNames = new String[24];
        for (int i = 1; i < 24; i++) {
            expectedSeqNames[i - 1] = String.valueOf(i);
        }
        expectedSeqNames[22] = "X";
        expectedSeqNames[23] = "Y";
        Assert.assertEquals(expectedSeqNames.length, sequenceNames.size());

        for (String s : expectedSeqNames) {
            Assert.assertTrue(sequenceNames.contains(s));
        }


    }

    /**
     * Test reading a local tabix file
     *
     * @throws java.io.IOException
     */
    @Test
    public void testLocalQuery() throws IOException {

         TabixIteratorLineReader lineReader = new TabixIteratorLineReader(
                tabixReader.query(tabixReader.mChr2tid.get("4"), 320, 330));

        int nRecords = 0;
        String nextLine;
        while ((nextLine = lineReader.readLine()) != null) {
            assertTrue(nextLine.startsWith("4"));
            nRecords++;
        }
        assertTrue(nRecords > 0);


    }

    /**
     * Test reading a tabix file over http
     *
     * @throws java.io.IOException
     */
    @Test
    public void testRemoteQuery() throws IOException {

        String tabixFile = "http://www.broadinstitute.org/igvdata/test/tabix/trioDup.vcf.gz";

        TabixReader tabixReader = new TabixReader(tabixFile);

        TabixIteratorLineReader lineReader = new TabixIteratorLineReader(
                tabixReader.query(tabixReader.mChr2tid.get("4"), 320, 330));

        int nRecords = 0;
        String nextLine;
        while ((nextLine = lineReader.readLine()) != null) {
            assertTrue(nextLine.startsWith("4"));
            nRecords++;
        }
        assertTrue(nRecords > 0);


    }
}
