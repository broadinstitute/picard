package org.broad.tribble.readers;

import org.broad.tribble.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * @author mccowan
 */
public class LongLineBufferedReaderTest {

    /**
     * Test that we read the correct number of lines
     * from a file
     * @throws Exception
     */
    @Test
    public void testReadLines() throws Exception {
        String filePath = TestUtils.DATA_DIR + "large.txt";
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(filePath)));
        LongLineBufferedReader testReader = new LongLineBufferedReader(new InputStreamReader(new FileInputStream(filePath)));
        String line;
        while((line = reader.readLine()) != null){
            Assert.assertEquals(testReader.readLine(), line);
        }
        Assert.assertNull(testReader.readLine());
    }
}
