package org.broad.tribble.readers;

import org.broad.tribble.TestUtils;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * User: jacob
 * Date: 2012/05/09
 */
public class PositionalBufferedStreamTest {

    InputStream FileIs;
    long expectedBytes;


    @BeforeMethod
    public void setUp() throws Exception {
        File fi =  new File(TestUtils.DATA_DIR + "test.bed");
        FileIs = new FileInputStream(fi);
        expectedBytes = fi.length();
    }

    @AfterMethod
    public void tearDown() throws Exception {
        if(FileIs != null){
            FileIs.close();
            FileIs = null;
        }
    }

    @Test
    public void testPeek() throws Exception{
        int trials = 10;
        PositionalBufferedStream is = new PositionalBufferedStream(FileIs);
        int bb = is.peek();
        for(int ii=0; ii < trials; ii++){
            Assert.assertEquals(is.peek(), bb);
            Assert.assertEquals(is.getPosition(), 0);
        }

        while((bb = is.peek()) >= 0){
            Assert.assertEquals(is.read(), bb);
        }
    }

    @Test
    public void testIsDone() throws Exception{
        PositionalBufferedStream is = new PositionalBufferedStream(FileIs);
        while(!is.isDone()){
            is.read();
        }
        Assert.assertTrue(is.isDone());
        Assert.assertEquals(is.getPosition(), expectedBytes);
    }

    @Test
    public void testReadCorrectNumberBytes() throws Exception{
        int[] bufSizes= new int[]{5, 20, 60, 120, 131, 150, 200, 1000, 10000, 20000, 512000, 2 << 20};
        for(Integer bufSize: bufSizes){
            setUp();
            tstReadCorrectNumberBytes(bufSize);
            tearDown();
        }
    }

    public void tstReadCorrectNumberBytes(int bufferSize) throws Exception{
        InputStream is = new PositionalBufferedStream(FileIs, bufferSize);
        long count = 0;
        while(is.read() >= 0){
            count++;
        }

        Assert.assertEquals(count, expectedBytes);
    }

    @DataProvider(name = "ReadBytesTestData")
    public Object[][] createReadBytesTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();

        for ( int byteReadSize : Arrays.asList(5, 10, 100, 255) )
            for ( int bufSize : Arrays.asList(1, 10, 100, 1000) )
                tests.add( new Object[]{ (Integer)byteReadSize, (Integer)bufSize });

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "ReadBytesTestData")
    public void testReadBytes(final int byteReadSize, final int bufsize) throws Exception {
        final byte[] bytes = new byte[255];
        for ( int i = 0; i < bytes.length; i++ ) bytes[i] = (byte)i;
        final ByteArrayInputStream bais = new ByteArrayInputStream(bytes);

        final byte[] readBytes = new byte[byteReadSize];
        final PositionalBufferedStream pbs = new PositionalBufferedStream(bais, bufsize);

        int i = 0;
        while ( i < 255 ) {
            final int expectedBytesToRead = Math.min(255 - i, readBytes.length);
            final int nBytesRead = pbs.read(readBytes);
            Assert.assertEquals(nBytesRead, expectedBytesToRead, "Didn't read as many bytes as expected from PBS");

            for ( int j = 0; j < nBytesRead; j++ )
                Assert.assertEquals(readBytes[j], bytes[i+j], "Bytes read not those expected");

            i += nBytesRead;
        }
    }
}
