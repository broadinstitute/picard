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
package net.sf.samtools.seekablestream;

import net.sf.samtools.seekablestream.SeekableFTPStream;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.IOException;
import java.net.URL;

/**
 * @author Jim Robinson
 * @since 10/3/11
 */
public class SeekableFTPStreamTest {


    static String urlString = "ftp://ftp.broadinstitute.org/pub/igv/TEST/test.txt";
    static long fileSize = 27;
    static byte[] expectedBytes = "abcdefghijklmnopqrstuvwxyz\n".getBytes();
    SeekableFTPStream stream;

    @BeforeMethod
    public void setUp() throws IOException {
        stream = new SeekableFTPStream(new URL(urlString));

    }

    @AfterMethod
    public void tearDown() throws IOException {
        stream.close();
    }

    @Test
    public void testLength() throws Exception {
        long length = stream.length();
        Assert.assertEquals(fileSize, length);
    }


    /**
     * Test a buffered read.  The buffer is much large than the file size,  assert that the desired # of bytes are read
     *
     * @throws Exception
     */
    @Test
    public void testBufferedRead() throws Exception {

        byte[] buffer = new byte[64000];
        int nRead = stream.read(buffer);
        Assert.assertEquals(fileSize, nRead);

    }

    /**
     * Test requesting a range that extends beyond the end of the file
     */

    @Test
    public void testRange() throws Exception {
        stream.seek(20);
        byte[] buffer = new byte[64000];
        int nRead = stream.read(buffer);
        Assert.assertEquals(fileSize - 20, nRead);

    }

    /**
     * Test requesting a range that begins beyond the end of the file
     */

    @Test
    public void testBadRange() throws Exception {
        stream.seek(30);
        byte[] buffer = new byte[64000];
        int nRead = stream.read(buffer);
        Assert.assertEquals(-1, nRead);
    }


}


