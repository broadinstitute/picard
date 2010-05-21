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

import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;

public class SeekableBufferedStreamTest {

    private final File BAM_INDEX_FILE = new File("testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam.bai");
    private final File BAM_FILE = new File("testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam");
    private final String BAM_URL_STRING = "http://picard.sourceforge.net/testdata/index_test.bam";

    /**
     * Test reading across a buffer boundary (buffer size is 512000).   The test first reads a range of
     * bytes using an unbuffered stream file stream,  then compares this to results from a buffered http stream.
     *
     * @throws IOException
     */
    @Test
    public void testRandomRead() throws IOException {

        int startPosition = 500000;
        int length = 50000;

        byte[] buffer1 = new byte[length];
        SeekableStream unBufferedStream = new SeekableFileStream(BAM_FILE);
        unBufferedStream.seek(startPosition);
        int bytesRead = unBufferedStream.read(buffer1, 0, length);
        assertEquals(length, bytesRead);

        byte[] buffer2 = new byte[length];
        SeekableStream bufferedStream = new SeekableBufferedStream(new SeekableHTTPStream(new URL(BAM_URL_STRING)));
        bufferedStream.seek(startPosition);
        bytesRead = bufferedStream.read(buffer2, 0, length);
        assertEquals(length, bytesRead);

        assertEquals(buffer1, buffer2);
    }

    /**
     * Test an attempt to read past the end of the file.  The test file is 594,149 bytes in length.  The test
     * attempts to read a 1000 byte block starting at position 594000.  A correct result would return 149 bytes.
     *
     * @throws IOException
     */
    @Test
    public void testEOF() throws IOException {

        int remainder = 149;
        long fileLength = BAM_FILE.length();
        long startPosition = fileLength - remainder;
        int length = 1000;


        byte[] buffer = new byte[length];
        SeekableStream bufferedStream = new SeekableBufferedStream(new SeekableHTTPStream(new URL(BAM_URL_STRING)));
        bufferedStream.seek(startPosition);
        int bytesRead = bufferedStream.read(buffer, 0, length);
        assertEquals(remainder, bytesRead);

        // Subsequent reads should return -1
        bytesRead = bufferedStream.read(buffer, 0, length);
        assertEquals(-1, bytesRead);

    }
}
