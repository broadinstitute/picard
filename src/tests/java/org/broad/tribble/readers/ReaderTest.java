package org.broad.tribble.readers;


import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Tests for streams and readers
 */
public class ReaderTest {
    @BeforeClass
    public void setup() throws IOException {
    }

    @AfterClass
    public void teardown() throws Exception {

    }

    @Test
    public void testMultipleLines() throws IOException {
        testStream("line 1\nline2\n");
    }

    @Test
    public void testSingleLine() throws IOException {
        testStream("line 1\n");
    }

    @Test
    public void testEmpty() throws IOException {
        testStream("");
    }


    @Test
    public void testLotsOfLines() throws IOException {
        final StringBuilder b = new StringBuilder();
        for ( int i = 0; i < 10000; i++ ) {
            b.append("line " + i + "\n");
        }
        testStream(b.toString());
    }

    @Test
    public void testMassiveLines() throws IOException {
        final StringBuilder b = new StringBuilder();
        for ( int i = 0; i < 10; i++ ) {
            for ( int j = 0; j < 1000000; j++) {
                b.append(i + "." + j);
            }
            b.append("\n");
        }
        testStream(b.toString());
    }

    @Test
    public void testSkip() throws IOException {
        for ( int skipSizeBase : Arrays.asList(0, 10, 100, 1000, 10000, 1000000)) {
            for ( int skipSizeAdd = 0; skipSizeAdd < 10; skipSizeAdd ++ ) {
                final int skipSize = skipSizeBase + skipSizeAdd;
                final byte[] bytes = new byte[skipSize+2];
                Arrays.fill(bytes, 0, skipSize, (byte)0);
                bytes[skipSize] = 1;
                bytes[skipSize+1] = 2;

                final InputStream is = new ByteArrayInputStream(bytes);
                final PositionalBufferedStream pbs = new PositionalBufferedStream(is);
                pbs.skip(skipSize);

                // first value is 1
                Assert.assertTrue(! pbs.isDone());
                Assert.assertEquals(pbs.getPosition(), skipSize);
                Assert.assertEquals(pbs.peek(), 1);
                Assert.assertEquals(pbs.read(), 1);

                Assert.assertTrue(! pbs.isDone());
                Assert.assertEquals(pbs.getPosition(), skipSize + 1);
                Assert.assertEquals(pbs.peek(), 2);
                Assert.assertEquals(pbs.read(), 2);

                Assert.assertTrue(pbs.isDone());
            }
        }
    }

    private void testStream(final String s) throws IOException {
        testStream(s.getBytes());
        testLineReader(s);
    }

    private void testStream(final byte[] bytes) throws IOException {
        final InputStream is = new ByteArrayInputStream(bytes);
        final PositionalBufferedStream pbs = new PositionalBufferedStream(is);

        int bytePos = 0;
        while ( ! pbs.isDone() ) {
            Assert.assertTrue(bytePos < bytes.length);

            // test position
            Assert.assertEquals(pbs.getPosition(), bytePos);

            // test peek
            final byte atPos = bytes[bytePos];
            Assert.assertEquals(toByte(pbs.peek()), atPos);
            // test position
            Assert.assertEquals(pbs.getPosition(), bytePos);

            // test read
            Assert.assertEquals(toByte(pbs.read()), atPos);
            bytePos++;
            // test position
            Assert.assertEquals(pbs.getPosition(), bytePos);

            // test repeek
            if ( bytePos < bytes.length ) {
                Assert.assertEquals(toByte(pbs.peek()), bytes[bytePos]);
                // test position
                Assert.assertEquals(pbs.getPosition(), bytePos);
            }
        }

        Assert.assertEquals(bytePos, bytes.length);
        pbs.close();
    }

    private void testLineReader(final String lines) throws IOException {
        // read all of the lines into the
        final BufferedReader br = new BufferedReader(new StringReader(lines));
        final List<String> eachLine = new ArrayList<String>();
        while (true) {
            final String line = br.readLine();
            if ( line == null ) break;
            eachLine.add(line);
        }

        final byte[] bytes = lines.getBytes();
        final InputStream is = new ByteArrayInputStream(bytes);
        final PositionalBufferedStream pbs = new PositionalBufferedStream(is);
        final AsciiLineReader alr = new AsciiLineReader(pbs);

        int bytePos = 0, linePos = 0;
        while ( ! pbs.isDone() ) {
            Assert.assertTrue(bytePos < bytes.length);

            // test position
            Assert.assertEquals(pbs.getPosition(), bytePos);

            // test read
            final String readLine = alr.readLine();
            Assert.assertEquals(readLine, eachLine.get(linePos));
            linePos++;

            bytePos += readLine.length() + 1; // 1 for the terminator
            // test position
            Assert.assertEquals(pbs.getPosition(), bytePos);
        }

        Assert.assertEquals(linePos, eachLine.size());
        Assert.assertEquals(bytePos, bytes.length);
        pbs.close();
    }

    private final byte toByte(int i) {
        return (byte)(i & 0xFF);
    }
}
