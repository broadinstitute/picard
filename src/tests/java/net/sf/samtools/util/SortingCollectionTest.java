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

import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Random;

public class SortingCollectionTest {
    // Create a separate directory for files so it is possible to confirm that the directory is emptied
    private final File tmpDir = new File(System.getProperty("java.io.tmpdir") + "/" + System.getProperty("user.name"), 
            "SortingCollectionTest");
    @BeforeTest void setup() {
        // Clear out any existing files if the directory exists
        if (tmpDir.exists()) {
            for (final File f : tmpDir.listFiles()) {
                f.delete();
            }
        }
        tmpDir.mkdirs();
    }

    @AfterTest void tearDown() {
        System.err.println("In SortingCollectionTest.tearDown.  tmpDir: " + tmpDir);
        for (final File f : tmpDir.listFiles()) {
            f.delete();
        }
        tmpDir.delete();
    }

    private boolean tmpDirIsEmpty() {
        System.err.println("In SortingCollectionTest.tmpDirIsEmpty.  tmpDir: " + tmpDir);
        return tmpDir.listFiles().length == 0;
    }

    @DataProvider(name = "test1")
    public Object[][] createTestData() {
        return new Object[][] {
                {"empty", 0, 100},
                {"less than threshold", 100, 200},
                {"greater than threshold", 550, 100},
                {"threshold multiple", 600, 100},
                {"threshold multiple plus one", 101, 100},
                {"exactly threshold", 100, 100},
        };
    }

    /**
     * Generate some strings, put into SortingCollection, confirm that the right number of
     * Strings come out, and in the right order.
     * @param numStringsToGenerate
     * @param maxRecordsInRam
     */
    @Test(dataProvider = "test1")
    public void testPositive(final String testName, final int numStringsToGenerate, final int maxRecordsInRam) {
        final String[] strings = new String[numStringsToGenerate];
        int numStringsGenerated = 0;
        final SortingCollection<String> sortingCollection = makeSortingCollection(maxRecordsInRam);
        for (final String s : new RandomStringGenerator(numStringsToGenerate)) {
            sortingCollection.add(s);
            strings[numStringsGenerated++] = s;
        }
        Arrays.sort(strings, new StringComparator());

        Assert.assertEquals(tmpDirIsEmpty(), numStringsToGenerate <= maxRecordsInRam);
        assertIteratorEqualsList(strings, sortingCollection);
        assertIteratorEqualsList(strings, sortingCollection);
        
        sortingCollection.cleanup();
        Assert.assertEquals(tmpDir.list().length, 0);
    }

    private void assertIteratorEqualsList(final String[] strings, final SortingCollection<String> sortingCollection) {
        int i = 0;
        for (final String s : sortingCollection) {
            Assert.assertEquals(s, strings[i++]);
        }
        Assert.assertEquals(i, strings.length);
    }

    private SortingCollection<String> makeSortingCollection(final int maxRecordsInRam) {
        return SortingCollection.newInstance(String.class, new StringCodec(), new StringComparator(),
                maxRecordsInRam, tmpDir);
    }

    /**
     * Generate pseudo-random Strings for testing
     */
    static class RandomStringGenerator implements Iterable<String>, Iterator<String> {
        Random random = new Random(0);
        int numElementsToGenerate;
        int numElementsGenerated = 0;

        /**
         * @param numElementsToGenerate Iteration ends after this many have been generated.
         */
        RandomStringGenerator(final int numElementsToGenerate) {
            this.numElementsToGenerate = numElementsToGenerate;
        }

        public Iterator<String> iterator() {
            return this;
        }

        public boolean hasNext() {
            return numElementsGenerated < numElementsToGenerate;
        }

        public String next() {
            ++numElementsGenerated;
            return Integer.toString(random.nextInt());
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    static class StringComparator implements Comparator<String> {

        public int compare(final String s, final String s1) {
            return s.compareTo(s1);
        }
    }

    static class StringCodec implements SortingCollection.Codec<String> {
        ByteBuffer byteBuffer = ByteBuffer.allocate(4);
        OutputStream os;
        InputStream is;

        public SortingCollection.Codec<String> clone() {
            return new StringCodec();
        }

        /**
         * Where to write encoded output
         *
         * @param os
         */
        public void setOutputStream(final OutputStream os) {
            this.os = os;
        }

        /**
         * Where to read encoded input from
         *
         * @param is
         */
        public void setInputStream(final InputStream is) {
            this.is = is;
        }

        /**
         * Write object to file
         *
         * @param val what to write
         */
        public void encode(final String val) {
            try {
                byteBuffer.clear();
                byteBuffer.putInt(val.length());
                os.write(byteBuffer.array());
                os.write(val.getBytes());
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        /**
         * Read the next record from the input stream and convert into a java object.
         *
         * @return null if no more records.  Should throw exception if EOF is encountered in the middle of
         *         a record.
         */
        public String decode() {
            try {
                byteBuffer.clear();
                int bytesRead = is.read(byteBuffer.array());
                if (bytesRead == -1) {
                    return null;
                }
                if (bytesRead != 4) {
                    throw new RuntimeException("Unexpected EOF in middle of record");
                }
                byteBuffer.limit(4);
                final int length = byteBuffer.getInt();
                final byte[] buf = new byte[length];
                bytesRead = is.read(buf);
                if (bytesRead != length) {
                    throw new RuntimeException("Unexpected EOF in middle of record");
                }
                return new String(buf);
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }
    }
}
