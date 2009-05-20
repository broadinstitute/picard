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

import org.testng.annotations.BeforeTest;
import org.testng.annotations.AfterTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.util.Arrays;
import java.util.Random;

/**
 * @author alecw@broadinstitute.org
 */
public class SortingLongCollectionTest {
    // Create a separate directory for files so it is possible to confirm that the directory is emptied
    private final File tmpDir = new File(System.getProperty("java.io.tmpdir"), "SortingCollectionTest");

    @BeforeTest
    void setup() {
        // Clear out any existing files if the directory exists
        if (tmpDir.exists()) {
            for (final File f : tmpDir.listFiles()) {
                f.delete();
            }
        }
        tmpDir.mkdir();
    }

    @AfterTest
    void tearDown() {
        if (!tmpDir.exists()) {
            // I don't know why it wouldn't exist, but sometimes it doesn't, and it causes the unit test
            // to fail.  AW 20-May-2009
            return;
        }
        for (final File f : tmpDir.listFiles()) {
            f.delete();
        }
        tmpDir.delete();
    }

    private boolean tmpDirIsEmpty() {
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
     * Generate some values, put into SortingLongCollection, confirm that the right number of
     * values come out, and in the right order.
     * @param numValuesToGenerate
     * @param maxValuesInRam
     */
    @Test(dataProvider = "test1")
    public void testPositive(final String testName, final int numValuesToGenerate, final int maxValuesInRam) {
        final long[] values = new long[numValuesToGenerate];
        int numStringsGenerated = 0;
        final SortingLongCollection sortingCollection = new SortingLongCollection(maxValuesInRam, tmpDir);
        final Random valueGenerator = new Random(123);
        for (int i = 0; i < numValuesToGenerate; ++i) {
            final long value = valueGenerator.nextLong();
            sortingCollection.add(value);
            values[numStringsGenerated++] = value;
        }
        Arrays.sort(values);

        Assert.assertEquals(tmpDirIsEmpty(), numValuesToGenerate <= maxValuesInRam);
        assertIteratorEqualsList(values, sortingCollection);

        sortingCollection.cleanup();
        Assert.assertTrue(tmpDirIsEmpty());
    }

    private void assertIteratorEqualsList(final long[] values, final SortingLongCollection sortingCollection) {
        int i = 0;
        sortingCollection.doneAddingStartIteration();
        while (sortingCollection.hasNext()) {
            Assert.assertEquals(sortingCollection.next(), values[i++], "values failed.  i: " + (i-1) + "; values[i]" + values[i-1]);
        }
        Assert.assertEquals(i, values.length);
    }
}
