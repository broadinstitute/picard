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
package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static net.sf.picard.util.CollectionUtil.makeList;

/**
* @author alecw@broadinstitute.org
*/
public class CycleIlluminaFileMapTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/L001");
    private static final File ZERO_LENGTH_TEST_DATA_DIR = new File("testdata/net/sf/picard/illumina/IlluminaTests/L002");

    private static String laneToDir(int lane) {
        String outStr = String.valueOf(lane);
        while(outStr.length() < 3) {
            outStr = "0" + outStr;
        }
        return "L" + outStr;
    }

    private static String constructPathString(int lane, int cycle) {
        return TEST_DATA_DIR + "/" + laneToDir(lane) + "/C" + cycle + ".1";
    }

    @DataProvider(name = "iteratorTestData")
    public Object [][] iteratorTestData() {
        return new Object[][] {
            new Object[] {
                TEST_DATA_DIR, 1, 1, ".cif",
                    makeList(new File(TEST_DATA_DIR + "/C1.1", "s_1_1.cif"),
                             new File(TEST_DATA_DIR + "/C2.1", "s_1_1.cif"),
                             new File(TEST_DATA_DIR + "/C3.1", "s_1_1.cif"),
                             new File(TEST_DATA_DIR + "/C4.1", "s_1_1.cif"))
            },
            new Object[] {
                TEST_DATA_DIR, 1, 2, ".cnf",
                    makeList(new File(TEST_DATA_DIR + "/C1.1", "s_1_2.cnf"),
                             new File(TEST_DATA_DIR + "/C2.1", "s_1_2.cnf"),
                             new File(TEST_DATA_DIR + "/C3.1", "s_1_2.cnf"),
                             new File(TEST_DATA_DIR + "/C4.1", "s_1_2.cnf"))
            },
            new Object[] {
                ZERO_LENGTH_TEST_DATA_DIR, 1, 1, ".cif", new ArrayList<File>()
            },
            new Object[] {
                TEST_DATA_DIR, 1, 3, ".cnf", new ArrayList<File>()
            },

            new Object[] {
                TEST_DATA_DIR, 2, 1, ".cnf", new ArrayList<File>()
            },
        };
    }

    @Test(dataProvider = "iteratorTestData")
    public void cycledFilesIteratorTest(final File parentDir, final int lane, final int tile, final String fileType, final List<File> expectedFiles) {
        final CycleFilesIterator toTest = new CycleFilesIterator(parentDir, lane, tile, fileType);
        Iterator<File> expectedIter = expectedFiles.iterator();
        while(expectedIter.hasNext()) {
            final File currentExpected = expectedIter.next();
            Assert.assertTrue(toTest.hasNext(), "CycledFileSetIterator is missing file: " + currentExpected.getAbsolutePath());
            Assert.assertEquals(toTest.next().getAbsolutePath(), currentExpected.getAbsolutePath());
        }
        
        Assert.assertTrue(!toTest.hasNext(), "CycledFilesIterator has extra files!");
    }

    @Test
    public void passingAssertCycledIlluminaFileMapTest() {
        final CycleIlluminaFileMap fileMap = new CycleIlluminaFileMap();
        fileMap.put(1, new CycleFilesIterator(TEST_DATA_DIR, 1, 1, ".cif"));
        fileMap.put(2, new CycleFilesIterator(TEST_DATA_DIR, 1, 2, ".cif"));
        fileMap.assertValid(makeList(1,2), 4);
    }

    @Test(expectedExceptions = PicardException.class)
    public void tileFailingAssertCycledIlluminaFileMapTest() {
        final CycleIlluminaFileMap fileMap = new CycleIlluminaFileMap();
        fileMap.put(1, new CycleFilesIterator(TEST_DATA_DIR, 1, 1, ".cif"));
        fileMap.put(2, new CycleFilesIterator(TEST_DATA_DIR, 1, 2, ".cif"));
        fileMap.assertValid(makeList(1,2,3), 4);
    }

    @Test(expectedExceptions = PicardException.class)
    public void cycleFailingAssertCycledIlluminaFileMapTest() {
        final CycleIlluminaFileMap fileMap = new CycleIlluminaFileMap();
        fileMap.put(1, new CycleFilesIterator(TEST_DATA_DIR, 1, 1, ".cif"));
        fileMap.put(2, new CycleFilesIterator(TEST_DATA_DIR, 1, 2, ".cif"));
        fileMap.assertValid(makeList(1,2), 5);
    }

}
