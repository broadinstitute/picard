/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import net.sf.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;

/**
 * Very basic test for scatter functionality in IntervalListTools
 */
public class IntervalListToolsTest {

    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/util/");

    @Test
    public void testScatter() throws IOException {
        final File scatterDir = new File(TEST_DATA_DIR, "scatter");

        try {
            final File listToScatter = new File(TEST_DATA_DIR, "scatterable.interval_list");
            final IntervalList gold = IntervalList.fromFile(listToScatter);
            Assert.assertEquals(gold.getUniqueBaseCount(),200, "Wrong unique base count");

            // First test the non-scattering option -- the file that compes out
            // should be the same as the one that went in.
            final File singleFileOutput = File.createTempFile("single", "interval_list");
            singleFileOutput.deleteOnExit();

            int result = new IntervalListTools().instanceMain(new String[] {
                    "INPUT="  + listToScatter.getAbsolutePath(),
                    "OUTPUT=" + singleFileOutput.getAbsolutePath() ,
                    "SCATTER_COUNT=" + 1
            });
            Assert.assertEquals(result, 0);

            final IntervalList single = IntervalList.fromFile(singleFileOutput);
            Assert.assertEquals(gold.size(), single.size());
            Assert.assertEquals(gold.getUniqueBaseCount(), single.getUniqueBaseCount());


            // Next we test scattering.  This should result in the first
            // scattered file having 2 intervals, the second of which is exactly
            // one base long.
            scatterDir.mkdir();
            result = new IntervalListTools().instanceMain(new String[] {
                    "INPUT="  + listToScatter.getAbsolutePath(),
                    "OUTPUT=" + scatterDir.getAbsolutePath() ,
                    "SCATTER_COUNT=" + 2
            });
            Assert.assertEquals(result, 0);
            Assert.assertEquals(scatterDir.listFiles(new FilenameFilter() {
                public boolean accept(File file, String s) {
                    return s.startsWith("temp_") && s.endsWith("_of_2");
                }
            }).length, 2, "Invalid # of scattered subdirectories found");
            IntervalList lists[] = getIntervalLists(scatterDir, 2);
            Assert.assertEquals(lists[0].size(), 2);
            Assert.assertEquals(lists[0].getIntervals().get(1).length(), 1, "Length of split interval list is wrong");
            Assert.assertEquals(lists[0].getUniqueBaseCount(), 100, "Unique base count in first interval is wrong.");
            final IntervalList second = IntervalList.fromFile(IntervalListTools.getScatteredFileName(scatterDir, 2, 1));
            Assert.assertEquals(lists[1].getUniqueBaseCount(), 100, "Unique base count in second interval is wrong.");

            // Test for when scattering breaks exactly at the end of an interval
            result = new IntervalListTools().instanceMain(new String[] {
                    "INPUT="  + listToScatter.getAbsolutePath(),
                    "OUTPUT=" + scatterDir.getAbsolutePath() ,
                    "SCATTER_COUNT=" + 4
            });
            Assert.assertEquals(result, 0);
            Assert.assertEquals(scatterDir.listFiles(new FilenameFilter() {
                public boolean accept(File file, String s) {
                    return s.startsWith("temp_") && s.endsWith("_of_4");
                }
            }).length, 4, "Invalid # of scattered subdirectories found");
            lists = getIntervalLists(scatterDir, 4);
            Assert.assertEquals(lists[0].size(), 1);
            Assert.assertEquals(lists[1].size(), 2);
            Assert.assertEquals(lists[2].size(), 2);
            Assert.assertEquals(lists[3].size(), 1);
            Assert.assertEquals(lists[3].getIntervals().get(0).getStart(), 30200);
            Assert.assertEquals(lists[0].getUniqueBaseCount() + lists[1].getUniqueBaseCount() +
                    lists[2].getUniqueBaseCount() + lists[3].getUniqueBaseCount(), 200);
        }
        finally {
            TestUtil.recursiveDelete(scatterDir);
        }


    }


    // Gets all the interval lists for a given scatter count in a directory
    final IntervalList[] getIntervalLists(final File scatterDir, final int scatterCount) {
        final IntervalList result[] = new IntervalList[scatterCount];
        for (int i = 0; i < scatterCount; i++) {
            result[i] = IntervalList.fromFile(IntervalListTools.getScatteredFileName(scatterDir, scatterCount, i+1));
        }
        return result;
    }
}
