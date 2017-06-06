/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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
package picard.sam.markduplicates;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Tests a few hand build sam files as they are.
 */
public class AsIsMarkDuplicatesTest {

    static final File TEST_DIR = new File("testdata/picard/sam/MarkDuplicates");

    @DataProvider
    public Object[][] testSameUnclipped5PrimeOppositeStrandData() {
        return new Object[][]{
                new Object[]{new File(TEST_DIR, "sameUnclipped5primeEndv1.sam")},
                new Object[]{new File(TEST_DIR, "sameUnclipped5primeEndv2.sam")},
                new Object[]{new File(TEST_DIR, "sameUnclipped5primeEndCoordinateSortedv1.sam")},
                new Object[]{new File(TEST_DIR, "sameUnclipped5primeEndCoordinateSortedv2.sam")},
                new Object[]{new File(TEST_DIR, "sameUnclipped5primeEndCoordinateSortedv3.sam")},
                new Object[]{new File(TEST_DIR, "sameUnclipped5primeEndCoordinateSortedv4.sam")}
        };
    }

    @Test(dataProvider = "testSameUnclipped5PrimeOppositeStrandData")
    public void testSameUnclipped5PrimeOppositeStrand(final File input) {
        doAsIsTest(input, null);
    }

    @Test
    public void testUnsortedInputWithAssumeSorted() {
        doAsIsTest(new File(TEST_DIR, "querynameGrouped.sam"), SAMFileHeader.SortOrder.queryname);
    }

    private void doAsIsTest(final File input, final SAMFileHeader.SortOrder assumeSortOrder) {

        final AbstractMarkDuplicatesCommandLineProgramTester tester = new BySumOfBaseQAndInOriginalOrderMDTester();
        final SamReader reader = SamReaderFactory.makeDefault().open(input);
        final SAMFileHeader header = reader.getFileHeader().clone();

        if (assumeSortOrder != null) {
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            tester.addArg("ASSUME_SORT_ORDER=" + assumeSortOrder.toString());
        }

        tester.setHeader(header);

        reader.iterator().stream().forEach(tester::addRecord);

        CloserUtil.close(reader);
        tester.setExpectedOpticalDuplicate(0);
        tester.runTest();
    }
}
