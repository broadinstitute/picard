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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Tests a few hand build sam files as they are.
 */
public class AsIsMarkDuplicatesTester {

    @DataProvider
    public Object[][] testSameUnclipped5PrimeOppositeStrandData() {
        final File TEST_DIR = new File("testdata/picard/sam/MarkDuplicates");
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

        final AbstractMarkDuplicatesCommandLineProgramTester tester = new BySumOfBaseQAndInOriginalOrderMDTester();

        final SamReader reader = SamReaderFactory.makeDefault().open(input);

        tester.setHeader(reader.getFileHeader());
        reader.iterator().stream().forEach(tester::addRecord);

        CloserUtil.close(reader);
        tester.setExpectedOpticalDuplicate(0);
        tester.runTest();
    }
}


