/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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
package picard.sam.util;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

public class SamComparisonTest {
    private static final File TEST_FILES_DIR = new File("testdata/picard/sam/CompareSAMs");

    @DataProvider(name="samComparisonTestData")
    public Object[][] samComparisonTestData() {
        return new Object[][]{
                { "genomic_sorted.sam", "unsorted.sam", 0, 0, 0, 0, 0, 0, 0, false },
                { "genomic_sorted.sam", "chr21.sam", 0, 0, 0, 0, 0, 0, 0, false },
                { "genomic_sorted.sam", "bigger_seq_dict.sam", 0, 0, 0, 0, 0, 0, 0, false },
                { "bigger_seq_dict.sam", "bigger_seq_dict.sam", 2, 0, 0, 0, 0, 0, 0, true },
                { "genomic_sorted.sam", "genomic_sorted.sam", 2, 0, 0, 0, 0, 0, 0, true },
                { "genomic_sorted.sam", "has_non_primary.sam", 2, 0, 0, 0, 0, 0, 0, true },
                { "genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", 3, 2, 0, 0, 0, 3, 0, false },
                { "genomic_sorted.sam", "genomic_sorted_sam_v1.6.sam", 2, 0, 0, 0, 0, 0, 0, false },
                { "group_same_coord.sam", "group_same_coord_diff_order.sam", 3, 0, 0, 0, 0, 1, 2, false },
                { "genomic_sorted_same_position.sam", "genomic_sorted_same_position.sam", 2, 0, 0, 0, 0, 0, 0, true },
                { "group_same_coord.sam", "diff_coords.sam", 0, 5, 0, 0, 0, 0, 0, false },
                { "genomic_sorted.sam", "unmapped_first.sam", 1, 0, 0, 0, 1, 0, 0, false },
                { "genomic_sorted.sam", "unmapped_second.sam", 1, 0, 0, 0, 1, 0, 0, false },
                { "unmapped_first.sam", "unmapped_second.sam", 0, 0, 0, 1, 1, 0, 0, false },
                { "unmapped_first.sam", "unmapped_first.sam", 1, 0, 1, 0, 0, 0, 0, true },
                { "unsorted.sam", "unsorted.sam", 2, 0, 0, 0, 0, 0, 0, true },
                { "unsorted.sam", "unsorted2.sam", 0, 1, 0, 0, 0, 0, 1, false }
        };
    }

    @Test(dataProvider="samComparisonTestData")
    public void testSamComparison(
            final String f1,
            final String f2,
            final int expectedMatch,
            final int expectedDiffer,
            final int expectedUnmappedBoth,
            final int expectedUnmappedLeft,
            final int expectedUnmappedRight,
            final int expectedMissingLeft,
            final int expectedMissingRight,
            final boolean areEqual) throws IOException
    {
        // compare forward
        testHelper(
                f1, f2,
                expectedMatch, expectedDiffer, expectedUnmappedBoth,
                expectedUnmappedLeft, expectedUnmappedRight,
                expectedMissingLeft, expectedMissingRight,
                areEqual
        );

        // compare reverse to validate that comparison commutes (swapping inputs and left and right expected values)
        testHelper(
                f2, f1,
                expectedMatch, expectedDiffer, expectedUnmappedBoth,
                expectedUnmappedRight, expectedUnmappedLeft,    // Swap left and right
                expectedMissingRight, expectedMissingLeft,      // Swap left and right
                areEqual
        );
    }

    private void testHelper(
            final String f1,
            final String f2,
            final int expectedMatch,
            final int expectedDiffer,
            final int expectedUnmappedBoth,
            final int expectedUnmappedLeft,
            final int expectedUnmappedRight,
            final int expectedMissingLeft,
            final int expectedMissingRight,
            final boolean areEqual) throws IOException
    {
        try (final SamReader samReader1 = SamReaderFactory.makeDefault().open(new File(TEST_FILES_DIR, f1));
             final SamReader samReader2 = SamReaderFactory.makeDefault().open(new File(TEST_FILES_DIR, f2)))
        {
            final SamComparison samComparison = new SamComparison(samReader1, samReader2);

            Assert.assertEquals(areEqual, samComparison.areEqual());
            Assert.assertEquals(expectedMatch, samComparison.getMappingsMatch());
            Assert.assertEquals(expectedDiffer, samComparison.getMappingsDiffer());
            Assert.assertEquals(expectedUnmappedBoth, samComparison.getUnmappedBoth());
            Assert.assertEquals(expectedUnmappedLeft, samComparison.getUnmappedLeft());
            Assert.assertEquals(expectedUnmappedRight, samComparison.getUnmappedRight());
            Assert.assertEquals(expectedMissingLeft, samComparison.getMissingLeft());
            Assert.assertEquals(expectedMissingRight, samComparison.getMissingRight());
        }
    }

}
