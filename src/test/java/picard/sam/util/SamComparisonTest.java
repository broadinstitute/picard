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
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.SamComparisonMetric;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class SamComparisonTest {
    private static final File TEST_FILES_DIR = new File("testdata/picard/sam/CompareSAMs");

    @DataProvider(name="samComparisonTestData")
    public Object[][] samComparisonTestData() {
        return new Object[][]{
                {"genomic_sorted.sam", "unsorted.sam", false, false, 0, 0, 0, 0, 0, 0, 0, 0, false},
                {"genomic_sorted.sam", "chr21.sam", false, false, 0, 0, 0, 0, 0, 0, 0, 0, false},
                {"genomic_sorted.sam", "bigger_seq_dict.sam", false, false, 0, 0, 0, 0, 0, 0, 0, 0, false},
                {"bigger_seq_dict.sam", "bigger_seq_dict.sam", false, false, 2, 0, 0, 0, 0, 0, 0, 0, true},
                {"genomic_sorted.sam", "genomic_sorted.sam", false, false, 2, 0, 0, 0, 0, 0, 0, 0, true},
                {"genomic_sorted.sam", "has_non_primary.sam", false, false, 2, 0, 0, 0, 0, 0, 0, 0, true},
                {"genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", false, false, 3, 2, 0, 0, 0, 3, 0, 0, false},
                {"genomic_sorted.sam", "genomic_sorted_sam_v1.6.sam", false, false, 2, 0, 0, 0, 0, 0, 0, 0, false},
                {"group_same_coord.sam", "group_same_coord_diff_order.sam", false, false, 3, 0, 0, 0, 0, 1, 2, 0, false},
                {"genomic_sorted_same_position.sam", "genomic_sorted_same_position.sam", false, false, 2, 0, 0, 0, 0, 0, 0, 0, true},
                {"group_same_coord.sam", "diff_coords.sam", false, false, 0, 5, 0, 0, 0, 0, 0, 0, false},
                {"genomic_sorted.sam", "unmapped_first.sam", false, false, 1, 0, 0, 0, 1, 0, 0, 0, false},
                {"genomic_sorted.sam", "unmapped_second.sam", false, false, 1, 0, 0, 0, 1, 0, 0, 1, false},
                {"unmapped_first.sam", "unmapped_second.sam", false, false, 0, 0, 0, 1, 1, 0, 0, 1, false},
                {"unmapped_first.sam", "unmapped_first.sam", false, false, 1, 0, 1, 0, 0, 0, 0, 0, true},
                {"unsorted.sam", "unsorted.sam", false, false, 2, 0, 0, 0, 0, 0, 0, 0, true},
                {"unsorted.sam", "unsorted2.sam", false, false, 0, 1, 0, 0, 0, 0, 1, 0, false},
                {"dup1.sam", "dup2.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 0, true},
                {"dup1.sam", "dup3.sam", true, false, 13, 1, 0, 0, 0, 0, 0, 4, false},
                {"dup1.sam", "dup4.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 2, false},
                {"dup1.sam", "dup5.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"dup1.sam", "dup2.sam", false, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"dup1_queryname.sam", "dup2_queryname.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 0, true},
                {"dup1_queryname.sam", "dup3_queryname.sam", true, false, 13, 1, 0, 0, 0, 0, 0, 4, false},
                {"dup1_queryname.sam", "dup4_queryname.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 2, false},
                {"dup1_queryname.sam", "dup5_queryname.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"dup1_queryname.sam", "dup2_queryname.sam", false, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"genomic_sorted.sam", "mq0_2.sam", false, true, 1, 1, 0, 0, 0, 0, 0, 0, false},
                {"mq0_1.sam", "mq0_2.sam", false, true, 2, 0, 0, 0, 0, 0, 0, 0, true},
                {"mq0_1.sam", "mq0_2.sam", false, false, 1, 1, 0, 0, 0, 0, 0, 0, false}
        };
    }

    @Test(dataProvider="samComparisonTestData")
    public void testSamComparison(
            final String f1,
            final String f2,
            final boolean lenientDup,
            final boolean mq0Match,
            final int expectedMatch,
            final int expectedDiffer,
            final int expectedUnmappedBoth,
            final int expectedUnmappedLeft,
            final int expectedUnmappedRight,
            final int expectedMissingLeft,
            final int expectedMissingRight,
            final int expectedDupDiffer,
            final boolean areEqual) throws IOException
    {
        // compare forward
        testHelper(
                f1, f2, lenientDup, mq0Match,
                expectedMatch, expectedDiffer, expectedUnmappedBoth,
                expectedUnmappedLeft, expectedUnmappedRight,
                expectedMissingLeft, expectedMissingRight, expectedDupDiffer,
                areEqual
        );

        // compare reverse to validate that comparison commutes (swapping inputs and left and right expected values)
        testHelper(
                f2, f1, lenientDup, mq0Match,
                expectedMatch, expectedDiffer, expectedUnmappedBoth,
                expectedUnmappedRight, expectedUnmappedLeft,    // Swap left and right
                expectedMissingRight, expectedMissingLeft, expectedDupDiffer,      // Swap left and right
                areEqual
        );
    }

    private void testHelper(
            final String f1,
            final String f2,
            final boolean lenientDup,
            final boolean mq0Match,
            final int expectedMatch,
            final int expectedDiffer,
            final int expectedUnmappedBoth,
            final int expectedUnmappedLeft,
            final int expectedUnmappedRight,
            final int expectedMissingLeft,
            final int expectedMissingRight,
            final int expectedDupDiffer,
            final boolean areEqual) throws IOException
    {
        try (final SamReader samReader1 = SamReaderFactory.makeDefault().open(new File(TEST_FILES_DIR, f1));
             final SamReader samReader2 = SamReaderFactory.makeDefault().open(new File(TEST_FILES_DIR, f2)))
        {
            final SamComparison samComparison = new SamComparison(samReader1, samReader2,
                    false, false, false, false, false, lenientDup, mq0Match, null, null);

            Assert.assertEquals(areEqual, samComparison.areEqual());
            Assert.assertEquals(expectedMatch, samComparison.getMappingsMatch());
            Assert.assertEquals(expectedDiffer, samComparison.getMappingsDiffer());
            Assert.assertEquals(expectedUnmappedBoth, samComparison.getUnmappedBoth());
            Assert.assertEquals(expectedUnmappedLeft, samComparison.getUnmappedLeft());
            Assert.assertEquals(expectedUnmappedRight, samComparison.getUnmappedRight());
            Assert.assertEquals(expectedMissingLeft, samComparison.getMissingLeft());
            Assert.assertEquals(expectedMissingRight, samComparison.getMissingRight());
            Assert.assertEquals(expectedDupDiffer, samComparison.getDuplicateMarkingsDiffer());

            final File outputFile = File.createTempFile("samComparison", ".txt");
            outputFile.deleteOnExit();
            samComparison.writeReport(outputFile);
            final MetricsFile<SamComparisonMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
            metricsOutput.read(new FileReader(outputFile));
            final SamComparisonMetric metric = metricsOutput.getMetrics().get(0);

            Assert.assertEquals(areEqual, metric.areEqual);
            Assert.assertEquals(expectedMatch, metric.mappingsMatch);
            Assert.assertEquals(expectedDiffer, metric.mappingsDiffer);
            Assert.assertEquals(expectedUnmappedBoth, metric.unmappedBoth);
            Assert.assertEquals(expectedUnmappedLeft, metric.unmappedLeft);
            Assert.assertEquals(expectedUnmappedRight, metric.unmappedRight);
            Assert.assertEquals(expectedMissingLeft, metric.missingLeft);
            Assert.assertEquals(expectedMissingRight, metric.missingRight);
            Assert.assertEquals(expectedDupDiffer, metric.duplicateMarkingsDiffer);

        }
    }

}
