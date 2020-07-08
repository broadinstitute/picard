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
                {"duplicate_base.sam", "duplicate_four_mismatch_strict.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 0, true},
                {"duplicate_base.sam", "duplicate_four_mismatch_lenient_one_align_differ.sam", true, false, 13, 1, 0, 0, 0, 0, 0, 4, false},
                {"duplicate_base.sam", "duplicate_two_mismatch_lenient.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 2, false},
                {"duplicate_base.sam", "duplicate_four_mismatch_lenient.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"duplicate_base.sam", "duplicate_four_mismatch_strict.sam", false, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_strict_queryname.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 0, true},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_lenient_one_align_differ_queryname.sam", true, false, 13, 1, 0, 0, 0, 0, 0, 4, false},
                {"duplicate_base_queryname.sam", "duplicate_two_mismatch_lenient_queryname.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 2, false},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_lenient_queryname.sam", true, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
                {"duplicate_base_queryname.sam", "duplicate_four_mismatch_strict_queryname.sam", false, false, 14, 0, 0, 0, 0, 0, 0, 4, false},
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
            final boolean lenientLowMQAlignment,
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
                f1, f2, lenientDup, lenientLowMQAlignment,
                expectedMatch, expectedDiffer, expectedUnmappedBoth,
                expectedUnmappedLeft, expectedUnmappedRight,
                expectedMissingLeft, expectedMissingRight, expectedDupDiffer,
                areEqual
        );

        // compare reverse to validate that comparison commutes (swapping inputs and left and right expected values)
        testHelper(
                f2, f1, lenientDup, lenientLowMQAlignment,
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
            final boolean lenientLowMQAlignment,
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
        try (final SamReader samReader1 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(new File(TEST_FILES_DIR, f1));
             final SamReader samReader2 = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(new File(TEST_FILES_DIR, f2)))
        {
            final SAMComparisonArgumentCollection argumentCollection = new SAMComparisonArgumentCollection();
            argumentCollection.LENIENT_DUP = lenientDup;
            argumentCollection.LENIENT_LOW_MQ_ALIGNMENT = lenientLowMQAlignment;
            final SamComparison samComparison = new SamComparison(samReader1, samReader2,
                    null, null, argumentCollection);

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

            Assert.assertEquals(areEqual, metric.ARE_EQUAL);
            Assert.assertEquals(expectedMatch, metric.MAPPINGS_MATCH);
            Assert.assertEquals(expectedDiffer, metric.MAPPINGS_DIFFER);
            Assert.assertEquals(expectedUnmappedBoth, metric.UNMAPPED_BOTH);
            Assert.assertEquals(expectedUnmappedLeft, metric.UNMAPPED_LEFT);
            Assert.assertEquals(expectedUnmappedRight, metric.UNMAPPED_RIGHT);
            Assert.assertEquals(expectedMissingLeft, metric.MISSING_LEFT);
            Assert.assertEquals(expectedMissingRight, metric.MISSING_RIGHT);
            Assert.assertEquals(expectedDupDiffer, metric.DUPLICATE_MARKINGS_DIFFER);

        }
    }

}
