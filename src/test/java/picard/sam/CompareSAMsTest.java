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
package picard.sam;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;

public class CompareSAMsTest extends CommandLineProgramTest {
    private static final File TEST_FILES_DIR = new File("testdata/picard/sam/CompareSAMs");

    public String getCommandLineProgramName() {
        return CompareSAMs.class.getSimpleName();
    }

    private void testHelper(final String f1, final String f2, final int expectedMatch, final int expectedDiffer,
                            final int expectedUnmappedBoth,
                            final int expectedUnmappedLeft, final int expectedUnmappedRight, final int expectedMissingLeft,
                            final int expectedMissingRight, final int expectedDifferQ0, final int expectedDuplicateMarkingsDiffer,
                            final int expectedDuplicateMarkingsDifferSophisticated, final int expectedMappingQualsDiffer, final boolean areEqual) {
        final String[] samFiles = {
                new File(TEST_FILES_DIR, f1).getAbsolutePath(),
                new File(TEST_FILES_DIR, f2).getAbsolutePath()
        };

        // TODO - Should switch over to using invocation via new PicardCommandLine() - BUT the test here is accessing class members directly.
        CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.instanceMain(samFiles);
        Assert.assertEquals(compareSAMs.areEqual(), areEqual);
        Assert.assertEquals(compareSAMs.getMappingsMatch(), expectedMatch);
        Assert.assertEquals(compareSAMs.getMappingsDiffer(), expectedDiffer);
        Assert.assertEquals(compareSAMs.getUnmappedBoth(), expectedUnmappedBoth);
        Assert.assertEquals(compareSAMs.getUnmappedLeft(), expectedUnmappedLeft);
        Assert.assertEquals(compareSAMs.getUnmappedRight(), expectedUnmappedRight);
        Assert.assertEquals(compareSAMs.getMissingLeft(), expectedMissingLeft);
        Assert.assertEquals(compareSAMs.getMissingRight(), expectedMissingRight);
        Assert.assertEquals(compareSAMs.getDifferQ0(),expectedDifferQ0);
        Assert.assertEquals(compareSAMs.getDuplicateMarkingsDiffer(),expectedDuplicateMarkingsDiffer);
        Assert.assertEquals(compareSAMs.getDuplicateMarkingsDifferSophisticated(),expectedDuplicateMarkingsDifferSophisticated);
        Assert.assertEquals(compareSAMs.getMappingQualsDiffer(),expectedMappingQualsDiffer);

        final String[] samFilesReversed = {
                new File(TEST_FILES_DIR, f2).getAbsolutePath(),
                new File(TEST_FILES_DIR, f1).getAbsolutePath()
        };
        compareSAMs = new CompareSAMs();
        compareSAMs.instanceMain(samFilesReversed);
        Assert.assertEquals(compareSAMs.areEqual(), areEqual);
        Assert.assertEquals(compareSAMs.getMappingsMatch(), expectedMatch);
        Assert.assertEquals(compareSAMs.getMappingsDiffer(), expectedDiffer);
        Assert.assertEquals(compareSAMs.getUnmappedBoth(), expectedUnmappedBoth);
        Assert.assertEquals(compareSAMs.getUnmappedLeft(), expectedUnmappedRight);
        Assert.assertEquals(compareSAMs.getUnmappedRight(), expectedUnmappedLeft);
        Assert.assertEquals(compareSAMs.getMissingLeft(), expectedMissingRight);
        Assert.assertEquals(compareSAMs.getMissingRight(), expectedMissingLeft);
        Assert.assertEquals(compareSAMs.getDifferQ0(),expectedDifferQ0);
        Assert.assertEquals(compareSAMs.getDuplicateMarkingsDiffer(),expectedDuplicateMarkingsDiffer);
        Assert.assertEquals(compareSAMs.getDuplicateMarkingsDifferSophisticated(),expectedDuplicateMarkingsDifferSophisticated);
        Assert.assertEquals(compareSAMs.getMappingQualsDiffer(),expectedMappingQualsDiffer);
    }

    @Test
    public void testSortsDifferent() {
        testHelper("genomic_sorted.sam", "unsorted.sam", 0, 0, 0, 0,
                0, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testSequenceDictionariesDifferent1() {
        testHelper("genomic_sorted.sam", "chr21.sam", 0, 0, 0, 0,
                0, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testSequenceDictionariesDifferent2() {
        testHelper("genomic_sorted.sam", "bigger_seq_dict.sam", 0, 0, 0, 0,
                0, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testBiggerSequenceDictionaries() {
        testHelper("bigger_seq_dict.sam", "bigger_seq_dict.sam", 2, 0, 0, 0,
                0, 0, 0, 0,0,0,0,true);
    }

    @Test
    public void testIdentical() {
        testHelper("genomic_sorted.sam", "genomic_sorted.sam", 2, 0, 0, 0,
                0, 0, 0, 0,0,0,0,true);
    }

    @Test
    public void testHasNonPrimary() {
        testHelper("genomic_sorted.sam", "has_non_primary.sam", 2, 0, 0, 0,
                0, 0, 0, 0,0,0,0,true);
    }

    @Test
    public void testMoreOnOneSide() {
        testHelper("genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", 3, 2, 0, 0,
                0, 3, 0, 0,0,0,0,false);
    }

    @Test
    public void testGroupWithSameCoordinate() {
        testHelper("group_same_coord.sam", "group_same_coord_diff_order.sam", 3, 0, 0, 0,
                0, 1, 2, 0,0,0,0,false);
    }

    @Test
    public void testGroupWithSameCoordinateSamePosition() {
        testHelper("genomic_sorted_same_position.sam", "genomic_sorted_same_position.sam", 2, 0, 0, 0,
                0, 0, 0, 0,0,0,0,true);
    }
    @Test
    public void testGroupWithSameCoordinateAndNoMatchInOther() {
        testHelper("group_same_coord.sam", "diff_coords.sam", 0, 5, 0, 0,
                0, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testUnmapped1() {
        testHelper("genomic_sorted.sam", "unmapped_first.sam", 1, 0, 0, 0,
                1, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testUnmapped2() {
        testHelper("genomic_sorted.sam", "unmapped_second.sam", 1, 0, 0, 0,
                1, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testUnmapped3() {
        testHelper("unmapped_first.sam", "unmapped_second.sam", 0, 0, 0, 1,
                1, 0, 0, 0,0,0,0,false);
    }

    @Test
    public void testUnmapped4() {
        testHelper("unmapped_first.sam", "unmapped_first.sam", 1, 0, 1, 0,
                0, 0, 0, 0,0,0,0,true);
    }

    @Test
    public void testUnsorted1() {
        testHelper("unsorted.sam", "unsorted.sam", 2, 0, 0, 0,
                0, 0, 0, 0,0,0,0,true);
    }

    @Test
    public void testUnsorted2() {
        testHelper("unsorted.sam", "unsorted2.sam", 0, 1, 0, 0,
                0, 0, 1, 0,0,0,0,false);
    }

    @Test
    public void testDuplicateMarking1() {
        testHelper("duplicate_marking_1.sam","duplicate_marking_2.sam",4,0,0,0,
                0,0,0,0,2,2,0,false);
    }

    @Test
    public void testDuplicateMarkingCoordinateOrder() {
        testHelper("duplicate_marking_1_CO_sort.sam","duplicate_marking_2_CO_sort.sam",4,0,0,0,
                0,0,0,0,2,2,0,false);
    }

    @Test
    public void testDuplicateMarkingSophisticated() {
        testHelper("duplicate_marking_sophisticated_1.sam","duplicate_marking_sophisticated_2.sam",6,0,0,0,
                0,0,0,0,4,0,0,true);
    }

    @Test
    public void testDuplicateMarkingSophisticatedCOSort() {
        testHelper("duplicate_marking_sophisticated_CO_sort_1.sam","duplicate_marking_sophisticated_CO_sort_2.sam",6,0,0,0,
                0,0,0,0,4,0,0,true);
    }

    @Test
    public void testDuplicateMarkingSophisticatedSumBaseQualsDiffer() {
        testHelper("duplicate_marking_sophisticated_sumBaseQualsDiff_1.sam","duplicate_marking_sophisticated_sumBaseQualsDiff_2.sam",6,0,0,0,
                0,0,0,0,4,4,0,false);
    }

    @Test
    public void testDuplicateMarkingSophisticatedRGDIff() {
        testHelper("duplicate_marking_sophisticated_RGDiff_1.sam","duplicate_marking_sophisticated_RGDiff_2.sam",6,0,0,0,
                0,0,0,0,4,4,0,false);
    }

}
