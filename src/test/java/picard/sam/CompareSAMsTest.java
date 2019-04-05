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
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;

public class CompareSAMsTest extends CommandLineProgramTest {
    private static final File TEST_FILES_DIR = new File("testdata/picard/sam/CompareSAMs");

    public String getCommandLineProgramName() {
        return CompareSAMs.class.getSimpleName();
    }

    @DataProvider(name="compareSams")
    public Object[][] compareSamsTestData() {
        return new Object[][] {
                { "genomic_sorted.sam", "unsorted.sam", false },
                { "genomic_sorted.sam", "chr21.sam", false },
                { "genomic_sorted.sam", "bigger_seq_dict.sam", false },
                { "bigger_seq_dict.sam", "bigger_seq_dict.sam", true },
                { "genomic_sorted.sam", "genomic_sorted.sam", true },
                { "genomic_sorted.sam", "has_non_primary.sam", true },
                { "genomic_sorted_5.sam", "genomic_sorted_5_plus.sam", false },
                { "group_same_coord.sam", "group_same_coord_diff_order.sam", false },
                { "genomic_sorted_same_position.sam", "genomic_sorted_same_position.sam", true },
                { "group_same_coord.sam", "diff_coords.sam", false },
                { "genomic_sorted.sam", "unmapped_first.sam", false },
                { "genomic_sorted.sam", "unmapped_second.sam", false },
                { "unmapped_first.sam", "unmapped_second.sam", false },
                { "unmapped_first.sam", "unmapped_first.sam", true },
                { "genomic_sorted.sam", "genomic_sorted_sam_v1.6.sam", false },
                { "unsorted.sam", "unsorted.sam", true },
                { "unsorted.sam", "unsorted2.sam", false},
                {"dup_first.sam", "dup_second.sam", true},
                {"dup_first.sam", "dup_third.sam", false},
                {"dup_first.sam", "dup_fourth.sam", false},
                {"dup_first.sam", "dup_fifth.sam", false},
                {"dup_first.sam", "dup_second.sam", "DUPLICATE_STRINGENCY=STRICT", false},
                {"dup_first_queryname.sam", "dup_second_queryname.sam", true},
                {"dup_first_queryname.sam", "dup_third_queryname.sam", false},
                {"dup_first_queryname.sam", "dup_fourth_queryname.sam", false},
                {"dup_first_queryname.sam", "dup_fifth_queryname.sam", false},
                {"dup_first_queryname.sam", "dup_second_queryname.sam", false}

        };
    }

    @Test(dataProvider="compareSams")
    public void testCompareSAMs(final String f1, final String f2, final boolean areEqual) {
        final String[] samFiles = {
                new File(TEST_FILES_DIR, f1).getAbsolutePath(),
                new File(TEST_FILES_DIR, f2).getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(samFiles) == 0, areEqual);

        final String[] samFilesReversed = {
                new File(TEST_FILES_DIR, f2).getAbsolutePath(),
                new File(TEST_FILES_DIR, f1).getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(samFilesReversed) == 0, areEqual);
    }

}
