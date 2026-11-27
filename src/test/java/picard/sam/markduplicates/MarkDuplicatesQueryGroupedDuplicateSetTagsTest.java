/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import static picard.sam.markduplicates.MarkDuplicatesQueryNameTestSupport.*;

/**
 * Tests TAG_DUPLICATE_SET_MEMBERS on query-grouped input. Every primary mapped read of a pair in a duplicate set
 * gets DS, the number of pairs in the set, and DI, the index in the file of the representative pair's read one,
 * the representative being the highest scoring pair and, among equal scores, the one whose group comes first.
 * Every query group here starts with its primary read one. Secondary records and pairs in no set get neither tag.
 */
public class MarkDuplicatesQueryGroupedDuplicateSetTagsTest {

    private static final String NO_OPTICAL = "READ_NAME_REGEX=null";

    @DataProvider(name = "inMemoryAndSpilled")
    public Object[][] inMemoryAndSpilled() {
        return new Object[][]{
                {"SO:queryname", new String[0]},
                {"SO:queryname", spillArgs()},
                {"SO:unsorted\tGO:query", new String[]{"ASSUME_SORT_ORDER=queryname"}},
        };
    }

    /**
     * Set A is d1 and d2 at 100F/309R, represented by d2 (index 2) although d1 comes first. Set B is d3, d4 and d5 at
     * 700F/909R, represented by d3 (index 5). d8 is in no set. Set C is d6 and d7 with equal scores, represented by
     * the first, d6 (index 11), and d2's secondary record (index 3) gets neither tag. The documentation does not
     * specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "inMemoryAndSpilled")
    public void testDuplicateSetTags(final String hdOrder, final String[] extraArgs) throws IOException {
        final SamText input = new SamText(hdOrder).readGroup("rgA", "libA")
                .pair("d1", "rgA", 100, false, 300, true, '5');
        input.record("d2", PAIRED | PROPER | FIRST | MATE_REVERSE, "chr1", 100, "10M", "chr1", 300, 'I', "MC:Z:10M", "RG:Z:rgA")
                .record("d2", PAIRED | FIRST | SECONDARY | MATE_REVERSE, "chr2", 500, "10M", "chr1", 300, 'I', "RG:Z:rgA")
                .record("d2", PAIRED | PROPER | SECOND | REVERSE, "chr1", 300, "10M", "chr1", 100, 'I', "MC:Z:10M", "RG:Z:rgA");
        input.pair("d3", "rgA", 700, false, 900, true, '?')
                .pair("d4", "rgA", 700, false, 900, true, '5')
                .pair("d5", "rgA", 700, false, 900, true, '5')
                .pair("d6", "rgA", 1500, false, 1700, true, 'I')
                .pair("d7", "rgA", 1500, false, 1700, true, 'I')
                .pair("d8", "rgA", 2500, false, 2700, true, 'I');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL, "TAG_DUPLICATE_SET_MEMBERS=true"));

        assertDuplicateFlags(result, "DD.....DDDD..DD..");
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_INDEX_TAG),
                Arrays.asList(2, 2, 2, null, 2, 5, 5, 5, 5, 5, 5, 11, 11, 11, 11, null, null));
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_SIZE_TAG),
                Arrays.asList(2, 2, 2, null, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, null, null));
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
    }

    /** Without TAG_DUPLICATE_SET_MEMBERS no record gets DI or DS. */
    @Test
    public void testNoDuplicateSetTagsByDefault() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("d1", "rgA", 100, false, 300, true, '5')
                .pair("d2", "rgA", 100, false, 300, true, 'I');

        final Result result = run(input.write(), NO_OPTICAL);

        assertDuplicateFlags(result, "DD..");
        final List<Object> none = Arrays.asList(null, null, null, null);
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_INDEX_TAG), none);
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_SIZE_TAG), none);
    }

    /**
     * Pairs at the same position in different libraries form separate sets: {e1, e2} in libA represented by e1
     * (index 0), and {e3, e4, e5} in libB represented by e4 (index 6).
     */
    @Test(dataProvider = "inMemoryAndSpilled")
    public void testDuplicateSetsAreSplitByLibrary(final String hdOrder, final String[] extraArgs) throws IOException {
        final SamText input = new SamText(hdOrder).readGroup("rgA", "libA").readGroup("rgB", "libB")
                .pair("e1", "rgA", 100, false, 300, true, 'I')
                .pair("e2", "rgA", 100, false, 300, true, '5')
                .pair("e3", "rgB", 100, false, 300, true, '5')
                .pair("e4", "rgB", 100, false, 300, true, '?')
                .pair("e5", "rgB", 100, false, 300, true, '5');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL, "TAG_DUPLICATE_SET_MEMBERS=true"));

        assertDuplicateFlags(result, "..DDDD..DD");
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_INDEX_TAG),
                Arrays.asList(0, 0, 0, 0, 6, 6, 6, 6, 6, 6));
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_SET_SIZE_TAG),
                Arrays.asList(2, 2, 2, 2, 3, 3, 3, 3, 3, 3));
    }
}
