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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static picard.sam.markduplicates.MarkDuplicatesQueryNameTestSupport.*;

/**
 * Tests MarkDuplicates on queryname-sorted and query-grouped input: the header modes that select the query-grouped
 * path, secondary, supplementary and unmapped records taking the duplicate status of their query group, fragments
 * against pairs, libraries and read groups, ties, multiple inputs, empty input, the DT tag and duplicate removal
 * applied to whole query groups, spilling the read-end collections to disk, the exact metrics, and that the output
 * keeps the input order and every field other than the duplicate flag and the tags MarkDuplicates writes.
 *
 * Reads are 10 bases long and every base has the same quality, so with the default SUM_OF_BASE_QUALITIES strategy
 * a read scores 10 times its quality: 'I' (Q40) 400, 'D' (Q35) 350, '?' (Q30) 300 and '5' (Q20) 200, and a pair
 * scores the sum of its two reads.
 */
public class MarkDuplicatesQueryGroupedTest {

    private static final String NO_OPTICAL = "READ_NAME_REGEX=null";

    @DataProvider(name = "queryGroupedHeaders")
    public Object[][] queryGroupedHeaders() {
        return new Object[][]{
                {"SO:queryname", new String[0], SAMFileHeader.SortOrder.queryname, SAMFileHeader.GroupOrder.none},
                {"SO:queryname", new String[]{"ASSUME_SORT_ORDER=queryname"}, SAMFileHeader.SortOrder.unknown, SAMFileHeader.GroupOrder.query},
                {"SO:unsorted\tGO:query", new String[]{"ASSUME_SORT_ORDER=queryname"}, SAMFileHeader.SortOrder.unknown, SAMFileHeader.GroupOrder.query},
                {"SO:unknown\tGO:query", new String[]{"ASSUME_SORT_ORDER=queryname"}, SAMFileHeader.SortOrder.unknown, SAMFileHeader.GroupOrder.query},
                // A header claiming coordinate order is overridden for content that is only query grouped.
                {"SO:coordinate", new String[]{"ASSUME_SORT_ORDER=queryname"}, SAMFileHeader.SortOrder.unknown, SAMFileHeader.GroupOrder.query},
        };
    }

    /**
     * Pairs b1 and b2 share both 5' ends and b2 scores higher, so b1 is the duplicate; the fragment b3 starts at
     * a pair's 5' end and is a duplicate whatever its score; a1 is alone and comes first although it maps later.
     */
    @Test(dataProvider = "queryGroupedHeaders")
    public void testQueryGroupedHeaders(final String hdOrder, final String[] extraArgs,
                                        final SAMFileHeader.SortOrder expectedSortOrder,
                                        final SAMFileHeader.GroupOrder expectedGroupOrder) throws IOException {
        final SamText input = new SamText(hdOrder).readGroup("rgA", "libA")
                .pair("a1", "rgA", 500, false, 700, true, 'I')
                .pair("b1", "rgA", 100, false, 300, true, '5')
                .pair("b2", "rgA", 100, false, 300, true, 'I')
                .fragment("b3", "rgA", 100, false, '?');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL));

        Assert.assertEquals(result.header.getSortOrder(), expectedSortOrder);
        Assert.assertEquals(result.header.getGroupOrder(), expectedGroupOrder);
        assertDuplicateFlags(result, "..DD..D");
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        Assert.assertEquals(result.metrics.size(), 1);
        assertMetrics(result.metricsFor("libA"), 1, 3, 0, 0, 1, 1, 0);
    }

    /**
     * Secondary and supplementary records are not read ends, but in query-grouped mode they take the duplicate
     * flag, and the DT tag, of the primary reads in their query group, including a secondary record that comes
     * before its primaries. s1 is kept over s2 and s3. The fragment s5 is not a duplicate although supplementary
     * records of s1 and s4 align at its 5' end.
     */
    @Test
    public void testSecondaryAndSupplementaryFollowTheirQueryGroup() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA");
        input.record("s1", PAIRED | PROPER | FIRST | MATE_REVERSE, "chr1", 100, "10M", "chr1", 300, 'I', "RG:Z:rgA")
                .record("s1", PAIRED | FIRST | SECONDARY | MATE_REVERSE, "chr2", 1000, "10M", "chr1", 300, 'I', "RG:Z:rgA")
                .record("s1", PAIRED | PROPER | SECOND | REVERSE, "chr1", 300, "10M", "chr1", 100, 'I', "RG:Z:rgA")
                .record("s1", PAIRED | SECOND | SUPPLEMENTARY, "chr1", 5000, "10M", "chr1", 100, 'I', "RG:Z:rgA");
        input.record("s2", PAIRED | PROPER | FIRST | MATE_REVERSE, "chr1", 100, "10M", "chr1", 300, '5', "RG:Z:rgA")
                .record("s2", PAIRED | PROPER | SECOND | REVERSE, "chr1", 300, "10M", "chr1", 100, '5', "RG:Z:rgA")
                .record("s2", PAIRED | SECOND | SECONDARY | REVERSE, "chr2", 2000, "10M", "chr1", 100, '5', "RG:Z:rgA")
                .record("s2", PAIRED | FIRST | SUPPLEMENTARY | MATE_REVERSE, "chr1", 6000, "10M", "chr1", 300, '5', "RG:Z:rgA");
        input.record("s3", PAIRED | FIRST | SECONDARY | MATE_REVERSE, "chr2", 3000, "10M", "chr1", 300, '?', "RG:Z:rgA")
                .record("s3", PAIRED | PROPER | FIRST | MATE_REVERSE, "chr1", 100, "10M", "chr1", 300, '?', "RG:Z:rgA")
                .record("s3", PAIRED | PROPER | SECOND | REVERSE, "chr1", 300, "10M", "chr1", 100, '?', "RG:Z:rgA");
        input.pair("s4", "rgA", 700, false, 900, true, 'I')
                .record("s4", PAIRED | FIRST | SUPPLEMENTARY | MATE_REVERSE, "chr1", 5000, "10M", "chr1", 900, 'I', "RG:Z:rgA");
        input.fragment("s5", "rgA", 5000, false, '5');

        final Result result = run(input.write(), NO_OPTICAL, "TAGGING_POLICY=All");

        assertDuplicateFlags(result, "....DDDDDDD....");
        for (final SAMRecord rec : result.records) {
            Assert.assertEquals(rec.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG),
                    rec.getDuplicateReadFlag() ? MarkDuplicates.DUPLICATE_TYPE_LIBRARY : null, rec.getSAMString());
        }
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        assertMetrics(result.metricsFor("libA"), 1, 4, 6, 0, 0, 2, 0);
    }

    /**
     * A pair with one unmapped read is judged by its mapped read as a fragment, and its unmapped read, whether it
     * comes first or second in the group, takes the mapped read's duplicate flag. The unplaced unmapped pair u0
     * comes first and is never a duplicate, and reading continues past it; u1 is kept as the highest scoring
     * fragment at 100F; u2 and the fragment u3 are duplicates.
     */
    private static SamText unmappedMatesInput() {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA");
        input.record("u0", PAIRED | FIRST | UNMAPPED | MATE_UNMAPPED, "*", 0, "*", "*", 0, 'I', "RG:Z:rgA")
                .record("u0", PAIRED | SECOND | UNMAPPED | MATE_UNMAPPED, "*", 0, "*", "*", 0, 'I', "RG:Z:rgA");
        return input.pairWithUnmappedMate("u1", "rgA", 100, false, false, 'I')
                .pairWithUnmappedMate("u2", "rgA", 100, false, true, '5')
                .fragment("u3", "rgA", 100, false, '?');
    }

    @Test
    public void testUnmappedMatesFollowTheirQueryGroup() throws IOException {
        final SamText input = unmappedMatesInput();

        final Result result = run(input.write(), NO_OPTICAL, "TAGGING_POLICY=All");

        assertDuplicateFlags(result, "....DDD");
        Assert.assertEquals(result.tagValues(MarkDuplicates.DUPLICATE_TYPE_TAG),
                Arrays.asList(null, null, null, null, "LB", "LB", "LB"));
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        assertMetrics(result.metricsFor("libA"), 3, 0, 0, 4, 2, 0, 0);
    }

    /** REMOVE_DUPLICATES drops the unmapped mate of a duplicate fragment along with it. */
    @Test
    public void testRemoveDuplicatesDropsUnmappedMates() throws IOException {
        final SamText input = unmappedMatesInput();

        final Result result = run(input.write(), NO_OPTICAL, "REMOVE_DUPLICATES=true");

        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines().subList(0, 4), result.recordLines);
        assertDuplicateFlags(result, "....");
        assertMetrics(result.metricsFor("libA"), 3, 0, 0, 4, 2, 0, 0);
    }

    /**
     * With a mapped pair u0 ending at 100F, every fragment there is a duplicate whatever its score, so both
     * records of u1 and u2 are marked and u0 is kept.
     */
    @Test
    public void testUnmappedMatesOfFragmentsBehindAPair() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("u0", "rgA", 100, false, 300, true, '5')
                .pairWithUnmappedMate("u1", "rgA", 100, false, false, 'I')
                .pairWithUnmappedMate("u2", "rgA", 100, false, true, 'I')
                .fragment("u3", "rgA", 100, false, '?');
        input.record("u4", PAIRED | FIRST | UNMAPPED | MATE_UNMAPPED, "*", 0, "*", "*", 0, 'I', "RG:Z:rgA")
                .record("u4", PAIRED | SECOND | UNMAPPED | MATE_UNMAPPED, "*", 0, "*", "*", 0, 'I', "RG:Z:rgA");

        final Result result = run(input.write(), NO_OPTICAL);

        assertDuplicateFlags(result, "..DDDDD..");
        assertMetrics(result.metricsFor("libA"), 3, 1, 0, 4, 3, 0, 0);
    }

    @DataProvider(name = "inMemoryAndSpilled")
    public Object[][] inMemoryAndSpilled() {
        return new Object[][]{{new String[0]}, {spillArgs()}};
    }

    /**
     * Pairs are duplicates when both 5' ends and the pair orientation match, whichever read is read one: p2 is a
     * duplicate of p1, while p3 (other read-two end) and p4 (FF) are not. A fragment is a duplicate when a pair
     * has an end at its 5' end and strand (f1 at 100F, f2 at 409R), or when a higher scoring fragment does (f4
     * against f5); f3 at 109R and f6 at 2009R are alone.
     */
    @Test(dataProvider = "inMemoryAndSpilled")
    public void testFragmentsAndPairs(final String[] extraArgs) throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .fragment("f1", "rgA", 100, false, 'I')
                .fragment("f2", "rgA", 400, true, 'I')
                .fragment("f3", "rgA", 100, true, 'I')
                .fragment("f4", "rgA", 2000, false, '5')
                .fragment("f5", "rgA", 2000, false, 'I')
                .fragment("f6", "rgA", 2000, true, '5')
                .pair("p1", "rgA", 100, false, 300, true, 'I')
                .pair("p2", "rgA", 300, true, 100, false, '5')
                .pair("p3", "rgA", 100, false, 400, true, '5')
                .pair("p4", "rgA", 100, false, 300, false, '5');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL));

        assertDuplicateFlags(result, "DD.D....DD....");
        assertMetrics(result.metricsFor("libA"), 6, 4, 0, 0, 3, 1, 0);
    }

    /**
     * Read groups rg1 and rg2 share libA, so m2 is a duplicate of m1 across read groups; libB (rg3) is judged
     * separately; rg4 has no LB and forms "Unknown Library", where the fragment m5 is kept although pairs of other
     * libraries end at 100F.
     */
    @Test(dataProvider = "inMemoryAndSpilled")
    public void testLibrariesAndReadGroups(final String[] extraArgs) throws IOException {
        final SamText input = new SamText("SO:queryname")
                .readGroup("rg1", "libA").readGroup("rg2", "libA").readGroup("rg3", "libB").readGroup("rg4", null)
                .pair("m1", "rg1", 100, false, 300, true, 'I')
                .pair("m2", "rg2", 100, false, 300, true, '5')
                .pair("m3", "rg3", 100, false, 300, true, '5')
                .pair("m4", "rg3", 100, false, 300, true, '?')
                .fragment("m5", "rg4", 100, false, 'I')
                .fragment("m6", "rg4", 100, false, '5');

        final Result result = run(input.write(), args(extraArgs, NO_OPTICAL));

        assertDuplicateFlags(result, "..DDDD...D");
        Assert.assertEquals(result.metrics.size(), 3);
        assertMetrics(result.metricsFor("libA"), 0, 2, 0, 0, 0, 1, 0);
        assertMetrics(result.metricsFor("libB"), 0, 2, 0, 0, 0, 1, 0);
        assertMetrics(result.metricsFor("Unknown Library"), 2, 0, 0, 0, 1, 0, 0);
    }

    /**
     * A pair scores the sum of its two reads. Read one / read two score 400/200 for pA, 350/350 for pB and 200/400
     * for pC, so pB (700) is kept; scoring by read one alone would keep pA and by read two alone pC.
     */
    @Test
    public void testPairScoreIsTheSumOfBothReads() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("pA", "rgA", 100, false, 300, true, 'I', '5')
                .pair("pB", "rgA", 100, false, 300, true, 'D', 'D')
                .pair("pC", "rgA", 100, false, 300, true, '5', 'I');

        final Result result = run(input.write(), NO_OPTICAL);

        assertDuplicateFlags(result, "DD..DD");
    }

    /**
     * Pairs with mates on different contigs are matched by read group and name across contigs: c2 (read one on chr2)
     * has the same ends as c1 (read one on chr1) and is a duplicate, while c3, whose read two is at the same
     * coordinate on chr1 instead of chr2, is not.
     */
    @Test
    public void testPairsAcrossContigs() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .record("c1", PAIRED | FIRST | MATE_REVERSE, "chr1", 100, "10M", "chr2", 500, 'I', "MC:Z:10M", "RG:Z:rgA")
                .record("c1", PAIRED | SECOND | REVERSE, "chr2", 500, "10M", "chr1", 100, 'I', "MC:Z:10M", "RG:Z:rgA")
                .record("c2", PAIRED | FIRST | REVERSE, "chr2", 500, "10M", "chr1", 100, '5', "MC:Z:10M", "RG:Z:rgA")
                .record("c2", PAIRED | SECOND | MATE_REVERSE, "chr1", 100, "10M", "chr2", 500, '5', "MC:Z:10M", "RG:Z:rgA")
                .pair("c3", "rgA", 100, false, 500, true, '5');

        final Result result = run(input.write(), NO_OPTICAL);

        assertDuplicateFlags(result, "..DD..");
        assertMetrics(result.metricsFor("libA"), 0, 3, 0, 0, 0, 1, 0);
    }

    /**
     * When both 5' ends of a pair are at one coordinate on opposite strands, the pair is FR whichever read is
     * reverse: o1 (read one reverse at 100, 5' end 109) and o2 (read one forward at 109) are duplicates, while o3
     * (both forward at 109) is FF and alone.
     */
    @Test
    public void testPairWithBothEndsAtOneCoordinate() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("o1", "rgA", 100, true, 109, false, 'I')
                .pair("o2", "rgA", 109, false, 100, true, '5')
                .pair("o3", "rgA", 109, false, 109, false, '5');

        final Result result = run(input.write(), NO_OPTICAL);

        assertDuplicateFlags(result, "..DD..");
    }

    @DataProvider(name = "tiedGroupOrders")
    public Object[][] tiedGroupOrders() {
        return new Object[][]{
                {Arrays.asList("t1", "t2", "t3"), Arrays.asList("t4", "t5")},
                {Arrays.asList("t3", "t1", "t2"), Arrays.asList("t5", "t4")},
        };
    }

    /**
     * With equal scores, the pair and the fragment kept are those whose query group comes first in the file. The
     * documentation does not specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "tiedGroupOrders")
    public void testTiedScoresKeepTheFirstQueryGroup(final List<String> pairNames, final List<String> fragmentNames)
            throws IOException {
        final SamText input = new SamText("SO:unsorted\tGO:query").readGroup("rgA", "libA");
        pairNames.forEach(name -> input.pair(name, "rgA", 100, false, 300, true, 'I'));
        fragmentNames.forEach(name -> input.fragment(name, "rgA", 2000, false, 'I'));

        final Result result = run(input.write(), NO_OPTICAL, "ASSUME_SORT_ORDER=queryname");

        assertDuplicateFlags(result, "..DDDD.D");
        Assert.assertEquals(result.records.get(0).getReadName(), pairNames.get(0));
        Assert.assertEquals(result.records.get(6).getReadName(), fragmentNames.get(0));
    }

    /**
     * Groups out of name order stay in input order; a stale duplicate flag and DT tag on the kept pair are
     * cleared; every other field and tag is unchanged.
     */
    @Test
    public void testOutputKeepsOrderAndFields() throws IOException {
        final SamText input = new SamText("SO:unsorted\tGO:query").readGroup("rgA", "libA");
        input.record("z1", PAIRED | PROPER | FIRST | MATE_REVERSE | DUPLICATE, "chr1", 100, "10M", "chr1", 300, 'I',
                        "MC:Z:10M", "RG:Z:rgA", "DT:Z:LB", "XA:Z:foo", "XI:i:7")
                .record("z1", PAIRED | PROPER | SECOND | REVERSE, "chr1", 300, "10M", "chr1", 100, 'I',
                        "MC:Z:10M", "RG:Z:rgA", "DT:Z:SQ", "XI:i:-3");
        input.pair("a1", "rgA", 100, false, 300, true, '5', "XA:Z:bar")
                .record("a1", PAIRED | FIRST | SECONDARY | MATE_REVERSE, "chr2", 4000, "4S6M", "chr1", 300, '5',
                        "RG:Z:rgA", "NM:i:1");
        input.fragment("m1", "rgA", 500, false, '?', "XZ:Z:a b");
        input.record("n1", PAIRED | FIRST | UNMAPPED | MATE_UNMAPPED, "*", 0, "*", "*", 0, 'I', "RG:Z:rgA")
                .record("n1", PAIRED | SECOND | UNMAPPED | MATE_UNMAPPED, "*", 0, "*", "*", 0, 'I', "RG:Z:rgA");

        final Result result = run(input.write(), NO_OPTICAL, "ASSUME_SORT_ORDER=queryname");

        assertDuplicateFlags(result, "..DDD...");
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        Assert.assertNull(result.records.get(0).getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG));
        Assert.assertNull(result.records.get(1).getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG));
    }

    @DataProvider(name = "opticalTaggingAndRemoval")
    public Object[][] opticalTaggingAndRemoval() {
        return new Object[][]{
                {"All", new String[0], "..SSSLLL..SS", new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}},
                {"OpticalOnly", new String[0], "..SSS.....SS", new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}},
                {"All", new String[]{"REMOVE_SEQUENCING_DUPLICATES=true"}, "..LLL..", new int[]{0, 1, 5, 6, 7, 8, 9}},
                {"All", new String[]{"REMOVE_DUPLICATES=true"}, "....", new int[]{0, 1, 8, 9}},
        };
    }

    /**
     * In query-grouped mode an optical duplicate pair's DT tag of SQ, and its removal by
     * REMOVE_SEQUENCING_DUPLICATES or REMOVE_DUPLICATES, cover its secondary and supplementary records too. Read
     * names carry tile:x:y, so RUN:1:1:100:150 lies 50 pixels from the kept RUN:1:1:100:100 (optical), and
     * RUN:1:1:5000:5000 lies far from it (library duplicate); RUN:1:1:6000:150 is optical against RUN:1:1:6000:100
     * in a second set.
     *
     * @param expectedTags one character per output record: 'S' for DT=SQ, 'L' for DT=LB, '.' for no DT tag
     * @param survivors    the input records expected in the output, in order
     */
    @Test(dataProvider = "opticalTaggingAndRemoval")
    public void testOpticalDuplicatesTagAndRemoveWholeQueryGroups(final String taggingPolicy, final String[] extraArgs,
                                                                  final String expectedTags, final int[] survivors)
            throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("RUN:1:1:100:100", "rgA", 100, false, 300, true, 'I')
                .pair("RUN:1:1:100:150", "rgA", 100, false, 300, true, '5')
                .record("RUN:1:1:100:150", PAIRED | FIRST | SECONDARY | MATE_REVERSE, "chr2", 1000, "10M", "chr1", 300, '5', "RG:Z:rgA")
                .pair("RUN:1:1:5000:5000", "rgA", 100, false, 300, true, '5')
                .record("RUN:1:1:5000:5000", PAIRED | SECOND | SUPPLEMENTARY, "chr1", 8000, "10M", "chr1", 100, '5', "RG:Z:rgA")
                .pair("RUN:1:1:6000:100", "rgA", 1000, false, 1200, true, 'I')
                .pair("RUN:1:1:6000:150", "rgA", 1000, false, 1200, true, '5');

        final Result result = run(input.write(), args(extraArgs, "TAGGING_POLICY=" + taggingPolicy));

        final List<String> expectedLines = new ArrayList<>();
        for (final int index : survivors) {
            expectedLines.add(input.recordLines().get(index));
        }
        assertOnlyDuplicateFlagAndTagsChanged(expectedLines, result.recordLines);

        final StringBuilder tags = new StringBuilder();
        for (final SAMRecord rec : result.records) {
            final Object dt = rec.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG);
            tags.append(dt == null ? '.' : MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(dt) ? 'S' : 'L');
            Assert.assertEquals(rec.getDuplicateReadFlag(),
                    !rec.getReadName().equals("RUN:1:1:100:100") && !rec.getReadName().equals("RUN:1:1:6000:100"));
        }
        Assert.assertEquals(tags.toString(), expectedTags);
        assertMetrics(result.metricsFor("libA"), 0, 5, 2, 0, 0, 3, 2);
    }

    @DataProvider(name = "multipleInputs")
    public Object[][] multipleInputs() {
        return new Object[][]{{new String[0]}, {new String[]{"ASSUME_SORT_ORDER=queryname"}}};
    }

    /**
     * Several queryname-sorted inputs are merged by name, and the pair kept is the best across all inputs; the
     * fragment q4 ends at 309R, a pair end, so it is a duplicate.
     */
    @Test(dataProvider = "multipleInputs")
    public void testMultipleInputsAreMergedByQueryName(final String[] extraArgs) throws IOException {
        final File first = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("q1", "rgA", 100, false, 300, true, '5')
                .pair("q3", "rgA", 100, false, 300, true, '?')
                .write();
        final File second = new SamText("SO:queryname").readGroup("rgA", "libA")
                .pair("q2", "rgA", 100, false, 300, true, 'I')
                .fragment("q4", "rgA", 300, true, 'I')
                .write();

        final Result result = run(Arrays.asList(first, second), args(extraArgs, NO_OPTICAL));

        final List<String> names = new ArrayList<>();
        result.records.forEach(r -> names.add(r.getReadName()));
        Assert.assertEquals(names, Arrays.asList("q1", "q1", "q2", "q2", "q3", "q3", "q4"));
        assertDuplicateFlags(result, "DD..DDD");
        assertMetrics(result.metricsFor("libA"), 1, 3, 0, 0, 1, 2, 0);
    }

    /** A header-only input writes no records and zeroed metrics for each library in the header. */
    @Test
    public void testEmptyInput() throws IOException {
        final SamText input = new SamText("SO:queryname").readGroup("rgA", "libA").readGroup("rgB", "libB");

        final Result result = run(input.write(), NO_OPTICAL, "ASSUME_SORT_ORDER=queryname");

        Assert.assertTrue(result.records.isEmpty());
        Assert.assertEquals(result.metrics.size(), 2);
        assertMetrics(result.metricsFor("libA"), 0, 0, 0, 0, 0, 0, 0);
        assertMetrics(result.metricsFor("libB"), 0, 0, 0, 0, 0, 0, 0);
    }

    @DataProvider(name = "singleGroups")
    public Object[][] singleGroups() {
        return new Object[][]{
                {new SamText("SO:queryname").readGroup("rgA", "libA").pair("s1", "rgA", 100, false, 300, true, 'I'), "..", 0, 1},
                {new SamText("SO:queryname").readGroup("rgA", "libA").fragment("s1", "rgA", 100, false, 'I'), ".", 1, 0},
                {new SamText("SO:queryname").readGroup("rgA", "libA").pairWithUnmappedMate("s1", "rgA", 100, true, false, 'I'), "..", 1, 0},
        };
    }

    @Test(dataProvider = "singleGroups")
    public void testSingleQueryGroup(final SamText input, final String expectedFlags, final long unpaired,
                                     final long pairs) throws IOException {
        final Result result = run(input.write(), NO_OPTICAL);

        assertDuplicateFlags(result, expectedFlags);
        assertMetrics(result.metricsFor("libA"), unpaired, pairs, 0, expectedFlags.length() - unpaired - 2 * pairs, 0, 0, 0);
    }
}
