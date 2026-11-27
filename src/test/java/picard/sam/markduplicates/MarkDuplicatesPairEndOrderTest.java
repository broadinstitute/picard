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

import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * Tests {@link MarkDuplicates} on coordinate-sorted pairs whose read with the lower 5' coordinate comes later in the
 * file than its mate. Pair "s" has its reverse read aligned at chr1:100 (50M, 5' end 149) and its forward read at
 * chr1:120 (50M, 5' end 120), so its forward read is the second of its records in the file although its 5' end is
 * the lower one. Pair "n" has the same 5' ends with its forward read first: forward at chr1:120 (50M) and reverse at
 * chr1:120 (20S30M, 5' end 149). The two pairs are duplicates of each other, and the duplicate is marked on both of
 * its records.
 *
 * The file holds s's reverse read (index 0), n's forward read (1), n's reverse read (2) and s's forward read (3).
 */
public class MarkDuplicatesPairEndOrderTest {
    private static final int S_FORWARD_INDEX = 3;
    private static final int N_FORWARD_INDEX = 1;

    /** Adds n, then s; s's first of pair is its reverse read unless sFirstIsForward. */
    private static MarkDuplicatesEndToEndRunner input(final String n, final String s, final boolean sFirstIsForward,
                                                      final int nQuality, final int sQuality) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        input.pair(n, "rg", 0, 120, false, "50M", 0, 120, true, "20S30M", nQuality);
        if (sFirstIsForward) {
            input.pair(s, "rg", 0, 120, false, "50M", 0, 100, true, "50M", sQuality);
        } else {
            input.pair(s, "rg", 0, 100, true, "50M", 0, 120, false, "50M", sQuality);
        }
        return input;
    }

    private static List<String> bothEnds(final String name) {
        return List.of(name + "/1", name + "/2");
    }

    @DataProvider(name = "scores")
    public Object[][] scores() {
        return new Object[][]{
                // s's first of pair is forward, quality of n's reads, quality of s's reads, whether s is the duplicate
                {false, 20, 15, true},
                {true, 20, 15, true},
                {false, 15, 20, false},
                {true, 15, 20, false},
                {false, 20, 20, true},
                {true, 20, 20, true},
        };
    }

    /**
     * The pair with the higher sum of base qualities is kept and both records of the other are marked, with DT:Z:LB,
     * DS:i:2 and the DI of the kept pair on all four records.
     *
     * The documentation does not specify this; the expected value follows the current implementation. With equal
     * scores n is kept, the pair whose lower-5' read comes first in the file (index 1 against 3). The DI of a pair is
     * the index of its read with the lower 5' coordinate, which is the later record of s.
     */
    @Test(dataProvider = "scores")
    public void testPairWithLaterLowerEndIsMarkedOnBothRecords(final boolean sFirstIsForward, final int nQuality,
                                                               final int sQuality, final boolean sIsDuplicate) {
        final MarkDuplicatesEndToEndRunner input = input("n", "s", sFirstIsForward, nQuality, sQuality);
        Assert.assertEquals(input.inputRecords().stream().map(MarkDuplicatesEndToEndRunner::key).collect(Collectors.toList()),
                sFirstIsForward ? List.of("s/2", "n/1", "n/2", "s/1") : List.of("s/1", "n/1", "n/2", "s/2"));

        final MarkDuplicatesEndToEndRunner.Result result = input.run("TAGGING_POLICY=All", "TAG_DUPLICATE_SET_MEMBERS=true");

        final String duplicate = sIsDuplicate ? "s" : "n";
        result.assertDuplicates(bothEnds(duplicate).toArray(new String[0]));
        final Map<String, Object> expectedDt = new TreeMap<>();
        bothEnds(duplicate).forEach(key -> expectedDt.put(key, MarkDuplicates.DUPLICATE_TYPE_LIBRARY));
        final Map<String, Object> dt = new TreeMap<>();
        final List<Object> di = new ArrayList<>();
        final List<Object> ds = new ArrayList<>();
        for (final SAMRecord rec : result.records) {
            if (rec.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG) != null) {
                dt.put(MarkDuplicatesEndToEndRunner.key(rec), rec.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG));
            }
            di.add(rec.getAttribute(MarkDuplicates.DUPLICATE_SET_INDEX_TAG));
            ds.add(rec.getAttribute(MarkDuplicates.DUPLICATE_SET_SIZE_TAG));
        }
        Assert.assertEquals(dt, expectedDt);
        final int representativeIndex = sIsDuplicate ? N_FORWARD_INDEX : S_FORWARD_INDEX;
        Assert.assertEquals(di, List.of(representativeIndex, representativeIndex, representativeIndex, representativeIndex));
        Assert.assertEquals(ds, List.of(2, 2, 2, 2));
        MarkDuplicatesEndToEndRunner.assertMetrics(result.onlyMetrics(), "lib", 0, 2, 0, 0, 0, 1, 0);
    }

    @DataProvider(name = "firstOfPairStrands")
    public Object[][] firstOfPairStrands() {
        return new Object[][]{{true}, {false}};
    }

    /**
     * With read names 50 pixels apart in one tile, s is an optical duplicate of n when the first of pair of both is
     * forward, and then both of its records carry DT:Z:SQ and are removed by REMOVE_SEQUENCING_DUPLICATES. When s's
     * first of pair is its reverse read the pairs are library duplicates only.
     */
    @Test(dataProvider = "firstOfPairStrands")
    public void testPairWithLaterLowerEndAsOpticalDuplicate(final boolean sFirstIsForward) {
        final String n = "n:1:1101:1000:1000";
        final String s = "s:1:1101:1050:1000";
        final MarkDuplicatesEndToEndRunner input = input(n, s, sFirstIsForward, 20, 15);

        final MarkDuplicatesEndToEndRunner.Result tagged = input.run("TAGGING_POLICY=All");
        tagged.assertDuplicates(bothEnds(s).toArray(new String[0]));
        final String expectedType = sFirstIsForward ? MarkDuplicates.DUPLICATE_TYPE_SEQUENCING : MarkDuplicates.DUPLICATE_TYPE_LIBRARY;
        for (final String key : bothEnds(s)) {
            Assert.assertEquals(tagged.record(key).getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG), expectedType, key);
        }
        for (final String key : bothEnds(n)) {
            Assert.assertNull(tagged.record(key).getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG), key);
        }
        MarkDuplicatesEndToEndRunner.assertMetrics(tagged.onlyMetrics(), "lib", 0, 2, 0, 0, 0, 1, sFirstIsForward ? 1 : 0);

        final MarkDuplicatesEndToEndRunner.Result removed = input.run("REMOVE_SEQUENCING_DUPLICATES=true");
        final List<String> expectedKeys = input.inputRecords().stream().map(MarkDuplicatesEndToEndRunner::key)
                .filter(key -> !sFirstIsForward || !key.startsWith(s + "/")).collect(Collectors.toList());
        Assert.assertEquals(removed.keys(), expectedKeys);
    }
}
