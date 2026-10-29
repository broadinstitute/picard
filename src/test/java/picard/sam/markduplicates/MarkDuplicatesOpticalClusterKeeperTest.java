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
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * Tests which pair of an optical duplicate cluster that does not contain the kept pair of its duplicate set is left
 * without DT:Z:SQ, and kept by REMOVE_SEQUENCING_DUPLICATES, when {@link MarkDuplicates} runs on coordinate-sorted
 * input. Both paths of {@link OpticalDuplicateFinder} are covered: sets of three pairs, which it compares pair by pair,
 * and sets of four or more, which it clusters through a graph, the latter also above
 * {@link OpticalDuplicateFinder#DEFAULT_BIG_DUPLICATE_SET_SIZE} pairs, where it logs its progress.
 *
 * The documentation does not specify this; the expected value follows the current implementation. Of two pairs within
 * the pixel distance of each other, the one with the lower x, then the lower y, keeps no DT:Z:SQ; of two pairs at the
 * same x and y, the one earlier in the file. Here that pair is never chosen by score, read name or file order alone:
 * the pair written first in the file ("a") has the higher score and the smaller name.
 *
 * Every set is a set of pairs with 5' ends chr1:100 forward and chr1:300 reverse. The kept pair K scores highest and
 * lies far from the cluster in the same tile; the other pairs ("f0", "f1", ...) lie more than 100 pixels from every
 * other pair and are library duplicates.
 */
public class MarkDuplicatesOpticalClusterKeeperTest {
    private static final int TILE = 1101;
    private static final String KEEPER = name("K", 5000, 5000);

    private static String name(final String prefix, final int x, final int y) {
        return prefix + ":1:" + TILE + ":" + x + ":" + y;
    }

    private static final int FAST_PATH_SET_SIZE = 3;
    private static final int GRAPH_PATH_SET_SIZE = 4;
    private static final int LOGGED_SET_SIZE = OpticalDuplicateFinder.DEFAULT_BIG_DUPLICATE_SET_SIZE + 1;

    @DataProvider(name = "clusters")
    public Object[][] clusters() {
        final int[][][] coordinates = {
                // {x, y} of "a" (first in the file), {x, y} of "z", and which of the two keeps no DT:Z:SQ (0 = a, 1 = z)
                {{1050, 1000}, {999, 1000}, {1}},
                {{1000, 1000}, {1000, 950}, {1}},
                {{1050, 1000}, {999, 1080}, {1}},
                {{999, 1000}, {1050, 1000}, {0}},
                {{1000, 1000}, {1000, 1000}, {0}},
        };
        final List<Object[]> rows = new ArrayList<>();
        for (final int setSize : new int[]{FAST_PATH_SET_SIZE, GRAPH_PATH_SET_SIZE, LOGGED_SET_SIZE}) {
            for (final int[][] c : coordinates) {
                rows.add(new Object[]{setSize, c[0][0], c[0][1], c[1][0], c[1][1], c[2][0] == 0});
            }
        }
        return rows.toArray(new Object[0][]);
    }

    @Test(dataProvider = "clusters")
    public void testWhichPairOfAClusterWithoutTheKeeperIsTagged(final int setSize, final int ax, final int ay,
                                                                final int zx, final int zy, final boolean aKeepsNoTag) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final String a = name("a", ax, ay);
        final String z = name("z", zx, zy);
        addPair(input, KEEPER, 30);
        addPair(input, a, 20);
        addPair(input, z, 15);
        final List<String> others = new ArrayList<>();
        for (int i = 0; i < setSize - 3; i++) {
            others.add(name("f" + i, 20000 + 101 * i, 1000));
            addPair(input, others.get(i), 20);
        }
        final String untagged = aKeepsNoTag ? a : z;
        final String tagged = aKeepsNoTag ? z : a;

        final MarkDuplicatesEndToEndRunner.Result tagging = input.run("TAGGING_POLICY=All");
        final Map<String, Object> expectedTags = new TreeMap<>();
        for (final String name : others) {
            putBothEnds(expectedTags, name, MarkDuplicates.DUPLICATE_TYPE_LIBRARY);
        }
        putBothEnds(expectedTags, untagged, MarkDuplicates.DUPLICATE_TYPE_LIBRARY);
        putBothEnds(expectedTags, tagged, MarkDuplicates.DUPLICATE_TYPE_SEQUENCING);
        Assert.assertEquals(duplicateTypeTags(tagging.records), expectedTags);
        Assert.assertEquals(tagging.duplicates(), expectedTags.keySet());
        final DuplicationMetrics metrics = tagging.onlyMetrics();
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, setSize);
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES, setSize - 1);
        Assert.assertEquals(metrics.READ_PAIR_OPTICAL_DUPLICATES, 1);

        final MarkDuplicatesEndToEndRunner.Result removing = input.run("REMOVE_SEQUENCING_DUPLICATES=true");
        final List<String> expectedKeys = input.inputRecords().stream()
                .map(MarkDuplicatesEndToEndRunner::key)
                .filter(key -> !key.startsWith(tagged + "/"))
                .collect(Collectors.toList());
        Assert.assertEquals(removing.keys(), expectedKeys);
        final List<String> expectedDuplicates = new ArrayList<>(expectedTags.keySet());
        expectedDuplicates.removeIf(key -> key.startsWith(tagged + "/"));
        removing.assertDuplicates(expectedDuplicates.toArray(new String[0]));
    }

    private static void addPair(final MarkDuplicatesEndToEndRunner input, final String name, final int quality) {
        input.pair(name, "rg", 0, 100, false, "10M", 0, 291, true, "10M", quality);
    }

    private static void putBothEnds(final Map<String, Object> tags, final String name, final String value) {
        tags.put(name + "/1", value);
        tags.put(name + "/2", value);
    }

    /** The DT tag of every output record that has one, by record key. */
    private static Map<String, Object> duplicateTypeTags(final List<SAMRecord> records) {
        final Map<String, Object> tags = new TreeMap<>();
        for (final SAMRecord rec : records) {
            final Object dt = rec.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG);
            if (dt != null) {
                tags.put(MarkDuplicatesEndToEndRunner.key(rec), dt);
            }
        }
        return tags;
    }
}
