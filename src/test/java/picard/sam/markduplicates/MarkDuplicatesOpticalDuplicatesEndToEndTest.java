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
import htsjdk.samtools.SAMTag;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.Result;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.TreeSet;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertMetrics;
import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertSameRecords;

/**
 * Tests optical duplicate detection end to end through {@link MarkDuplicates} on coordinate-sorted input:
 * READ_NAME_REGEX (default parsing of five and seven field names, a custom regex, null), OPTICAL_DUPLICATE_PIXEL_DISTANCE
 * and MAX_OPTICAL_DUPLICATE_SET_SIZE at and one past their limits, the requirement of the same read group and tile,
 * clusters linked through intermediate reads, the separation of pairs by which end was sequenced first, and the
 * READ_PAIR_OPTICAL_DUPLICATES, ESTIMATED_LIBRARY_SIZE, set-size histograms and DT:Z:SQ tags that result.
 *
 * Every duplicate set here is a set of pairs with 5' ends chr1:100 forward and chr1:300 reverse, the first pair
 * scoring highest; optical duplicates are only ever pairs, never fragments.
 */
public class MarkDuplicatesOpticalDuplicatesEndToEndTest {
    private static final int CHR1 = 0, CHR2 = 1;
    private static final boolean F = false, R = true;

    /** A read name with the five ':' separated fields the default READ_NAME_REGEX parses. */
    private static String illuminaName(final String prefix, final int tile, final int x, final int y) {
        return prefix + ":1:" + tile + ":" + x + ":" + y;
    }

    private static void addPair(final MarkDuplicatesEndToEndRunner input, final String name, final String readGroup, final int quality) {
        input.pair(name, readGroup, CHR1, 100, F, "10M", CHR1, 291, R, "10M", quality);
    }

    /** The keys of the records tagged DT:Z:SQ. */
    private static TreeSet<String> sequencingTagged(final Result result) {
        return result.records.stream()
                .filter(r -> MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(r.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG)))
                .map(MarkDuplicatesEndToEndRunner::key)
                .collect(Collectors.toCollection(TreeSet::new));
    }

    private static TreeSet<String> bothEnds(final String... names) {
        final TreeSet<String> keys = new TreeSet<>();
        for (final String name : names) {
            keys.add(name + "/1");
            keys.add(name + "/2");
        }
        return keys;
    }

    private static String[] concat(final String[] args, final String... more) {
        final List<String> all = new ArrayList<>(Arrays.asList(args));
        all.addAll(Arrays.asList(more));
        return all.toArray(new String[0]);
    }

    /**
     * Runs a single duplicate set of pairs in library "lib", the first one the keeper, and checks the duplicate
     * flags, the metrics, the set-size histograms and that exactly expectedOptical pairs carry DT:Z:SQ.
     */
    private static Result runOneSet(final MarkDuplicatesEndToEndRunner input, final List<String> names,
                                    final int expectedOptical, final boolean checkHistograms, final String... args) {
        final Result result = input.run(concat(args, "TAGGING_POLICY=OpticalOnly"));
        assertSameRecords(input.inputRecords(), result.records, SAMTag.PG.name(), MarkDuplicates.DUPLICATE_TYPE_TAG);
        result.assertDuplicates(bothEnds(names.subList(1, names.size()).toArray(new String[0])).toArray(new String[0]));
        final int n = names.size();
        assertMetrics(result.onlyMetrics(), "lib", 0, n, 0, 0, 0, n - 1, expectedOptical);
        Assert.assertEquals(sequencingTagged(result).size(), 2 * expectedOptical);
        Assert.assertFalse(sequencingTagged(result).contains(names.get(0) + "/1"));
        if (checkHistograms) {
            result.assertSetSizeHistograms(new double[]{n, 1},
                    expectedOptical == 0 ? new double[0] : new double[]{expectedOptical + 1, 1},
                    new double[]{n - expectedOptical, 1});
        }
        return result;
    }

    @DataProvider(name = "pixelDistances")
    public Object[][] pixelDistances() {
        return new Object[][]{
                // the default distance is 100
                {null, 100, 0, true},
                {null, 101, 0, false},
                {null, 0, 100, true},
                {null, 0, 101, false},
                {null, 100, 100, true},
                {null, -100, -100, true},
                {null, -101, 0, false},
                {null, 0, -101, false},
                {250, 250, 0, true},
                {250, 251, 0, false},
                {250, 0, 250, true},
                {250, 0, 251, false},
                {0, 0, 0, true},
                {0, 1, 0, false},
                {0, 0, 1, false},
        };
    }

    @Test(dataProvider = "pixelDistances")
    public void testOpticalDuplicatePixelDistance(final Integer distance, final int dx, final int dy, final boolean optical) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final String keeper = illuminaName("K", 1101, 5000, 5000);
        final String other = illuminaName("O", 1101, 5000 + dx, 5000 + dy);
        addPair(input, keeper, "rg", 30);
        addPair(input, other, "rg", 20);

        final Result result = distance == null ?
                runOneSet(input, Arrays.asList(keeper, other), optical ? 1 : 0, true) :
                runOneSet(input, Arrays.asList(keeper, other), optical ? 1 : 0, true, "OPTICAL_DUPLICATE_PIXEL_DISTANCE=" + distance);
        Assert.assertEquals(sequencingTagged(result), optical ? bothEnds(other) : new TreeSet<>());
    }

    @DataProvider(name = "readGroupsAndTiles")
    public Object[][] readGroupsAndTiles() {
        return new Object[][]{
                {new String[]{"close"}, 1},
                {new String[]{"otherTile"}, 0},
                {new String[]{"otherReadGroup"}, 0},
                {new String[]{"otherTile", "close"}, 1},
                {new String[]{"otherReadGroup", "close"}, 1},
                {new String[]{"otherReadGroup", "otherReadGroupClose"}, 1},
                {new String[]{"otherTile", "otherReadGroup", "close"}, 1},
                {new String[]{"otherTile", "otherReadGroup", "otherReadGroupClose", "close"}, 2},
        };
    }

    /** Optical duplicates must share the read group and the tile, even within one library. */
    @Test(dataProvider = "readGroupsAndTiles")
    public void testOpticalDuplicatesNeedTheSameReadGroupAndTile(final String[] members, final int expectedOptical) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner()
                .readGroup("rgA", "lib")
                .readGroup("rgB", "lib");
        final List<String> names = new ArrayList<>();
        names.add(illuminaName("keeper", 1101, 1000, 1000));
        addPair(input, names.get(0), "rgA", 30);
        for (final String member : members) {
            final String name;
            final String readGroup;
            switch (member) {
                case "close":
                    name = illuminaName(member, 1101, 1010, 1010);
                    readGroup = "rgA";
                    break;
                case "otherTile":
                    name = illuminaName(member, 1102, 1000, 1000);
                    readGroup = "rgA";
                    break;
                case "otherReadGroup":
                    name = illuminaName(member, 1101, 1000, 1000);
                    readGroup = "rgB";
                    break;
                case "otherReadGroupClose":
                    name = illuminaName(member, 1101, 1005, 1005);
                    readGroup = "rgB";
                    break;
                default:
                    throw new IllegalArgumentException(member);
            }
            names.add(name);
            addPair(input, name, readGroup, 20);
        }
        runOneSet(input, names, expectedOptical, true);
    }

    @DataProvider(name = "readNames")
    public Object[][] readNames() {
        final String semicolonRegex = "(?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*$";
        return new Object[][]{
                {null, "K:1:1101:1000:1000", "O:1:1101:1050:1050", 1},
                // seven fields: the last three are tile, x and y
                {null, "MA:55:FCA:1:1101:1000:1000", "MB:56:FCB:2:1101:1050:1050", 1},
                {null, "MA:55:FCA:1:1101:1000:1000", "MB:56:FCB:2:1102:1050:1050", 0},
                // only the leading digits of a field are parsed
                {null, "K:1:1101:1000:1000#0/1", "O:1:1101:1050:1050#0/1", 1},
                // six, four and three fields are not parsed
                {null, "a:b:K:1101:1000:1000", "a:b:O:1101:1050:1050", 0},
                {null, "K:1101:1000:1000", "O:1101:1050:1050", 0},
                {null, "1101:1000:1000", "1101:1050:1050", 0},
                // a field that does not start with a digit
                {null, "K:1:1101:x1000:1000", "O:1:1101:x1050:1050", 0},
                {semicolonRegex, "K;1101;1000;1000", "O;1101;1050;1050", 1},
                {semicolonRegex, "K;1101;1000;1000", "O;1101;1101;1000", 0},
                {semicolonRegex, "K;1101;1000;1000", "O;1102;1050;1050", 0},
                {semicolonRegex, "K:1:1101:1000:1000", "O:1:1101:1050:1050", 0},
        };
    }

    @Test(dataProvider = "readNames")
    public void testReadNameRegex(final String regex, final String keeper, final String other, final int expectedOptical) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        addPair(input, keeper, "rg", 30);
        addPair(input, other, "rg", 20);
        if (regex == null) {
            runOneSet(input, Arrays.asList(keeper, other), expectedOptical, true);
        } else {
            runOneSet(input, Arrays.asList(keeper, other), expectedOptical, true, "READ_NAME_REGEX=" + regex);
        }
    }

    @Test
    public void testNullReadNameRegexFindsNoOpticalDuplicates() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final List<String> names = Arrays.asList(illuminaName("K", 1101, 1000, 1000), illuminaName("O", 1101, 1000, 1000));
        addPair(input, names.get(0), "rg", 30);
        addPair(input, names.get(1), "rg", 20);
        runOneSet(input, names, 0, false, "READ_NAME_REGEX=null");
    }

    @DataProvider(name = "maxSetSizes")
    public Object[][] maxSetSizes() {
        return new Object[][]{
                {2, null, 1},
                {2, 2, 1},
                {2, 1, 0},
                {3, null, 2},
                {3, 3, 2},
                {3, 2, 0},
                {4, 5, 3},
                {4, 4, 3},
                {4, 3, 0},
        };
    }

    /** Sets larger than MAX_OPTICAL_DUPLICATE_SET_SIZE are not examined for optical duplicates. */
    @Test(dataProvider = "maxSetSizes")
    public void testMaxOpticalDuplicateSetSize(final int setSize, final Integer max, final int expectedOptical) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final List<String> names = new ArrayList<>();
        for (int i = 0; i < setSize; i++) {
            names.add(illuminaName("p" + i, 1101, 1000 + 10 * i, 1000));
            addPair(input, names.get(i), "rg", i == 0 ? 30 : 20);
        }
        if (max == null) {
            runOneSet(input, names, expectedOptical, true);
        } else {
            runOneSet(input, names, expectedOptical, true, "MAX_OPTICAL_DUPLICATE_SET_SIZE=" + max);
        }
    }

    @DataProvider(name = "clusters")
    public Object[][] clusters() {
        return new Object[][]{
                // x positions on one tile with y fixed; the first is the keeper
                {new int[]{1000, 1100, 1200}, 2},
                {new int[]{1000, 3000, 3050}, 1},
                {new int[]{1000, 3000, 3200}, 0},
                {new int[]{1000, 1100, 1200, 3000, 3050}, 3},
                {new int[]{1000, 3000, 3100, 3200, 3300}, 3},
                {new int[]{1000, 1050, 3000, 3050}, 2},
                {new int[]{1000, 1101, 1202, 1303}, 0},
        };
    }

    /**
     * Reads within the pixel distance of each other form a cluster, also through intermediate reads; each cluster
     * keeps one read that is not an optical duplicate, which is the keeper of the set when it is in the cluster.
     */
    @Test(dataProvider = "clusters")
    public void testClustersOfOpticalDuplicates(final int[] xs, final int expectedOptical) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final List<String> names = new ArrayList<>();
        for (int i = 0; i < xs.length; i++) {
            names.add(illuminaName("p" + i, 1101, xs[i], 1000));
            addPair(input, names.get(i), "rg", i == 0 ? 30 : 20);
        }
        runOneSet(input, names, expectedOptical, true);
    }

    @Test
    public void testClusterLinkedThroughTheKeeper() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final List<String> names = new ArrayList<>();
        final int[] xs = {1000, 1100, 1200, 3000, 3050};
        for (int i = 0; i < xs.length; i++) {
            names.add(illuminaName("p" + i, 1101, xs[i], 1000));
            addPair(input, names.get(i), "rg", i == 0 ? 30 : 20);
        }
        final Result result = runOneSet(input, names, 3, true);
        final TreeSet<String> tagged = sequencingTagged(result);
        Assert.assertTrue(tagged.containsAll(bothEnds(names.get(1), names.get(2))));
        // one of the two reads of the second cluster stays untagged
        Assert.assertEquals(tagged.stream().filter(k -> k.startsWith("p3:") || k.startsWith("p4:")).count(), 2);
    }

    @DataProvider(name = "firstOfPairStrands")
    public Object[][] firstOfPairStrands() {
        return new Object[][]{
                {new String[]{"sameStrand"}, 1},
                {new String[]{"otherStrand"}, 0},
                {new String[]{"otherStrand", "sameStrand"}, 1},
                {new String[]{"otherStrand", "otherStrand2"}, 1},
                {new String[]{"otherStrand", "otherStrand2", "sameStrand"}, 2},
        };
    }

    /**
     * Pairs whose first of pair is on different strands are duplicates of each other but never optical duplicates
     * of each other.
     */
    @Test(dataProvider = "firstOfPairStrands")
    public void testOpticalDuplicatesNeedTheSameFirstOfPairStrand(final String[] members, final int expectedOptical) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        final List<String> names = new ArrayList<>();
        names.add(illuminaName("keeper", 1101, 1000, 1000));
        addPair(input, names.get(0), "rg", 30);
        int offset = 0;
        for (final String member : members) {
            offset += 10;
            final String name = illuminaName(member, 1101, 1000 + offset, 1000);
            names.add(name);
            if (member.startsWith("sameStrand")) {
                addPair(input, name, "rg", 20);
            } else {
                // the first of pair is the reverse read
                input.pair(name, "rg", CHR1, 291, R, "10M", CHR1, 100, F, "10M", 20);
            }
        }
        runOneSet(input, names, expectedOptical, true);
    }

    @Test
    public void testFragmentsAreNeverOpticalDuplicates() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        input.fragment(illuminaName("f1", 1101, 1000, 1000), "rg", CHR2, 100, F, "10M", 30);
        input.fragment(illuminaName("f2", 1101, 1000, 1000), "rg", CHR2, 100, F, "10M", 20);
        final Result result = input.run("TAGGING_POLICY=All");
        result.assertDuplicates(illuminaName("f2", 1101, 1000, 1000));
        assertMetrics(result.onlyMetrics(), "lib", 2, 0, 0, 0, 1, 0, 0);
        final SAMRecord dup = result.record(illuminaName("f2", 1101, 1000, 1000));
        Assert.assertEquals(dup.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG), MarkDuplicates.DUPLICATE_TYPE_LIBRARY);
    }

    @Test
    public void testOpticalDuplicatesPerLibrary() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner()
                .readGroup("rgA", "libA")
                .readGroup("rgB", "libB");
        // libA: 3 pairs, 1 optical; libB: 3 pairs, 2 optical, 1 more pair that is not a duplicate
        addPair(input, illuminaName("a0", 1101, 1000, 1000), "rgA", 30);
        addPair(input, illuminaName("a1", 1101, 1050, 1000), "rgA", 20);
        addPair(input, illuminaName("a2", 1101, 9000, 1000), "rgA", 20);
        addPair(input, illuminaName("b0", 1101, 1000, 1000), "rgB", 30);
        addPair(input, illuminaName("b1", 1101, 1050, 1000), "rgB", 20);
        addPair(input, illuminaName("b2", 1101, 1100, 1000), "rgB", 20);
        input.pair(illuminaName("b3", 1101, 1000, 1000), "rgB", CHR2, 100, F, "10M", CHR2, 291, R, "10M", 20);

        final Result result = input.run("TAGGING_POLICY=OpticalOnly");
        assertMetrics(result.metrics("libA"), "libA", 0, 3, 0, 0, 0, 2, 1);
        assertMetrics(result.metrics("libB"), "libB", 0, 4, 0, 0, 0, 2, 2);
        Assert.assertEquals(sequencingTagged(result), bothEnds(illuminaName("a1", 1101, 1050, 1000),
                illuminaName("b1", 1101, 1050, 1000), illuminaName("b2", 1101, 1100, 1000)));
        // the histograms pool the libraries
        result.assertSetSizeHistograms(new double[]{1, 1, 3, 2}, new double[]{2, 1, 3, 1}, new double[]{1, 2, 2, 1});
        Assert.assertNull(result.histogram("CoverageMult"));
    }
}
