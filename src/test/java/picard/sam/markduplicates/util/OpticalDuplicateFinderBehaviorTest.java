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

package picard.sam.markduplicates.util;

import htsjdk.samtools.util.Log;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.util.PhysicalLocation;
import picard.sam.util.ReadNameParser;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Tests {@link OpticalDuplicateFinder#findOpticalDuplicates} on hand-built location lists: the pixel distance at and
 * one past the threshold in x and y, tiles and read groups, reads without a location, the keeper argument (in the
 * list, absent, without a location, or only equal by value), the comparison used for two or three reads and the
 * clustering used for larger sets, the maximum set size, the progress logging for large sets, and the read name
 * parsing the finder inherits.
 */
public class OpticalDuplicateFinderBehaviorTest {

    private static final int TILE = 1101;

    /** Builds a location in the given read group and tile. */
    private static ReadEndsForMarkDuplicates loc(final int readGroup, final int tile, final int x, final int y) {
        final ReadEndsForMarkDuplicates loc = new ReadEndsForMarkDuplicates();
        loc.setReadGroup((short) readGroup);
        loc.setTile((short) tile);
        loc.setX(x);
        loc.setY(y);
        return loc;
    }

    private static ReadEndsForMarkDuplicates loc(final int x, final int y) {
        return loc(1, TILE, x, y);
    }

    /** A read whose name could not be parsed: tile stays at -1. */
    private static ReadEndsForMarkDuplicates noLoc(final int x, final int y) {
        final ReadEndsForMarkDuplicates loc = new ReadEndsForMarkDuplicates();
        loc.setReadGroup((short) 1);
        loc.setX(x);
        loc.setY(y);
        Assert.assertFalse(loc.hasLocation());
        return loc;
    }

    private static OpticalDuplicateFinder finder(final int distance) {
        return new OpticalDuplicateFinder(ReadNameParser.DEFAULT_READ_NAME_REGEX, distance, null);
    }

    private static void assertFlags(final boolean[] actual, final boolean... expected) {
        Assert.assertEquals(actual.length, expected.length, "flag array length");
        if (!Arrays.equals(actual, expected)) {
            Assert.fail("expected " + Arrays.toString(expected) + " but was " + Arrays.toString(actual));
        }
    }

    private static int countTrue(final boolean[] flags) {
        int count = 0;
        for (final boolean flag : flags) {
            if (flag) count++;
        }
        return count;
    }

    /** Runs the action with the global log stream redirected and returns what was written. */
    private static String captureLog(final Runnable action) {
        final PrintStream originalStream = Log.getGlobalPrintStream();
        final Log.LogLevel originalLevel = Log.getGlobalLogLevel();
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream capture = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            Log.setGlobalPrintStream(capture);
            Log.setGlobalLogLevel(Log.LogLevel.INFO);
            action.run();
        } finally {
            Log.setGlobalPrintStream(originalStream);
            Log.setGlobalLogLevel(originalLevel);
        }
        return bytes.toString(StandardCharsets.UTF_8);
    }

    @DataProvider(name = "pixelDistanceData")
    public Object[][] pixelDistanceData() {
        return new Object[][]{
                // distance, dx, dy, within distance
                {100, 0, 0, true},
                {100, 100, 0, true},
                {100, 101, 0, false},
                {100, -100, 0, true},
                {100, -101, 0, false},
                {100, 0, 100, true},
                {100, 0, 101, false},
                {100, 0, -100, true},
                {100, 0, -101, false},
                {100, 100, 100, true},
                {100, 100, -100, true},
                {100, -100, 100, true},
                {100, 101, 100, false},
                {100, 100, 101, false},
                {100, -101, -101, false},
                {100, 1, 1000, false},
                {0, 0, 0, true},
                {0, 1, 0, false},
                {0, 0, 1, false},
                {0, -1, -1, false},
                {2500, 2500, -2500, true},
                {2500, 2501, 0, false},
                {2500, 0, -2501, false},
        };
    }

    /** Two reads: the comparison used when there is no keeper, with the second read flagged. */
    @Test(dataProvider = "pixelDistanceData")
    public void testPixelDistanceTwoReads(final int distance, final int dx, final int dy, final boolean close) {
        final PhysicalLocation base = loc(5000, 5000);
        final PhysicalLocation other = loc(5000 + dx, 5000 + dy);
        assertFlags(finder(distance).findOpticalDuplicates(Arrays.asList(base, other), null), false, close);
    }

    /** Two reads with one of them the keeper: only the other one can be flagged. */
    @Test(dataProvider = "pixelDistanceData")
    public void testPixelDistanceTwoReadsWithKeeper(final int distance, final int dx, final int dy, final boolean close) {
        final PhysicalLocation base = loc(5000, 5000);
        final PhysicalLocation other = loc(5000 + dx, 5000 + dy);
        final List<PhysicalLocation> list = Arrays.asList(base, other);
        assertFlags(finder(distance).findOpticalDuplicates(list, base), false, close);
        assertFlags(finder(distance).findOpticalDuplicates(list, other), close, false);
    }

    /** Four reads, two of them in other tiles at the same coordinates: the clustering path keeps the smaller (x, y). */
    @Test(dataProvider = "pixelDistanceData")
    public void testPixelDistanceClustered(final int distance, final int dx, final int dy, final boolean close) {
        final PhysicalLocation base = loc(5000, 5000);
        final PhysicalLocation other = loc(5000 + dx, 5000 + dy);
        final List<PhysicalLocation> list = Arrays.asList(base, other, loc(1, TILE + 1, 5000, 5000), loc(1, TILE + 2, 5000, 5000));
        final boolean[] flags = finder(distance).findOpticalDuplicates(list, null);

        if (!close) {
            assertFlags(flags, false, false, false, false);
        } else if (dx == 0 && dy == 0) {
            // which of two identical locations is kept is not specified
            Assert.assertEquals(countTrue(flags), 1);
            Assert.assertFalse(flags[2]);
            Assert.assertFalse(flags[3]);
        } else {
            final boolean otherIsLarger = dx > 0 || (dx == 0 && dy > 0);
            assertFlags(flags, !otherIsLarger, otherIsLarger, false, false);
        }
    }

    @Test
    public void testNegativeCoordinates() {
        final OpticalDuplicateFinder finder = finder(100);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(-50, -60), loc(50, 40)), null), false, true);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(-50, -60), loc(51, 40)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(50, 40), loc(-50, -60), loc(1, TILE + 1, 0, 0)), null),
                true, false, false);
    }

    @Test
    public void testTilesAndReadGroupsSeparateLocations() {
        final OpticalDuplicateFinder finder = finder(100);

        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(1, 1101, 100, 100), loc(1, 1102, 100, 100)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(1, 1101, 100, 100), loc(2, 1101, 100, 100)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(1, 1101, 100, 100), loc(1, 1101, 100, 100)), null), false, true);

        // Distinct (read group, tile) pairs, chosen so that packing both into one small key would collide.
        final List<PhysicalLocation> list = new ArrayList<>(Arrays.asList(
                loc(0, 1, 100, 100),
                loc(1, 1, 100, 100),
                loc(0, 2, 100, 100),
                loc(2, 1, 100, 100),
                loc(0, 256, 100, 100),
                loc(1, 0, 100, 100),
                loc(0, 0, 100, 100),
                loc(256, 0, 100, 100),
                loc(0, 32767, 100, 100),
                loc(32767, 1, 100, 100),
                loc(1, 32767, 100, 100),
                loc(0, 10001, 100, 100),
                loc(1, 1, 100, 100 + 101)
        ));
        assertFlags(finder.findOpticalDuplicates(list, null), new boolean[list.size()]);

        list.add(loc(0, 256, 150, 50));
        final boolean[] expected = new boolean[list.size()];
        expected[list.size() - 1] = true;
        assertFlags(finder.findOpticalDuplicates(list, null), expected);
    }

    @Test
    public void testReadsWithoutLocationAreNeverFlagged() {
        final OpticalDuplicateFinder finder = finder(100);

        assertFlags(finder.findOpticalDuplicates(Arrays.asList(noLoc(100, 100), noLoc(100, 100)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(noLoc(100, 100), loc(100, 100)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(100, 100), noLoc(100, 100)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(
                Arrays.asList(noLoc(100, 100), noLoc(100, 100), loc(100, 100), noLoc(100, 100)), null),
                false, false, false, false);

        final PhysicalLocation keeper = loc(100, 100);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(noLoc(100, 100), keeper), keeper), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(noLoc(100, 100), keeper, noLoc(100, 100)), keeper),
                false, false, false);
    }

    /** With two reads and no keeper the later read of a close pair is flagged, whatever its coordinates. */
    @Test
    public void testTwoReadsWithoutKeeperFlagSecond() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation larger = loc(200, 0);
        final PhysicalLocation smaller = loc(100, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(larger, smaller), null), false, true);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(smaller, larger), null), false, true);
    }

    /**
     * Three reads take the clustering path when the keeper is null, but the pairwise path when a keeper is given
     * that is absent from the list or has no location.
     */
    @Test
    public void testThreeReadsPathDependsOnKeeperArgument() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation a = loc(200, 0);
        final PhysicalLocation b = loc(100, 0);
        final PhysicalLocation c = loc(1, TILE + 1, 0, 0);
        final List<PhysicalLocation> list = Arrays.asList(a, b, c);

        assertFlags(finder.findOpticalDuplicates(list, null), true, false, false);
        assertFlags(finder.findOpticalDuplicates(list, loc(150, 0)), false, true, false);

        final PhysicalLocation unlocatedKeeper = noLoc(150, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(a, b, unlocatedKeeper), unlocatedKeeper), false, true, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(a, b, unlocatedKeeper, c), unlocatedKeeper), true, false, false, false);
    }

    @Test
    public void testKeeperWithThreeReads() {
        final OpticalDuplicateFinder finder = finder(100);

        PhysicalLocation k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(k, loc(100, 0), loc(200, 0)), k), false, true, true);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(k, loc(300, 0), loc(200, 0)), k), false, false, true);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(200, 0), loc(300, 0), k), k), false, true, false);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(k, loc(500, 0), loc(100, 0)), k), false, false, true);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(100, 0), k, loc(200, 0)), k), true, false, true);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(200, 0), loc(100, 0), k), k), true, true, false);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(k, loc(50, 50), loc(-50, -50)), k), false, true, true);

        k = loc(0, 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(k, loc(1, TILE + 1, 0, 0), loc(2, TILE, 0, 0)), k), false, false, false);
    }

    @Test
    public void testKeeperWithClustering() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation a = loc(300, 0);
        final PhysicalLocation k = loc(400, 0);
        final PhysicalLocation b = loc(250, 0);
        final PhysicalLocation c = loc(2000, 0);
        final PhysicalLocation d = loc(2050, 0);
        final List<PhysicalLocation> list = Arrays.asList(a, k, b, c, d);

        // the keeper's cluster flags every other member, even ones with smaller coordinates
        assertFlags(finder.findOpticalDuplicates(list, k), true, false, true, false, true);
        assertFlags(finder.findOpticalDuplicates(list, null), true, true, false, false, true);
        // equal by value but not the same object, so it is not in the list
        assertFlags(finder.findOpticalDuplicates(list, loc(400, 0)), true, true, false, false, true);
        assertFlags(finder.findOpticalDuplicates(list, c), true, true, false, false, true);
        assertFlags(finder.findOpticalDuplicates(list, d), true, true, false, true, false);
    }

    @Test
    public void testTransitiveChainIsOneCluster() {
        final OpticalDuplicateFinder finder = finder(10);
        final PhysicalLocation a = loc(0, 0);
        final PhysicalLocation b = loc(30, 0);
        final PhysicalLocation c = loc(10, 0);
        final PhysicalLocation d = loc(20, 0);
        final List<PhysicalLocation> list = Arrays.asList(a, b, c, d);

        assertFlags(finder.findOpticalDuplicates(list, null), false, true, true, true);
        assertFlags(finder.findOpticalDuplicates(list, b), true, false, true, true);
        assertFlags(finder.findOpticalDuplicates(list, d), true, true, true, false);

        final OpticalDuplicateFinder narrower = finder(9);
        assertFlags(narrower.findOpticalDuplicates(list, null), false, false, false, false);
    }

    @Test
    public void testMinimumXThenMinimumYIsKept() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation pad = loc(1, TILE + 1, 0, 0);

        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(100, 200), loc(100, 150), loc(100, 250), pad), null),
                true, false, true, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(150, 0), loc(100, 90), pad, loc(1, TILE + 2, 0, 0)), null),
                true, false, false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(5, -40), loc(5, -41), loc(6, -100), pad), null),
                true, false, true, false);
    }

    @Test
    public void testIdenticalLocationsKeepExactlyOne() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation pad = loc(1, TILE + 1, 100, 100);

        final boolean[] two = finder.findOpticalDuplicates(Arrays.asList(loc(100, 100), loc(100, 100), pad, loc(1, TILE + 2, 100, 100)), null);
        Assert.assertEquals(countTrue(two), 1);
        Assert.assertFalse(two[2]);
        Assert.assertFalse(two[3]);

        final boolean[] three = finder.findOpticalDuplicates(Arrays.asList(loc(100, 100), loc(100, 100), loc(100, 100)), null);
        Assert.assertEquals(countTrue(three), 2);

        // a smaller location in the cluster is kept over the identical ones
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(100, 100), loc(100, 100), loc(99, 150), pad), null),
                true, true, false, false);
    }

    @Test
    public void testTwoClustersAndUnlocatedReads() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation k = loc(1000, 1000);
        final List<PhysicalLocation> list = Arrays.asList(
                noLoc(0, 0), loc(0, 0), k, loc(50, 50), noLoc(50, 50), loc(1050, 950));

        assertFlags(finder.findOpticalDuplicates(list, k), false, false, false, true, false, true);
        assertFlags(finder.findOpticalDuplicates(list, null), false, false, false, true, false, true);
        assertFlags(finder.findOpticalDuplicates(list, list.get(5)), false, false, true, true, false, false);
        assertFlags(finder.findOpticalDuplicates(list, list.get(0)), false, false, false, true, false, true);
    }

    @Test
    public void testEmptyAndSingleElementLists() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation only = loc(100, 100);

        assertFlags(finder.findOpticalDuplicates(Collections.emptyList(), null));
        assertFlags(finder.findOpticalDuplicates(Collections.emptyList(), only));
        assertFlags(finder.findOpticalDuplicates(Collections.singletonList(only), null), false);
        assertFlags(finder.findOpticalDuplicates(Collections.singletonList(only), only), false);
        assertFlags(finder.findOpticalDuplicates(Collections.singletonList(only), loc(100, 100)), false);
    }

    @Test
    public void testInputListIsNotModified() {
        final OpticalDuplicateFinder finder = finder(100);
        final PhysicalLocation a = loc(300, 0);
        final PhysicalLocation k = loc(400, 0);
        final PhysicalLocation b = loc(250, 0);
        final List<PhysicalLocation> list = new ArrayList<>(Arrays.asList(a, k, b, loc(1, TILE + 1, 0, 0)));
        final List<PhysicalLocation> copy = new ArrayList<>(list);

        finder.findOpticalDuplicates(list, k);
        Assert.assertEquals(list.size(), copy.size());
        for (int i = 0; i < list.size(); i++) {
            Assert.assertSame(list.get(i), copy.get(i));
        }
        Assert.assertEquals(a.getX(), 300);
        Assert.assertEquals(b.getX(), 250);
        Assert.assertEquals(k.getTile(), TILE);
    }

    @DataProvider(name = "maxSetSizeData")
    public Object[][] maxSetSizeData() {
        return new Object[][]{
                // max set size, list size, expected flags for a list of mutually close reads with no keeper
                {1L, 1, new boolean[]{false}},
                {1L, 2, new boolean[]{false, false}},
                {2L, 2, new boolean[]{false, true}},
                {2L, 3, new boolean[]{false, false, false}},
                {3L, 3, new boolean[]{false, true, true}},
                {3L, 4, new boolean[]{false, false, false, false}},
                {4L, 4, new boolean[]{false, true, true, true}},
                {4L, 5, new boolean[]{false, false, false, false, false}},
                {5L, 5, new boolean[]{false, true, true, true, true}},
                {Long.MAX_VALUE, 5, new boolean[]{false, true, true, true, true}},
        };
    }

    @Test(dataProvider = "maxSetSizeData")
    public void testMaxDuplicateSetSize(final long maxSetSize, final int listSize, final boolean[] expected) {
        final List<PhysicalLocation> list = new ArrayList<>();
        for (int i = 0; i < listSize; i++) {
            list.add(loc(10 * i, 0));
        }

        final OpticalDuplicateFinder fromConstructor =
                new OpticalDuplicateFinder(ReadNameParser.DEFAULT_READ_NAME_REGEX, 100, maxSetSize, null);
        assertFlags(fromConstructor.findOpticalDuplicates(list, null), expected);

        final OpticalDuplicateFinder fromSetter = new OpticalDuplicateFinder();
        fromSetter.setMaxDuplicateSetSize(maxSetSize);
        assertFlags(fromSetter.findOpticalDuplicates(list, null), expected);
    }

    @Test
    public void testMaxDuplicateSetSizeWithKeeper() {
        final PhysicalLocation k = loc(0, 0);
        final List<PhysicalLocation> list = Arrays.asList(loc(30, 0), k, loc(10, 0), loc(20, 0));

        final OpticalDuplicateFinder atMax = new OpticalDuplicateFinder(ReadNameParser.DEFAULT_READ_NAME_REGEX, 100, 4, null);
        assertFlags(atMax.findOpticalDuplicates(list, k), true, false, true, true);

        final OpticalDuplicateFinder belowMax = new OpticalDuplicateFinder(ReadNameParser.DEFAULT_READ_NAME_REGEX, 100, 3, null);
        assertFlags(belowMax.findOpticalDuplicates(list, k), false, false, false, false);
    }

    @Test(timeOut = 120_000)
    public void testDefaultMaxDuplicateSetSize() {
        Assert.assertEquals(OpticalDuplicateFinder.DEFAULT_MAX_DUPLICATE_SET_SIZE, 300000);
        final List<PhysicalLocation> tooMany = Collections.nCopies(OpticalDuplicateFinder.DEFAULT_MAX_DUPLICATE_SET_SIZE + 1, loc(100, 100));

        final boolean[] fromDefault = new OpticalDuplicateFinder().findOpticalDuplicates(tooMany, null);
        Assert.assertEquals(fromDefault.length, tooMany.size());
        Assert.assertEquals(countTrue(fromDefault), 0);

        final boolean[] fromThreeArgs = finder(100).findOpticalDuplicates(tooMany, null);
        Assert.assertEquals(fromThreeArgs.length, tooMany.size());
        Assert.assertEquals(countTrue(fromThreeArgs), 0);
    }

    @Test
    public void testDefaultConstructorDistance() {
        Assert.assertEquals(OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE, 100);
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder();
        Assert.assertEquals(finder.opticalDuplicatePixelDistance, 100);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(0, 0), loc(100, 100)), null), false, true);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(0, 0), loc(101, 0)), null), false, false);
    }

    @Test
    public void testPixelDistanceFieldIsReadAtCallTime() {
        final OpticalDuplicateFinder finder = finder(100);
        finder.opticalDuplicatePixelDistance = 10;
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(0, 0), loc(11, 0)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(0, 0), loc(10, 0)), null), false, true);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(0, 0), loc(11, 0), loc(10, 10)), null), false, true, true);
    }

    @Test
    public void testNullRegexFindsNoOpticalDuplicates() {
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(null, 100, null);
        final PhysicalLocation k = loc(100, 100);

        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(100, 100), loc(100, 100)), null), false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(100, 100), loc(100, 100), loc(100, 100)), null), false, false, false);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(k, loc(100, 100), loc(100, 100), loc(100, 100)), k), false, false, false, false);
        Assert.assertFalse(finder.addLocationInformation("RUN:1:1101:100:100", new ReadEndsForMarkDuplicates()));
    }

    @Test
    public void testDefaultLargeSetLogging() {
        Assert.assertEquals(OpticalDuplicateFinder.DEFAULT_BIG_DUPLICATE_SET_SIZE, 1000);
        final OpticalDuplicateFinder finder = finder(100);

        final List<PhysicalLocation> atThreshold = new ArrayList<>();
        for (int i = 0; i < 1000; i++) {
            atThreshold.add(noLoc(100, 100));
        }
        final String quiet = captureLog(() -> Assert.assertEquals(countTrue(finder.findOpticalDuplicates(atThreshold, null)), 0));
        Assert.assertFalse(quiet.contains("Large duplicate set"), quiet);

        final List<PhysicalLocation> pastThreshold = new ArrayList<>(atThreshold);
        pastThreshold.add(noLoc(100, 100));
        final String logged = captureLog(() -> Assert.assertEquals(countTrue(finder.findOpticalDuplicates(pastThreshold, null)), 0));
        Assert.assertTrue(logged.contains("Large duplicate set. size = 1001"), logged);
    }

    /** Lowering the large-set size turns on progress logging; the flags must not change. */
    @Test
    public void testLargeSetLoggingKeepsResults() {
        final OpticalDuplicateFinder finder = finder(100);

        finder.setBigDuplicateSetSize(2);
        final String twoReads = captureLog(() ->
                assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(200, 0), loc(100, 0)), null), false, true));
        Assert.assertFalse(twoReads.contains("Large duplicate set"), twoReads);

        final PhysicalLocation k3 = loc(0, 0);
        final String threeReads = captureLog(() ->
                assertFlags(finder.findOpticalDuplicates(Arrays.asList(k3, loc(300, 0), loc(200, 0)), k3), false, false, true));
        Assert.assertTrue(threeReads.contains("Large duplicate set. size = 3"), threeReads);

        final PhysicalLocation k5 = loc(400, 0);
        final List<PhysicalLocation> five = Arrays.asList(loc(300, 0), k5, loc(250, 0), loc(2000, 0), loc(2050, 0));
        final String fiveReads = captureLog(() -> {
            assertFlags(finder.findOpticalDuplicates(five, k5), true, false, true, false, true);
            assertFlags(finder.findOpticalDuplicates(five, null), true, true, false, false, true);
        });
        Assert.assertTrue(fiveReads.contains("Large duplicate set. size = 5"), fiveReads);

        finder.setBigDuplicateSetSize(1);
        final String twoReadsLogged = captureLog(() ->
                assertFlags(finder.findOpticalDuplicates(Arrays.asList(loc(200, 0), loc(100, 0)), null), false, true));
        Assert.assertTrue(twoReadsLogged.contains("Large duplicate set. size = 2"), twoReadsLogged);
    }

    /** The lane field of a seven-field name is not parsed; only the read group can separate lanes. */
    @Test
    public void testParsedNamesFromDifferentLanes() {
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder();
        final ReadEndsForMarkDuplicates lane1 = new ReadEndsForMarkDuplicates();
        final ReadEndsForMarkDuplicates lane2 = new ReadEndsForMarkDuplicates();
        Assert.assertTrue(finder.addLocationInformation("M01234:7:FC1:1:1101:1000:2000", lane1));
        Assert.assertTrue(finder.addLocationInformation("M01234:7:FC1:2:1101:1000:2000", lane2));
        Assert.assertEquals(lane1.getTile(), 1101);
        Assert.assertEquals(lane2.getTile(), 1101);
        Assert.assertEquals(lane1.getX(), 1000);
        Assert.assertEquals(lane2.getX(), 1000);
        Assert.assertEquals(lane1.getY(), 2000);
        Assert.assertEquals(lane2.getY(), 2000);

        lane1.setReadGroup((short) 0);
        lane2.setReadGroup((short) 0);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(lane1, lane2), null), false, true);

        lane2.setReadGroup((short) 1);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(lane1, lane2), null), false, false);
    }

    @Test
    public void testParsedNamesAtThreshold() {
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(ReadNameParser.DEFAULT_READ_NAME_REGEX, 2500, null);
        final String[] names = {
                "A00123:8:HABCDEFXX:1:2101:10000:20000",
                "A00123:8:HABCDEFXX:1:2101:12500:17500",
                "A00123:8:HABCDEFXX:1:2101:12501:20000",
                "A00123:8:HABCDEFXX:1:2102:10000:20000",
        };
        final List<PhysicalLocation> list = new ArrayList<>();
        for (final String name : names) {
            final ReadEndsForMarkDuplicates loc = new ReadEndsForMarkDuplicates();
            Assert.assertTrue(finder.addLocationInformation(name, loc), name);
            loc.setReadGroup((short) 0);
            list.add(loc);
        }
        // 0-1 at exactly 2500 in x and y, 1-2 one apart in x and 2500 in y, 0-2 one past in x, 3 in another tile
        assertFlags(finder.findOpticalDuplicates(list, null), false, true, true, false);
        assertFlags(finder.findOpticalDuplicates(list, list.get(2)), true, true, false, false);
        assertFlags(finder.findOpticalDuplicates(list.subList(0, 2), null), false, true);
        assertFlags(finder.findOpticalDuplicates(Arrays.asList(list.get(0), list.get(2)), null), false, false);
    }

    /**
     * Compares the clustering path against an independent reference on random lists: reads are linked when they share
     * read group and tile, both have a location and lie within the distance in x and in y; each linked group keeps the
     * keeper if it holds it, otherwise the read with the smallest x and then y. Coordinates are distinct so the kept
     * read is unique.
     */
    @Test
    public void testClusteringMatchesReference() {
        final int distance = 40;
        final OpticalDuplicateFinder finder = finder(distance);
        final Random random = new Random(42);

        for (int trial = 0; trial < 300; trial++) {
            final int n = 3 + random.nextInt(38);
            final List<PhysicalLocation> list = new ArrayList<>();
            final Set<Long> used = new HashSet<>();
            while (list.size() < n) {
                final int x = random.nextInt(400) - 100;
                final int y = random.nextInt(400) - 100;
                if (!used.add(((long) x << 32) ^ (y & 0xffffffffL))) continue;
                final int kind = random.nextInt(10);
                if (kind == 0) {
                    list.add(noLoc(x, y));
                } else {
                    list.add(loc(random.nextInt(2), TILE + random.nextInt(2), x, y));
                }
            }

            final int keeperChoice = random.nextInt(n + 2);
            final PhysicalLocation keeper;
            if (keeperChoice < n) {
                keeper = list.get(keeperChoice);
            } else if (keeperChoice == n) {
                keeper = loc(0, 0);
            } else {
                keeper = null;
            }
            // three reads with a non-null keeper use the pairwise path, which this reference does not describe
            if (keeper != null && n < 4) continue;

            final boolean[] expected = referenceFlags(list, keeper, distance);
            final boolean[] actual = finder.findOpticalDuplicates(list, keeper);
            if (!Arrays.equals(actual, expected)) {
                Assert.fail("trial " + trial + ": expected " + Arrays.toString(expected) + " but was " + Arrays.toString(actual));
            }
        }
    }

    private static boolean linked(final PhysicalLocation a, final PhysicalLocation b, final int distance) {
        return a.hasLocation() && b.hasLocation()
                && a.getReadGroup() == b.getReadGroup()
                && a.getTile() == b.getTile()
                && Math.abs((long) a.getX() - b.getX()) <= distance
                && Math.abs((long) a.getY() - b.getY()) <= distance;
    }

    private static boolean[] referenceFlags(final List<PhysicalLocation> list, final PhysicalLocation keeper, final int distance) {
        final int n = list.size();
        final boolean[] flags = new boolean[n];
        final boolean[] seen = new boolean[n];
        for (int start = 0; start < n; start++) {
            if (seen[start]) continue;
            final List<Integer> component = new ArrayList<>();
            final Deque<Integer> queue = new ArrayDeque<>();
            queue.add(start);
            seen[start] = true;
            while (!queue.isEmpty()) {
                final int i = queue.poll();
                component.add(i);
                for (int j = 0; j < n; j++) {
                    if (!seen[j] && linked(list.get(i), list.get(j), distance)) {
                        seen[j] = true;
                        queue.add(j);
                    }
                }
            }
            int kept = -1;
            for (final int i : component) {
                if (list.get(i) == keeper && keeper.hasLocation()) kept = i;
            }
            if (kept < 0) {
                for (final int i : component) {
                    if (kept < 0) {
                        kept = i;
                    } else {
                        final PhysicalLocation candidate = list.get(i);
                        final PhysicalLocation best = list.get(kept);
                        if (candidate.getX() < best.getX() || (candidate.getX() == best.getX() && candidate.getY() < best.getY())) {
                            kept = i;
                        }
                    }
                }
            }
            for (final int i : component) {
                flags[i] = i != kept;
            }
        }
        return flags;
    }
}
