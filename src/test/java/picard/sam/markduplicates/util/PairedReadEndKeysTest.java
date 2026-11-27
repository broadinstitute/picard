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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * Tests for {@link PairedReadEndKeys}: which ends {@code add} takes, {@code contains} against a set of the four key
 * fields (library id, read1 reference index, read1 coordinate, orientation) in any query order, the overflow rule and
 * the ends it hands over, and {@link PairedReadEndKeys#splitRuns()}, after which every query answers as before, on
 * key sets whose sorted runs are shorter than, as long as and longer than the 32,768 keys of a split array.
 */
public class PairedReadEndKeysTest {

    /** The four fields a key holds, ordered field by field. */
    private static final class Key implements Comparable<Key> {
        final int library;
        final int reference;
        final int coordinate;
        final int orientation;

        Key(final int library, final int reference, final int coordinate, final int orientation) {
            this.library = library;
            this.reference = reference;
            this.coordinate = coordinate;
            this.orientation = orientation;
        }

        static Key of(final ReadEndsForMarkDuplicates end) {
            return new Key(end.libraryId, end.read1ReferenceIndex, end.read1Coordinate, end.orientation);
        }

        ReadEndsForMarkDuplicates toEnd() {
            final ReadEndsForMarkDuplicates end = new ReadEndsForMarkDuplicates();
            end.libraryId = (short) library;
            end.read1ReferenceIndex = reference;
            end.read1Coordinate = coordinate;
            end.orientation = (byte) orientation;
            end.read2ReferenceIndex = reference;
            end.read2Coordinate = coordinate + 150;
            return end;
        }

        @Override
        public int compareTo(final Key other) {
            int c = Integer.compare(library, other.library);
            if (c == 0) {
                c = Integer.compare(reference, other.reference);
            }
            if (c == 0) {
                c = Integer.compare(coordinate, other.coordinate);
            }
            return c != 0 ? c : Integer.compare(orientation, other.orientation);
        }

        @Override
        public boolean equals(final Object o) {
            return o instanceof Key && compareTo((Key) o) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(library, reference, coordinate, orientation);
        }

        @Override
        public String toString() {
            return library + ":" + reference + ":" + coordinate + ":" + orientation;
        }
    }

    private static final int REFERENCES = 25;

    /**
     * Keys in coordinate order, as a coordinate-sorted input adds them: each a few positions past the one before, with
     * some repeated at once and some repeated later, over the given libraries.
     */
    private static List<Key> sortedKeys(final Random random, final int count, final int libraries) {
        final List<Key> keys = new ArrayList<>(count);
        int reference = 0;
        int coordinate = 1;
        while (keys.size() < count) {
            coordinate += random.nextInt(4);
            if (random.nextInt(100_000) == 0 && reference + 1 < REFERENCES) {
                reference++;
                coordinate = 1;
            }
            final Key key = new Key(random.nextInt(libraries), reference, coordinate, random.nextInt(2));
            keys.add(key);
            if (random.nextInt(5) == 0) {
                keys.add(key);
            }
            if (random.nextInt(7) == 0 && keys.size() > 10) {
                keys.add(keys.get(keys.size() - 1 - random.nextInt(10)));
            }
        }
        return keys;
    }

    private static List<Key> randomKeys(final Random random, final int count) {
        final List<Key> keys = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            keys.add(new Key(random.nextInt(3), random.nextInt(REFERENCES), random.nextInt(1_000_000) - 1000,
                    random.nextInt(2)));
        }
        return keys;
    }

    private static PairedReadEndKeys build(final List<Key> keys, final int maxKeysInRam, final int chunkSize,
                                           final List<Key> handedOver) {
        final PairedReadEndKeys set = new PairedReadEndKeys(REFERENCES, maxKeysInRam,
                end -> handedOver.add(Key.of(end)), chunkSize);
        for (final Key key : keys) {
            Assert.assertTrue(set.add(key.toEnd()));
        }
        set.doneAdding();
        return set;
    }

    /** Queries in ascending order, each added key with its neighbours, then all of them in a random order. */
    private static List<Key> queries(final Random random, final Set<Key> added) {
        final TreeSet<Key> ascending = new TreeSet<>();
        for (final Key key : added) {
            ascending.add(key);
            ascending.add(new Key(key.library, key.reference, key.coordinate - 1, key.orientation));
            ascending.add(new Key(key.library, key.reference, key.coordinate, 1 - key.orientation));
            ascending.add(new Key(key.library, key.reference + 1, key.coordinate, key.orientation));
        }
        final List<Key> queries = new ArrayList<>(ascending);
        final List<Key> shuffled = new ArrayList<>(ascending);
        Collections.shuffle(shuffled, random);
        queries.addAll(shuffled);
        queries.add(new Key(0, 0, Integer.MIN_VALUE, 0));
        queries.add(new Key(2, REFERENCES - 1, Integer.MAX_VALUE, 1));
        return queries;
    }

    @DataProvider(name = "keySets")
    public Object[][] keySets() {
        // Chunk sizes that make sorted runs of 32,767, 32,768, 32,769 and 65,537 keys, and the default 2^20.
        return new Object[][]{
                {"sorted, chunks of 32767", 32_767, 1, true, 150_000},
                {"sorted, chunks of 32768", 32_768, 1, true, 150_000},
                {"sorted, chunks of 32769", 32_769, 1, true, 150_000},
                {"sorted, chunks of 65537", 65_537, 1, true, 150_000},
                {"sorted, default chunks", 1 << 20, 1, true, 150_000},
                {"sorted, two libraries, chunks of 32769", 32_769, 2, true, 150_000},
                {"random, chunks of 32768", 32_768, 1, false, 120_000},
                {"random, default chunks", 1 << 20, 1, false, 120_000},
                {"sorted, chunks of 7", 7, 1, true, 5_000},
        };
    }

    @Test(dataProvider = "keySets")
    public void testContainsAnswersAsASetDoes(final String name, final int chunkSize, final int libraries,
                                             final boolean sorted, final int count) {
        final Random random = new Random(42 + chunkSize + count);
        final List<Key> keys = sorted ? sortedKeys(random, count, libraries) : randomKeys(random, count);
        final Set<Key> added = new HashSet<>(keys);
        final List<Key> handedOver = new ArrayList<>();
        final PairedReadEndKeys set = build(keys, Integer.MAX_VALUE, chunkSize, handedOver);
        Assert.assertTrue(handedOver.isEmpty());
        for (final Key query : queries(random, added)) {
            Assert.assertEquals(set.contains(query.toEnd()), added.contains(query), name + ": " + query);
        }
    }

    @Test(dataProvider = "keySets")
    public void testSplitRunsAnswersAsBefore(final String name, final int chunkSize, final int libraries,
                                             final boolean sorted, final int count) {
        final Random random = new Random(7 + chunkSize + count);
        final List<Key> keys = sorted ? sortedKeys(random, count, libraries) : randomKeys(random, count);
        final Set<Key> added = new HashSet<>(keys);
        final List<Key> queries = queries(random, added);
        for (final int maxKeys : new int[]{PairedReadEndKeys.SPLIT_KEYS, 1, 32_769}) {
            final PairedReadEndKeys whole = build(keys, Integer.MAX_VALUE, chunkSize, new ArrayList<>());
            final PairedReadEndKeys split = build(keys, Integer.MAX_VALUE, chunkSize, new ArrayList<>());
            // Split part-way through an ascending pass, so that the split moves a cursor that is in use.
            final int splitAt = queries.size() / 5;
            for (int i = 0; i < queries.size(); i++) {
                if (i == splitAt) {
                    if (maxKeys == PairedReadEndKeys.SPLIT_KEYS) {
                        split.splitRuns();
                    } else {
                        split.splitRuns(maxKeys);
                    }
                }
                final ReadEndsForMarkDuplicates end = queries.get(i).toEnd();
                final boolean expected = added.contains(queries.get(i));
                Assert.assertEquals(whole.contains(end), expected, name + ": " + queries.get(i));
                Assert.assertEquals(split.contains(end), expected,
                        name + ", split into " + maxKeys + ": " + queries.get(i));
            }
        }
    }

    @Test
    public void testSplitRunsBeforeAnyQueryAndTwice() {
        final Random random = new Random(11);
        final List<Key> keys = sortedKeys(random, 100_000, 1);
        final Set<Key> added = new HashSet<>(keys);
        final PairedReadEndKeys set = build(keys, Integer.MAX_VALUE, 32_768, new ArrayList<>());
        set.splitRuns();
        set.splitRuns(5);
        for (final Key query : queries(random, added)) {
            Assert.assertEquals(set.contains(query.toEnd()), added.contains(query), query.toString());
        }
    }

    @Test
    public void testSplitRunsOfAnEmptySet() {
        final PairedReadEndKeys set = build(Collections.emptyList(), Integer.MAX_VALUE, 32_768, new ArrayList<>());
        set.splitRuns();
        Assert.assertFalse(set.contains(new Key(0, 0, 1, 0).toEnd()));
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testSplitRunsBeforeDoneAdding() {
        new PairedReadEndKeys(REFERENCES, 100, end -> { }).splitRuns();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testContainsBeforeDoneAdding() {
        new PairedReadEndKeys(REFERENCES, 100, end -> { }).contains(new Key(0, 0, 1, 0).toEnd());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddAfterDoneAdding() {
        final PairedReadEndKeys set = new PairedReadEndKeys(REFERENCES, 100, end -> { });
        set.doneAdding();
        set.add(new Key(0, 0, 1, 0).toEnd());
    }

    @DataProvider(name = "unpackable")
    public Object[][] unpackable() {
        // 4 references take 2 bits, which leave 28 for the library; 2^30 references leave none, 2^31 - 1 fewer.
        return new Object[][]{
                {4, new Key(-1, 0, 1, ReadEnds.F)},
                {4, new Key(0, -1, 1, ReadEnds.F)},
                {4, new Key(0, 4, 1, ReadEnds.F)},
                {4, new Key(0, 0, 1, ReadEnds.FR)},
                {4, new Key(0, 0, 1, ReadEnds.RF)},
                {4, new Key(0, 0, 1, ReadEnds.FF)},
                {4, new Key(0, 0, 1, ReadEnds.RR)},
                {1 << 30, new Key(1, 0, 1, ReadEnds.F)},
                {Integer.MAX_VALUE, new Key(0, 0, 1, ReadEnds.F)},
        };
    }

    @Test(dataProvider = "unpackable")
    public void testAddRefusesWhatCannotBePacked(final int references, final Key key) {
        final List<Key> handedOver = new ArrayList<>();
        final PairedReadEndKeys set = new PairedReadEndKeys(references, 100, end -> handedOver.add(Key.of(end)));
        Assert.assertFalse(set.add(key.toEnd()));
        set.doneAdding();
        set.splitRuns();
        Assert.assertFalse(set.contains(key.toEnd()));
        Assert.assertTrue(handedOver.isEmpty());
    }

    @DataProvider(name = "packable")
    public Object[][] packable() {
        return new Object[][]{
                {4, new Key(0, 3, 1, ReadEnds.R)},
                {4, new Key(Short.MAX_VALUE, 0, Integer.MIN_VALUE, ReadEnds.F)},
                {4, new Key(0, 0, Integer.MAX_VALUE, ReadEnds.R)},
                {1, new Key(Short.MAX_VALUE, 0, -1, ReadEnds.F)},
                {1 << 30, new Key(0, (1 << 30) - 1, 0, ReadEnds.R)},
        };
    }

    @Test(dataProvider = "packable")
    public void testAddTakesWhatFits(final int references, final Key key) {
        final PairedReadEndKeys set = new PairedReadEndKeys(references, 100, end -> { });
        Assert.assertTrue(set.add(key.toEnd()));
        set.doneAdding();
        Assert.assertTrue(set.contains(key.toEnd()));
        Assert.assertFalse(set.contains(new Key(key.library, key.reference, key.coordinate,
                1 - key.orientation).toEnd()));
        set.splitRuns();
        Assert.assertTrue(set.contains(key.toEnd()));
    }

    /**
     * The overflow rule as the constructor documents it: a key repeated at once is not held again; when the keys held
     * reach maxKeysInRam and more than half of them are distinct, each distinct key is handed over, in key order, and
     * dropped; otherwise only the distinct keys stay held.
     */
    private static List<Key> expectedHandedOver(final List<Key> keys, final int maxKeysInRam, final Set<Key> kept) {
        final List<Key> handedOver = new ArrayList<>();
        final TreeSet<Key> distinct = new TreeSet<>();
        int held = 0;
        Key last = null;
        for (final Key key : keys) {
            if (held > 0 && key.equals(last)) {
                continue;
            }
            if (held == maxKeysInRam) {
                if (distinct.size() > maxKeysInRam / 2) {
                    handedOver.addAll(distinct);
                    distinct.clear();
                }
                held = distinct.size();
            }
            distinct.add(key);
            held++;
            last = key;
        }
        kept.addAll(distinct);
        return handedOver;
    }

    @DataProvider(name = "overflow")
    public Object[][] overflow() {
        return new Object[][]{
                {"distinct keys, chunks of 32767", 32_767, 100_000, 1_000_000, 350_000},
                {"distinct keys, chunks of 32768", 32_768, 100_000, 1_000_000, 350_000},
                {"distinct keys, chunks of 32769", 32_769, 100_000, 1_000_000, 350_000},
                {"distinct keys, default chunks", 1 << 20, 100_000, 1_000_000, 350_000},
                {"distinct keys, held whole", 1 << 20, 80_000, 1_000_000, 350_000},
                {"few distinct keys, chunks of 32769", 32_769, 100_000, 5_000, 350_000},
                {"some distinct keys, chunks of 1000", 1000, 20_000, 15_000, 200_000},
                {"tiny", 3, 10, 12, 1000},
        };
    }

    @Test
    public void testOverflowByHand() {
        final Key a = new Key(0, 0, 10, 0), b = new Key(0, 0, 20, 0), c = new Key(0, 1, 5, 1), d = new Key(0, 0, 15, 1),
                e = new Key(0, 0, 1, 0);
        final List<Key> handedOver = new ArrayList<>();
        final PairedReadEndKeys set = new PairedReadEndKeys(REFERENCES, 4, end -> handedOver.add(Key.of(end)), 2);
        // A, B, A, B fill the 4 keys held with 2 distinct, not more than half, so C keeps them; D brings 4 held again,
        // and E finds 4 distinct, which go over in key order.
        for (final Key key : new Key[]{a, b, a, b, c}) {
            set.add(key.toEnd());
        }
        Assert.assertTrue(handedOver.isEmpty());
        set.add(d.toEnd());
        set.add(e.toEnd());
        Assert.assertEquals(handedOver, List.of(a, d, b, c));
        Assert.assertEquals(set.keysHandedOver(), 4);
        set.doneAdding();
        Assert.assertTrue(set.contains(e.toEnd()));
        for (final Key key : new Key[]{a, b, c, d}) {
            Assert.assertFalse(set.contains(key.toEnd()));
        }
    }

    @Test(dataProvider = "overflow")
    public void testOverflowHandsOverAsDocumented(final String name, final int chunkSize, final int maxKeysInRam,
                                                  final int pool, final int count) {
        final Random random = new Random(3 + chunkSize + maxKeysInRam + pool);
        final List<Key> keys = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            keys.add(i > 0 && random.nextInt(6) == 0 ? keys.get(i - 1) :
                    new Key(0, random.nextInt(2), random.nextInt(pool), random.nextInt(2)));
        }
        final Set<Key> kept = new HashSet<>();
        final List<Key> expected = expectedHandedOver(keys, maxKeysInRam, kept);

        final List<Key> handedOver = new ArrayList<>();
        final List<ReadEndsForMarkDuplicates> ends = new ArrayList<>();
        final PairedReadEndKeys set = new PairedReadEndKeys(REFERENCES, maxKeysInRam, end -> {
            handedOver.add(Key.of(end));
            ends.add(end);
        }, chunkSize);
        for (final Key key : keys) {
            Assert.assertTrue(set.add(key.toEnd()));
        }
        set.doneAdding();
        Assert.assertEquals(handedOver, expected, name);
        Assert.assertEquals(set.keysHandedOver(), expected.size(), name);
        for (final ReadEndsForMarkDuplicates end : ends) {
            Assert.assertTrue(end.isPaired());
            Assert.assertEquals(end.read2ReferenceIndex, end.read1ReferenceIndex);
        }
        final Set<Key> queried = new HashSet<>(keys);
        queried.addAll(randomKeys(random, 10_000));
        for (final Key query : queried) {
            Assert.assertEquals(set.contains(query.toEnd()), kept.contains(query), name + ": " + query);
        }
        set.splitRuns();
        for (final Key query : queried) {
            Assert.assertEquals(set.contains(query.toEnd()), kept.contains(query), name + ", split: " + query);
        }
    }
}
