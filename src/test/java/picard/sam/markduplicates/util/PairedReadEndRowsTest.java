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
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import java.util.function.Consumer;

/**
 * Tests for {@link PairedReadEndRows}: every field of an end comes back from {@code get} as it was added; {@code add}
 * refuses what its Javadoc lists; {@code drainTo} hands the rows over in the order added; and {@code sort} orders the
 * rows as MarkDuplicates' comparator without barcodes orders ends (library id, read1 reference index and coordinate,
 * orientation, read2 reference index and coordinate, tile, x, y, read1 and read2 index in the file, then the order
 * added), marks where each duplicate set starts, and {@code getAt} gives each position's end.
 */
public class PairedReadEndRowsTest {
    private static final int LOW = -(1 << 30);
    private static final int HIGH = (1 << 30) - 1;

    private static ReadEndsForMarkDuplicates end(final Random random, final int coordinates, final boolean wide) {
        final ReadEndsForMarkDuplicates end = new ReadEndsForMarkDuplicates();
        end.libraryId = (short) random.nextInt(3);
        end.read1ReferenceIndex = random.nextInt(3);
        end.read2ReferenceIndex = end.read1ReferenceIndex + random.nextInt(2);
        if (wide) {
            end.read1Coordinate = LOW + random.nextInt(Integer.MAX_VALUE);
            end.read2Coordinate = random.nextInt(4) == 0 ? end.read1Coordinate : LOW + random.nextInt(Integer.MAX_VALUE);
            end.x = LOW + random.nextInt(Integer.MAX_VALUE);
            end.y = LOW + random.nextInt(Integer.MAX_VALUE);
        } else {
            end.read1Coordinate = random.nextInt(coordinates);
            end.read2Coordinate = end.read1Coordinate + random.nextInt(3);
            end.x = random.nextInt(4);
            end.y = random.nextInt(4);
        }
        end.orientation = (byte) (ReadEnds.FF + random.nextInt(4));
        end.orientationForOpticalDuplicates = (byte) (ReadEnds.FF + random.nextInt(4));
        end.tile = (short) (random.nextInt(3) - 1);
        end.readGroup = (short) (random.nextInt(4) - 1);
        end.score = (short) (random.nextInt(2000) - 100);
        end.read1IndexInFile = random.nextInt(4) == 0 ? 5 : random.nextInt(1 << 20) + ((long) random.nextInt(3) << 32);
        end.read2IndexInFile = random.nextInt(4) == 0 ? end.read1IndexInFile : end.read1IndexInFile + random.nextInt(100);
        return end;
    }

    private static void assertSameFields(final ReadEndsForMarkDuplicates actual, final ReadEndsForMarkDuplicates expected,
                                         final String message) {
        Assert.assertEquals(actual.getClass(), ReadEndsForMarkDuplicates.class, message);
        Assert.assertEquals(actual.libraryId, expected.libraryId, message);
        Assert.assertEquals(actual.read1ReferenceIndex, expected.read1ReferenceIndex, message);
        Assert.assertEquals(actual.read1Coordinate, expected.read1Coordinate, message);
        Assert.assertEquals(actual.read2ReferenceIndex, expected.read2ReferenceIndex, message);
        Assert.assertEquals(actual.read2Coordinate, expected.read2Coordinate, message);
        Assert.assertEquals(actual.orientation, expected.orientation, message);
        Assert.assertEquals(actual.orientationForOpticalDuplicates, expected.orientationForOpticalDuplicates, message);
        Assert.assertEquals(actual.tile, expected.tile, message);
        Assert.assertEquals(actual.x, expected.x, message);
        Assert.assertEquals(actual.y, expected.y, message);
        Assert.assertEquals(actual.readGroup, expected.readGroup, message);
        Assert.assertEquals(actual.score, expected.score, message);
        Assert.assertEquals(actual.read1IndexInFile, expected.read1IndexInFile, message);
        Assert.assertEquals(actual.read2IndexInFile, expected.read2IndexInFile, message);
        Assert.assertEquals(actual.duplicateSetSize, -1, message);
        Assert.assertFalse(actual.isOpticalDuplicate, message);
    }

    @Test
    public void testGetGivesBackEveryField() {
        final Random random = new Random(1);
        final List<ReadEndsForMarkDuplicates> ends = new ArrayList<>();
        final PairedReadEndRows rows = new PairedReadEndRows(100_000);
        for (int i = 0; i < 100_000; i++) {
            final ReadEndsForMarkDuplicates end = end(random, 1000, i % 2 == 0);
            ends.add(end);
            Assert.assertTrue(rows.add(end));
        }
        Assert.assertEquals(rows.size(), ends.size());
        for (int row = 0; row < ends.size(); row++) {
            assertSameFields(rows.get(row), ends.get(row), "row " + row);
        }
    }

    @DataProvider(name = "refused")
    public Object[][] refused() {
        final List<Object[]> cases = new ArrayList<>();
        cases.add(new Object[]{"barcodes", (Consumer<ReadEndsForMarkDuplicates>) end -> { }, true});
        cases.add(new Object[]{"duplicate set size", (Consumer<ReadEndsForMarkDuplicates>) end -> end.duplicateSetSize = 2, false});
        cases.add(new Object[]{"optical duplicate", (Consumer<ReadEndsForMarkDuplicates>) end -> end.isOpticalDuplicate = true, false});
        final String[] fields = {"read1ReferenceIndex", "read1Coordinate", "read2ReferenceIndex", "read2Coordinate", "x", "y"};
        for (final String field : fields) {
            for (final int value : new int[]{LOW - 1, HIGH + 1, Integer.MIN_VALUE, Integer.MAX_VALUE}) {
                cases.add(new Object[]{field + " " + value, (Consumer<ReadEndsForMarkDuplicates>) end -> set(end, field, value), false});
            }
        }
        return cases.toArray(new Object[0][]);
    }

    private static void set(final ReadEndsForMarkDuplicates end, final String field, final int value) {
        switch (field) {
            case "read1ReferenceIndex": end.read1ReferenceIndex = value; break;
            case "read1Coordinate": end.read1Coordinate = value; break;
            case "read2ReferenceIndex": end.read2ReferenceIndex = value; break;
            case "read2Coordinate": end.read2Coordinate = value; break;
            case "x": end.x = value; break;
            case "y": end.y = value; break;
            default: throw new IllegalArgumentException(field);
        }
    }

    @Test(dataProvider = "refused")
    public void testAddRefuses(final String name, final Consumer<ReadEndsForMarkDuplicates> change, final boolean barcodes) {
        final PairedReadEndRows rows = new PairedReadEndRows(10);
        final ReadEndsForMarkDuplicates end = barcodes ? new ReadEndsForMarkDuplicatesWithBarcodes() :
                new ReadEndsForMarkDuplicates();
        end.read1ReferenceIndex = 0;
        end.read2ReferenceIndex = 0;
        end.read1Coordinate = 100;
        end.read2Coordinate = 200;
        change.accept(end);
        Assert.assertFalse(rows.add(end), name);
        Assert.assertEquals(rows.size(), 0, name);
    }

    @Test
    public void testAddTakesTheEdgesOfTheRange() {
        final String[] fields = {"read1ReferenceIndex", "read1Coordinate", "read2ReferenceIndex", "read2Coordinate", "x", "y"};
        final PairedReadEndRows rows = new PairedReadEndRows(100);
        final List<ReadEndsForMarkDuplicates> ends = new ArrayList<>();
        for (final String field : fields) {
            for (final int value : new int[]{LOW, HIGH}) {
                final ReadEndsForMarkDuplicates end = new ReadEndsForMarkDuplicates();
                set(end, field, value);
                Assert.assertTrue(rows.add(end), field + " " + value);
                ends.add(end);
            }
        }
        for (int row = 0; row < ends.size(); row++) {
            assertSameFields(rows.get(row), ends.get(row), "row " + row);
        }
    }

    @Test
    public void testAddRefusesPastTheCapacity() {
        final PairedReadEndRows rows = new PairedReadEndRows(3);
        for (int i = 0; i < 3; i++) {
            Assert.assertTrue(rows.add(new ReadEndsForMarkDuplicates()));
        }
        Assert.assertFalse(rows.add(new ReadEndsForMarkDuplicates()));
        Assert.assertEquals(rows.size(), 3);
        Assert.assertFalse(new PairedReadEndRows(0).add(new ReadEndsForMarkDuplicates()));
    }

    @Test
    public void testDrainToHandsOverInTheOrderAdded() {
        final Random random = new Random(2);
        final PairedReadEndRows rows = new PairedReadEndRows(50_000);
        final List<ReadEndsForMarkDuplicates> ends = new ArrayList<>();
        for (int i = 0; i < 20_000; i++) {
            ends.add(end(random, 100, false));
            Assert.assertTrue(rows.add(ends.get(i)));
        }
        final List<ReadEndsForMarkDuplicates> drained = new ArrayList<>();
        rows.drainTo(drained::add);
        Assert.assertEquals(rows.size(), 0);
        Assert.assertEquals(drained.size(), ends.size());
        for (int i = 0; i < ends.size(); i++) {
            assertSameFields(drained.get(i), ends.get(i), "end " + i);
        }
        // Emptied, it takes ends again from row 0.
        final ReadEndsForMarkDuplicates again = end(random, 100, false);
        Assert.assertTrue(rows.add(again));
        assertSameFields(rows.get(0), again, "after drainTo");
    }

    private static final Comparator<ReadEndsForMarkDuplicates> ORDER =
            Comparator.<ReadEndsForMarkDuplicates>comparingInt(e -> e.libraryId)
                    .thenComparingInt(e -> e.read1ReferenceIndex)
                    .thenComparingInt(e -> e.read1Coordinate)
                    .thenComparingInt(e -> e.orientation)
                    .thenComparingInt(e -> e.read2ReferenceIndex)
                    .thenComparingInt(e -> e.read2Coordinate)
                    .thenComparingInt(e -> e.tile)
                    .thenComparingInt(e -> e.x)
                    .thenComparingInt(e -> e.y)
                    .thenComparingLong(e -> e.read1IndexInFile)
                    .thenComparingLong(e -> e.read2IndexInFile);

    private static boolean sameSet(final ReadEndsForMarkDuplicates a, final ReadEndsForMarkDuplicates b) {
        return a.libraryId == b.libraryId && a.read1ReferenceIndex == b.read1ReferenceIndex &&
                a.read1Coordinate == b.read1Coordinate && a.orientation == b.orientation &&
                a.read2ReferenceIndex == b.read2ReferenceIndex && a.read2Coordinate == b.read2Coordinate;
    }

    @DataProvider(name = "sorts")
    public Object[][] sorts() {
        // Over 2 * 65,536 rows the tied prefixes are sorted in more than one task; wide fields drop prefix bits.
        return new Object[][]{
                {"empty", 0, 10, false},
                {"one", 1, 10, false},
                {"few coordinates", 2_000, 20, false},
                {"many ties", 200_000, 3_000, false},
                {"wide fields", 150_000, 0, true},
        };
    }

    @Test(dataProvider = "sorts")
    public void testSortOrdersAndMarksSets(final String name, final int n, final int coordinates, final boolean wide) {
        final Random random = new Random(n + coordinates);
        final PairedReadEndRows rows = new PairedReadEndRows(n + 10);
        final List<ReadEndsForMarkDuplicates> ends = new ArrayList<>(n);
        final List<Integer> order = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            ends.add(end(random, coordinates, wide));
            // Some exact copies, so that only the order added tells them apart.
            if (i > 0 && random.nextInt(10) == 0) {
                final ReadEndsForMarkDuplicates copy = new ReadEndsForMarkDuplicates(ends.get(i - 1));
                copy.score = (short) (ends.get(i).score + 1);
                ends.set(i, copy);
            }
            Assert.assertTrue(rows.add(ends.get(i)));
            order.add(i);
        }
        // List.sort is stable: ends equal under the comparator stay in the order added.
        order.sort((a, b) -> ORDER.compare(ends.get(a), ends.get(b)));

        rows.sort();
        Assert.assertEquals(rows.size(), n);
        for (int position = 0; position < n; position++) {
            Assert.assertEquals(rows.rowAt(position), (int) order.get(position), name + " position " + position);
            final boolean starts = position == 0 ||
                    !sameSet(ends.get(order.get(position - 1)), ends.get(order.get(position)));
            Assert.assertEquals(rows.startsDuplicateSet(position), starts, name + " position " + position);
            if (position > 0) {
                Assert.assertEquals(rows.sameDuplicateSet(order.get(position - 1), order.get(position)), !starts);
            }
            assertSameFields(rows.getAt(position), ends.get(order.get(position)), name + " position " + position);
        }
        // Out of sequence, getAt still gives the end at the position.
        for (int i = 0; i < Math.min(n, 1000); i++) {
            final int position = random.nextInt(n);
            assertSameFields(rows.getAt(position), ends.get(order.get(position)), name + " position " + position);
        }
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddAfterSort() {
        final PairedReadEndRows rows = new PairedReadEndRows(10);
        rows.sort();
        rows.add(new ReadEndsForMarkDuplicates());
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testSortTwice() {
        final PairedReadEndRows rows = new PairedReadEndRows(10);
        rows.sort();
        rows.sort();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testDrainToAfterSort() {
        final PairedReadEndRows rows = new PairedReadEndRows(10);
        rows.add(new ReadEndsForMarkDuplicates());
        rows.sort();
        rows.drainTo(end -> { });
    }
}
