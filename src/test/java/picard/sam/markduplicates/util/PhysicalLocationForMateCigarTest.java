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

import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.util.HashSet;
import java.util.Set;

/**
 * Tests for {@link PhysicalLocationForMateCigar} and {@link PhysicalLocationForMateCigarSet}: copying a location,
 * equality and hashing over library, read group, tile, x and y, and the set's add, remove, replace and size, which
 * keep one read end per physical location.
 */
public class PhysicalLocationForMateCigarTest {

    private static ReadEndsForMarkDuplicates location(final int libraryId, final int readGroup, final int tile, final int x, final int y) {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.libraryId = (short) libraryId;
        ends.readGroup = (short) readGroup;
        ends.tile = (short) tile;
        ends.x = x;
        ends.y = y;
        return ends;
    }

    @Test
    public void testCopiesLocation() {
        final PhysicalLocationForMateCigar loc = new PhysicalLocationForMateCigar(location(2, 3, 1101, -5, Integer.MAX_VALUE));
        Assert.assertEquals(loc.getLibraryId(), 2);
        Assert.assertEquals(loc.getReadGroup(), 3);
        Assert.assertEquals(loc.getTile(), 1101);
        Assert.assertEquals(loc.getX(), -5);
        Assert.assertEquals(loc.getY(), Integer.MAX_VALUE);

        final PhysicalLocationForMateCigar copy = new PhysicalLocationForMateCigar(loc);
        Assert.assertEquals(copy, loc);
        Assert.assertEquals(copy.getLibraryId(), 2);
        Assert.assertEquals(copy.getReadGroup(), 3);
        Assert.assertEquals(copy.getTile(), 1101);
        Assert.assertEquals(copy.getX(), -5);
        Assert.assertEquals(copy.getY(), Integer.MAX_VALUE);
    }

    @Test
    public void testSetters() {
        final PhysicalLocationForMateCigar loc = new PhysicalLocationForMateCigar(location(0, 0, 0, 0, 0));
        loc.setLibraryId(Short.MAX_VALUE);
        loc.setReadGroup(Short.MIN_VALUE);
        loc.setTile((short) 7);
        loc.setX(8);
        loc.setY(9);
        Assert.assertEquals(loc.getLibraryId(), Short.MAX_VALUE);
        Assert.assertEquals(loc.getReadGroup(), Short.MIN_VALUE);
        Assert.assertEquals(loc.getTile(), 7);
        Assert.assertEquals(loc.getX(), 8);
        Assert.assertEquals(loc.getY(), 9);
    }

    @DataProvider(name = "pairs")
    public Object[][] pairs() {
        return new Object[][]{
                {location(1, 2, 1101, 300, 400), location(1, 2, 1101, 300, 400), true},
                {location(1, 2, 1101, 300, 400), location(2, 2, 1101, 300, 400), false},
                {location(1, 2, 1101, 300, 400), location(1, 3, 1101, 300, 400), false},
                {location(1, 2, 1101, 300, 400), location(1, 2, 1102, 300, 400), false},
                {location(1, 2, 1101, 300, 400), location(1, 2, 1101, 301, 400), false},
                {location(1, 2, 1101, 300, 400), location(1, 2, 1101, 300, 399), false},
                {location(1, 2, 1101, 300, 400), location(1, 2, 1101, 400, 300), false},
                {location(1, 2, 1101, Integer.MIN_VALUE, 0), location(1, 2, 1101, 0, 0), false},
                {location(1, 2, 1101, 0, Integer.MIN_VALUE), location(1, 2, 1101, 0, Integer.MAX_VALUE), false},
                {location(1, 2, 1101, Integer.MAX_VALUE, Integer.MIN_VALUE), location(1, 2, 1101, Integer.MAX_VALUE, Integer.MIN_VALUE), true},
                {location(Short.MIN_VALUE, Short.MAX_VALUE, Short.MIN_VALUE, 1, 1), location(Short.MAX_VALUE, Short.MAX_VALUE, Short.MIN_VALUE, 1, 1), false},
                {location(1, -1, -1, -1, -1), location(1, -1, -1, -1, -1), true}
        };
    }

    @Test(dataProvider = "pairs")
    public void testEqualsAndHashCode(final ReadEndsForMarkDuplicates a, final ReadEndsForMarkDuplicates b, final boolean expected) {
        final PhysicalLocationForMateCigar locA = new PhysicalLocationForMateCigar(a);
        final PhysicalLocationForMateCigar locB = new PhysicalLocationForMateCigar(b);
        Assert.assertEquals(locA.equals(locB), expected);
        Assert.assertEquals(locB.equals(locA), expected);
        Assert.assertTrue(locA.equals(locA));
        if (expected) {
            Assert.assertEquals(locA.hashCode(), locB.hashCode());
        }
        final Set<PhysicalLocationForMateCigar> set = new HashSet<>();
        set.add(locA);
        Assert.assertEquals(set.contains(locB), expected);
    }

    @Test
    public void testNotEqualToOtherTypes() {
        final ReadEndsForMarkDuplicates ends = location(1, 2, 1101, 300, 400);
        final PhysicalLocationForMateCigar loc = new PhysicalLocationForMateCigar(ends);
        Assert.assertFalse(loc.equals(ends));
        Assert.assertFalse(loc.equals(null));
        Assert.assertFalse(loc.equals("1:2:1101:300:400"));
    }

    private static ReadEndsForMateCigar mateCigarEnds(final MateCigarTestRecords records, final String name,
                                                      final int libraryId, final int readGroup, final int tile, final int x, final int y) {
        final SAMRecord rec = records.pairEnd(name, "rgA", 0, 100, false, "10M", true, 0, 300, true, "10M", 30);
        final ReadEndsForMateCigar ends = new ReadEndsForMateCigar(records.header, new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, 0),
                new OpticalDuplicateFinder(), (short) libraryId);
        ends.readGroup = (short) readGroup;
        ends.tile = (short) tile;
        ends.x = x;
        ends.y = y;
        return ends;
    }

    private static void assertSet(final PhysicalLocationForMateCigarSet set, final ReadEnds... expected) {
        Assert.assertEquals(set.size(), expected.length);
        Assert.assertEquals(set.getReadEnds().size(), expected.length);
        for (final ReadEnds ends : expected) {
            Assert.assertTrue(set.getReadEnds().contains(ends));
        }
    }

    @Test
    public void testSetKeepsOneEndPerLocation() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final ReadEndsForMateCigar a = mateCigarEnds(records, "a", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar sameAsA = mateCigarEnds(records, "a2", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar otherX = mateCigarEnds(records, "b", 1, 0, 1101, 11, 20);
        final ReadEndsForMateCigar otherLibrary = mateCigarEnds(records, "c", 2, 0, 1101, 10, 20);
        final ReadEndsForMateCigar otherReadGroup = mateCigarEnds(records, "d", 1, 1, 1101, 10, 20);
        final ReadEndsForMateCigar otherTile = mateCigarEnds(records, "e", 1, 0, 1102, 10, 20);

        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        assertSet(set);
        set.add(a);
        assertSet(set, a);
        set.add(a);
        assertSet(set, a);
        set.add(sameAsA);
        assertSet(set, a);
        Assert.assertFalse(set.getReadEnds().contains(sameAsA));
        set.add(otherX);
        set.add(otherLibrary);
        set.add(otherReadGroup);
        set.add(otherTile);
        assertSet(set, a, otherX, otherLibrary, otherReadGroup, otherTile);
    }

    @Test
    public void testSetRemove() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final ReadEndsForMateCigar a = mateCigarEnds(records, "a", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar b = mateCigarEnds(records, "b", 1, 0, 1101, 30, 40);
        final ReadEndsForMateCigar absent = mateCigarEnds(records, "c", 1, 0, 1101, 50, 60);

        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        set.add(a);
        set.add(b);
        set.remove(absent);
        assertSet(set, a, b);
        set.remove(a);
        assertSet(set, b);
        set.remove(a);
        assertSet(set, b);
        set.remove(b);
        assertSet(set);

        // a location freed by remove can be taken again
        set.add(a);
        assertSet(set, a);
    }

    @Test
    public void testSetReplace() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final ReadEndsForMateCigar current = mateCigarEnds(records, "a", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar other = mateCigarEnds(records, "b", 1, 0, 1101, 30, 40);
        final ReadEndsForMateCigar third = mateCigarEnds(records, "c", 1, 0, 1101, 50, 60);

        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        set.add(current);
        set.add(third);
        set.replace(current, other);
        assertSet(set, other, third);
        Assert.assertFalse(set.getReadEnds().contains(current));
    }

    /** The mate of a pair shares its read name and so its location: replacing one end with the other keeps one entry. */
    @Test
    public void testSetReplaceWithMateAtSameLocation() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final ReadEndsForMateCigar end1 = mateCigarEnds(records, "p", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar end2 = mateCigarEnds(records, "p", 1, 0, 1101, 10, 20);
        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        set.add(end1);
        set.replace(end1, end2);
        assertSet(set, end2);
        Assert.assertFalse(set.getReadEnds().contains(end1));
    }

    @Test
    public void testSetReplaceWithOtherAtAnOccupiedLocation() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final ReadEndsForMateCigar current = mateCigarEnds(records, "a", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar third = mateCigarEnds(records, "c", 1, 0, 1101, 50, 60);
        final ReadEndsForMateCigar other = mateCigarEnds(records, "b", 1, 0, 1101, 50, 60);

        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        set.add(current);
        set.add(third);
        set.replace(current, other);
        assertSet(set, third);
    }

    @Test
    public void testSetReplaceRequiresCurrent() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final ReadEndsForMateCigar a = mateCigarEnds(records, "a", 1, 0, 1101, 10, 20);
        final ReadEndsForMateCigar b = mateCigarEnds(records, "b", 1, 0, 1101, 30, 40);
        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        set.add(a);
        final PicardException e = Assert.expectThrows(PicardException.class, () -> set.replace(b, a));
        Assert.assertTrue(e.getMessage().contains("not in the set"), e.getMessage());
        assertSet(set, a);
    }
}
