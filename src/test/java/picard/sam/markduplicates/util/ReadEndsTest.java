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

/**
 * Tests for {@link ReadEnds}: the orientation constants (written to disk by the codecs and ordered by the
 * duplicate-marking comparators), {@link ReadEnds#getOrientationByte}, {@link ReadEnds#isPaired()}, the state of
 * new read ends, and the read group and library accessors.
 */
public class ReadEndsTest {

    @Test
    public void testOrientationConstants() {
        Assert.assertEquals(ReadEnds.F, 0);
        Assert.assertEquals(ReadEnds.R, 1);
        Assert.assertEquals(ReadEnds.FF, 2);
        Assert.assertEquals(ReadEnds.FR, 3);
        Assert.assertEquals(ReadEnds.RR, 4);
        Assert.assertEquals(ReadEnds.RF, 5);
    }

    @DataProvider(name = "orientationBytes")
    public Object[][] orientationBytes() {
        return new Object[][]{
                {false, false, ReadEnds.FF},
                {false, true, ReadEnds.FR},
                {true, true, ReadEnds.RR},
                {true, false, ReadEnds.RF}
        };
    }

    @Test(dataProvider = "orientationBytes")
    public void testGetOrientationByte(final boolean read1NegativeStrand, final boolean read2NegativeStrand, final byte expected) {
        Assert.assertEquals(ReadEnds.getOrientationByte(read1NegativeStrand, read2NegativeStrand), expected);
    }

    @DataProvider(name = "pairedness")
    public Object[][] pairedness() {
        return new Object[][]{
                {-1, -1, false},
                {-1, 500, false},
                {0, -1, true},
                {0, 500, true},
                {3, 500, true},
                {Integer.MAX_VALUE, 0, true},
                {-2, 500, true}
        };
    }

    /** Pairedness is decided by the mate's reference index alone, not its coordinate. */
    @Test(dataProvider = "pairedness")
    public void testIsPaired(final int read2ReferenceIndex, final int read2Coordinate, final boolean expected) {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.read2ReferenceIndex = read2ReferenceIndex;
        ends.read2Coordinate = read2Coordinate;
        Assert.assertEquals(ends.isPaired(), expected);
    }

    /**
     * New read ends are unpaired until a mate is recorded, have no location until a read name parses, and are not
     * optical duplicates, which is what MarkDuplicates reports for every duplicate when READ_NAME_REGEX is null.
     */
    @Test
    public void testNewReadEnds() {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        Assert.assertFalse(ends.isPaired());
        Assert.assertFalse(ends.hasLocation());
        Assert.assertFalse(ends.isOpticalDuplicate);
    }

    @DataProvider(name = "shortValues")
    public Object[][] shortValues() {
        return new Object[][]{{(short) 0}, {(short) 1}, {(short) -1}, {Short.MAX_VALUE}, {Short.MIN_VALUE}};
    }

    @Test(dataProvider = "shortValues")
    public void testReadGroupAccessors(final short value) {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.libraryId = 9;
        ends.setReadGroup(value);
        Assert.assertEquals(ends.getReadGroup(), value);
        Assert.assertEquals(ends.readGroup, value);
        Assert.assertEquals(ends.libraryId, 9);

        ends.readGroup = 17;
        Assert.assertEquals(ends.getReadGroup(), 17);
    }

    @Test(dataProvider = "shortValues")
    public void testLibraryIdAccessors(final short value) {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.readGroup = 9;
        ends.setLibraryId(value);
        Assert.assertEquals(ends.getLibraryId(), value);
        Assert.assertEquals(ends.libraryId, value);
        Assert.assertEquals(ends.readGroup, 9);

        ends.libraryId = 23;
        Assert.assertEquals(ends.getLibraryId(), 23);
    }

    @Test
    public void testLocationSettersWriteFields() {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.setTile((short) 2203);
        ends.setX(Integer.MAX_VALUE);
        ends.setY(Integer.MIN_VALUE);
        Assert.assertEquals(ends.tile, 2203);
        Assert.assertEquals(ends.x, Integer.MAX_VALUE);
        Assert.assertEquals(ends.y, Integer.MIN_VALUE);
        Assert.assertTrue(ends.hasLocation());
    }
}
