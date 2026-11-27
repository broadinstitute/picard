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

import java.util.Set;

/**
 * Tests for {@link ReadEndsForMateCigar}: the 5' unclipped coordinates of the read and of its mate (from the MC
 * tag), the pair orientation and the first-then-second orientation used for optical duplicates, reads whose mate
 * is unmapped, read groups kept apart and tile, x and y parsed from the read name, the copy constructor, and the
 * accessors for the record and the optical-duplicate location set.
 */
public class ReadEndsForMateCigarTest {

    private static ReadEndsForMateCigar ends(final MateCigarTestRecords records, final SAMRecord rec, final long ordinal, final short libraryId) {
        return new ReadEndsForMateCigar(records.header, new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, ordinal),
                new OpticalDuplicateFinder(), libraryId);
    }

    /**
     * A null orientationForOpticalDuplicates is not checked: it is used only for pairs with both ends mapped. The
     * read2 fields of a fragment are checked as -1 because MarkQueue finds a fragment at a pair's 5' end with a probe
     * whose read2 fields are -1.
     */
    private static void assertCoordinates(final ReadEndsForMateCigar ends, final int read1ReferenceIndex, final int read1Coordinate,
                                          final int read2ReferenceIndex, final int read2Coordinate, final byte orientation,
                                          final Byte orientationForOpticalDuplicates, final byte hasUnmapped,
                                          final boolean isPaired) {
        Assert.assertEquals(ends.read1ReferenceIndex, read1ReferenceIndex, "read1ReferenceIndex");
        Assert.assertEquals(ends.read1Coordinate, read1Coordinate, "read1Coordinate");
        Assert.assertEquals(ends.read2ReferenceIndex, read2ReferenceIndex, "read2ReferenceIndex");
        Assert.assertEquals(ends.read2Coordinate, read2Coordinate, "read2Coordinate");
        Assert.assertEquals(ends.orientation, orientation, "orientation");
        if (orientationForOpticalDuplicates != null) {
            Assert.assertEquals(ends.orientationForOpticalDuplicates, orientationForOpticalDuplicates.byteValue(),
                    "orientationForOpticalDuplicates");
        }
        Assert.assertEquals(ends.hasUnmapped, hasUnmapped, "hasUnmapped");
        Assert.assertEquals(ends.isPaired(), isPaired, "isPaired");
    }

    @DataProvider(name = "fragments")
    public Object[][] fragments() {
        return new Object[][]{
                // forward reads start at the unclipped start
                {false, 100, "10M", 100},
                {false, 100, "5S10M", 95},
                {false, 100, "2H3S10M", 95},
                {false, 100, "10M5S", 100},
                // reverse reads start at the unclipped end
                {true, 100, "10M", 109},
                {true, 100, "10M5S", 114},
                {true, 100, "5S10M", 109},
                {true, 100, "4M2D4M1I1M3S", 100 + 4 + 2 + 4 + 1 - 1 + 3},
                {true, 1, "10M", 10}
        };
    }

    @Test(dataProvider = "fragments")
    public void testFragment(final boolean negativeStrand, final int start, final String cigar, final int expectedCoordinate) {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.fragment("frag", "rgA", 1, start, negativeStrand, cigar, 30);
        final ReadEndsForMateCigar ends = ends(records, rec, 7, (short) 4);
        assertCoordinates(ends, 1, expectedCoordinate, -1, -1, negativeStrand ? ReadEnds.R : ReadEnds.F, null, (byte) 0,
                false);
        Assert.assertEquals(ends.libraryId, 4);
        Assert.assertEquals(ends.getLibraryId(), 4);
        Assert.assertSame(ends.getRecord(), rec);
        Assert.assertEquals(ends.getRecordReadName(), "frag");
        Assert.assertEquals(ends.getSamRecordIndex().getRecordOrdinal(), 7);
        Assert.assertSame(ends.getSamRecordIndex(), ends.samRecordWithOrdinal);
        Assert.assertNull(ends.getLocationSet());
    }

    /**
     * Each row is one end of a pair: its strand, whether it is the first end, and its mate's strand. The mate maps to
     * chr2:500 with cigar 20M4S, so its 5' end is 500 forward and 523 reversed; this end maps to chr1:100 with 3S20M,
     * so its 5' end is 97 forward and 119 reversed.
     */
    @DataProvider(name = "pairEnds")
    public Object[][] pairEnds() {
        return new Object[][]{
                {false, true, false, 97, 500, ReadEnds.FF, ReadEnds.FF},
                {false, true, true, 97, 523, ReadEnds.FR, ReadEnds.FR},
                {true, true, false, 119, 500, ReadEnds.RF, ReadEnds.RF},
                {true, true, true, 119, 523, ReadEnds.RR, ReadEnds.RR},
                // for the second end the optical orientation is still the first end's strand, then the second's
                {false, false, false, 97, 500, ReadEnds.FF, ReadEnds.FF},
                {false, false, true, 97, 523, ReadEnds.FR, ReadEnds.RF},
                {true, false, false, 119, 500, ReadEnds.RF, ReadEnds.FR},
                {true, false, true, 119, 523, ReadEnds.RR, ReadEnds.RR}
        };
    }

    @Test(dataProvider = "pairEnds")
    public void testPairEnd(final boolean negativeStrand, final boolean firstOfPair, final boolean mateNegativeStrand,
                            final int expectedRead1Coordinate, final int expectedRead2Coordinate,
                            final byte expectedOrientation, final byte expectedOpticalOrientation) {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.pairEnd("pair", "rgA", 0, 100, negativeStrand, "3S20M", firstOfPair,
                1, 500, mateNegativeStrand, "20M4S", 30);
        final ReadEndsForMateCigar ends = ends(records, rec, 0, (short) 1);
        assertCoordinates(ends, 0, expectedRead1Coordinate, 1, expectedRead2Coordinate, expectedOrientation,
                expectedOpticalOrientation, (byte) 0, true);
    }

    @DataProvider(name = "mateUnmapped")
    public Object[][] mateUnmapped() {
        return new Object[][]{{false, true}, {true, true}, {false, false}, {true, false}};
    }

    /** A read with an unmapped mate is a fragment for orientation and coordinates, but still counts as paired. */
    @Test(dataProvider = "mateUnmapped")
    public void testMateUnmapped(final boolean negativeStrand, final boolean firstOfPair) {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.mateUnmapped("halfMapped", "rgB", 0, 200, negativeStrand, "10M", firstOfPair, 30);
        final ReadEndsForMateCigar ends = ends(records, rec, 0, (short) 2);
        assertCoordinates(ends, 0, negativeStrand ? 209 : 200, -1, -1, negativeStrand ? ReadEnds.R : ReadEnds.F,
                null, (byte) 1, true);
    }

    @Test
    public void testUnmappedReadIsRejected() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.fragment("unmapped", "rgA", 0, 200, false, "10M", 30);
        rec.setReadUnmappedFlag(true);
        final PicardException e = Assert.expectThrows(PicardException.class, () -> ends(records, rec, 0, (short) 1));
        Assert.assertTrue(e.getMessage().contains("unmapped"), e.getMessage());
    }

    /** Tile, x and y parsed from the read name, whatever the read group; a name without them gives no location. */
    @DataProvider(name = "locations")
    public Object[][] locations() {
        return new Object[][]{
                {"M1:7:FC:2:1101:2000:3000", "rgA", 1101, 2000, 3000},
                {"M1:7:FC:2:1101:2000:3000", "rgB", 1101, 2000, 3000},
                {"M1:7:FC:2:2203:15:70000", "rgC", 2203, 15, 70000},
                {"M1:7:FC:2:1101:2000:3000", null, 1101, 2000, 3000},
                {"M1:7:FC:2:1101:2000:3000", "rgNotInHeader", 1101, 2000, 3000},
                {"noLocation", "rgB", -1, -1, -1}
        };
    }

    @Test(dataProvider = "locations")
    public void testLocation(final String name, final String readGroup, final int expectedTile, final int expectedX,
                             final int expectedY) {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.pairEnd(name, readGroup, 0, 100, false, "10M", true, 0, 300, true, "10M", 30);
        final ReadEndsForMateCigar ends = ends(records, rec, 0, (short) 1);
        Assert.assertEquals(ends.getTile(), expectedTile);
        Assert.assertEquals(ends.getX(), expectedX);
        Assert.assertEquals(ends.getY(), expectedY);
        Assert.assertEquals(ends.hasLocation(), expectedTile != -1);
    }

    /** Optical duplicates must share a read group: ends of one read group share an index that others do not. */
    @Test
    public void testReadGroupIndexSeparatesReadGroups() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final short[] indices = new short[4];
        final String[] readGroups = {"rgA", "rgB", "rgC", "rgA"};
        for (int i = 0; i < readGroups.length; i++) {
            final SAMRecord rec = records.pairEnd("M1:7:FC:2:1101:2000:" + (3000 + i), readGroups[i], 0, 100, false,
                    "10M", true, 0, 300, true, "10M", 30);
            indices[i] = ends(records, rec, i, (short) 1).getReadGroup();
        }
        Assert.assertEquals(indices[3], indices[0]);
        Assert.assertNotEquals(indices[1], indices[0]);
        Assert.assertNotEquals(indices[2], indices[0]);
        Assert.assertNotEquals(indices[2], indices[1]);
    }

    @Test
    public void testCopyConstructor() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.pairEnd("M1:7:FC:2:1101:2000:3000", "rgB", 0, 100, true, "10M", true, 1, 300, false, "10M", 30);
        final ReadEndsForMateCigar original = ends(records, rec, 3, (short) 2);
        original.setLocationSet(new PhysicalLocationForMateCigarSet());

        final SAMRecord other = records.fragment("other", "rgA", 1, 5, false, "10M", 20);
        final SamRecordWithOrdinalAndSetDuplicateReadFlag otherIndex = new SamRecordWithOrdinalAndSetDuplicateReadFlag(other, 9);
        final ReadEndsForMateCigar copy = new ReadEndsForMateCigar(original, otherIndex);

        Assert.assertEquals(copy.readGroup, original.readGroup);
        Assert.assertEquals(copy.tile, 1101);
        Assert.assertEquals(copy.x, 2000);
        Assert.assertEquals(copy.y, 3000);
        Assert.assertEquals(copy.read1ReferenceIndex, 0);
        Assert.assertEquals(copy.read1Coordinate, 109);
        Assert.assertEquals(copy.read2ReferenceIndex, 1);
        Assert.assertEquals(copy.read2Coordinate, 300);
        Assert.assertEquals(copy.hasUnmapped, 0);
        Assert.assertEquals(copy.orientation, ReadEnds.RF);
        Assert.assertEquals(copy.libraryId, 2);
        // the record is the one given, not the original's
        Assert.assertSame(copy.getSamRecordIndex(), otherIndex);
        Assert.assertSame(copy.getRecord(), other);
        Assert.assertEquals(copy.getRecordReadName(), "other");
        Assert.assertFalse(copy.isPaired());
    }

    @Test
    public void testLocationSetAccessors() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMRecord rec = records.pairEnd("p", "rgA", 0, 100, false, "10M", true, 0, 300, true, "10M", 30);
        final ReadEndsForMateCigar ends = ends(records, rec, 0, (short) 1);
        Assert.assertNull(ends.getLocationSet());
        Assert.assertNull(ends.removeLocationSet());

        final PhysicalLocationForMateCigarSet set = new PhysicalLocationForMateCigarSet();
        set.add(ends);
        ends.setLocationSet(set);
        Assert.assertSame(ends.getLocationSet(), set);
        Assert.assertSame(ends.getLocationSet(), set);
        Assert.assertSame(ends.removeLocationSet(), set);
        Assert.assertNull(ends.getLocationSet());

        ends.setLocationSet(set);
        final Set<ReadEnds> readEnds = ends.getReadEndSetForOpticalDuplicates();
        Assert.assertSame(readEnds, set.getReadEnds());
        Assert.assertEquals(readEnds.size(), 1);
        Assert.assertTrue(readEnds.contains(ends));
        Assert.assertNull(ends.getLocationSet());

        final PicardException e = Assert.expectThrows(PicardException.class, ends::getReadEndSetForOpticalDuplicates);
        Assert.assertTrue(e.getMessage().contains("Already called"), e.getMessage());
    }
}
