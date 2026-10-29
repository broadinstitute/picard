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

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.util.LibraryIdGenerator;
import picard.sam.markduplicates.util.ReadEnds;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicates;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Consumer;

/**
 * Tests the parts of {@link MarkDuplicates} that turn a record into read ends and group read ends into duplicate
 * sets: {@link MarkDuplicates#buildReadEnds}, {@link MarkDuplicates#getReadDuplicateScore}, the ordering of
 * {@link MarkDuplicates.ReadEndsMDComparator} without barcodes, {@link MarkDuplicates#areComparableForDuplicates}
 * without barcodes, and the codes of {@link MarkDuplicates.DuplicateType}.
 */
public class MarkDuplicatesReadEndsTest {
    private static final int CHR1 = 0, CHR3 = 2;

    private static MarkDuplicatesEndToEndRunner header() {
        return new MarkDuplicatesEndToEndRunner()
                .readGroup("rg1", "libA")
                .readGroup("rg2", "libB")
                .readGroup("rg3", "libA");
    }

    private static MarkDuplicates markDuplicates(final MarkDuplicatesEndToEndRunner input) {
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.libraryIdGenerator = new LibraryIdGenerator(input.header);
        return markDuplicates;
    }

    @Test
    public void testForwardFragment() {
        final MarkDuplicatesEndToEndRunner input = header();
        final SAMRecord rec = input.fragment("RUN:1:1101:30:40", "rg2", CHR1, 102, false, "2S8M", 30);
        final ReadEndsForMarkDuplicates ends = markDuplicates(input).buildReadEnds(input.header, 17, rec, false);

        Assert.assertEquals(ends.read1ReferenceIndex, CHR1);
        Assert.assertEquals(ends.read1Coordinate, 100);
        Assert.assertEquals(ends.orientation, ReadEnds.F);
        Assert.assertEquals(ends.read1IndexInFile, 17);
        Assert.assertEquals(ends.score, 300);
        Assert.assertEquals(ends.read2ReferenceIndex, -1);
        Assert.assertFalse(ends.isPaired());
        Assert.assertEquals(ends.getTile(), 1101);
        Assert.assertEquals(ends.getX(), 30);
        Assert.assertEquals(ends.getY(), 40);
    }

    @Test
    public void testReverseFragmentUsesTheUnclippedEnd() {
        final MarkDuplicatesEndToEndRunner input = header();
        final SAMRecord rec = input.fragment("RUN:1:7:8:9", "rg3", CHR1, 100, true, "8M2S3H", 20);
        final ReadEndsForMarkDuplicates ends = markDuplicates(input).buildReadEnds(input.header, 0, rec, false);

        Assert.assertEquals(ends.read1Coordinate, 112);
        Assert.assertEquals(ends.orientation, ReadEnds.R);
        Assert.assertEquals(ends.score, 200);
        Assert.assertEquals(ends.getTile(), 7);
    }

    @Test
    public void testPairedReadsRecordTheMateContig() {
        final MarkDuplicatesEndToEndRunner input = header();
        final SAMRecord[] pair = input.pair("p", "rg1", CHR1, 100, false, "10M", CHR3, 500, true, "10M", 30);
        final SAMRecord[] mateUnmapped = input.pairWithUnmappedMate("u", "rg1", CHR1, 100, false, "10M", 30);
        final MarkDuplicates markDuplicates = markDuplicates(input);

        final ReadEndsForMarkDuplicates first = markDuplicates.buildReadEnds(input.header, 0, pair[0], false);
        Assert.assertEquals(first.read2ReferenceIndex, CHR3);
        Assert.assertTrue(first.isPaired());
        final ReadEndsForMarkDuplicates second = markDuplicates.buildReadEnds(input.header, 1, pair[1], false);
        Assert.assertEquals(second.read1ReferenceIndex, CHR3);
        Assert.assertEquals(second.read1Coordinate, 509);
        Assert.assertEquals(second.read2ReferenceIndex, CHR1);
        final ReadEndsForMarkDuplicates single = markDuplicates.buildReadEnds(input.header, 2, mateUnmapped[0], false);
        Assert.assertEquals(single.read2ReferenceIndex, -1);
        Assert.assertFalse(single.isPaired());
    }

    /** Read groups of one library share a library id, which differs from that of another library. */
    @Test
    public void testLibraryIdsFollowTheLibrary() {
        final MarkDuplicatesEndToEndRunner input = header();
        final MarkDuplicates markDuplicates = markDuplicates(input);
        final SAMRecord libB = input.fragment("b", "rg2", CHR1, 100, false, "10M", 30);
        final SAMRecord libA = input.fragment("a", "rg3", CHR1, 100, false, "10M", 30);
        final SAMRecord libAOtherReadGroup = input.fragment("c", "rg1", CHR1, 100, false, "10M", 30);
        final short libBId = markDuplicates.buildReadEnds(input.header, 0, libB, false).libraryId;
        final short libAId = markDuplicates.buildReadEnds(input.header, 1, libA, false).libraryId;
        Assert.assertNotEquals(libAId, libBId);
        Assert.assertEquals(markDuplicates.buildReadEnds(input.header, 2, libAOtherReadGroup, false).libraryId, libAId);
        Assert.assertEquals(markDuplicates.buildReadEnds(input.header, 3, libB, false).libraryId, libBId);
    }

    @Test
    public void testNamesWithoutLocationOrNullRegex() {
        final MarkDuplicatesEndToEndRunner input = header();
        final SAMRecord unparsable = input.fragment("noLocation", "rg2", CHR1, 100, false, "10M", 30);
        final ReadEndsForMarkDuplicates ends = markDuplicates(input).buildReadEnds(input.header, 0, unparsable, false);
        Assert.assertFalse(ends.hasLocation());

        final MarkDuplicates nullRegex = new MarkDuplicates();
        nullRegex.READ_NAME_REGEX = null;
        nullRegex.setupOpticalDuplicateFinder();
        nullRegex.libraryIdGenerator = new LibraryIdGenerator(input.header);
        final SAMRecord parsable = input.fragment("RUN:1:1101:30:40", "rg2", CHR1, 100, false, "10M", 30);
        final ReadEndsForMarkDuplicates noLocation = nullRegex.buildReadEnds(input.header, 0, parsable, false);
        Assert.assertFalse(noLocation.hasLocation());
    }

    @DataProvider(name = "scores")
    public Object[][] scores() {
        return new Object[][]{
                {ScoringStrategy.SUM_OF_BASE_QUALITIES, "5M10D5M", 30, false, 300},
                {ScoringStrategy.SUM_OF_BASE_QUALITIES, "10M", 14, false, 0},
                {ScoringStrategy.SUM_OF_BASE_QUALITIES, "10M", 15, false, 150},
                {ScoringStrategy.SUM_OF_BASE_QUALITIES, "10M", 30, true, 300 - 16384},
                {ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, "5M10D5M", 30, false, 20},
                {ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, "3S5M2I5N", 30, false, 10},
                {ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, "10M", 30, true, 10 - 16384},
        };
    }

    @Test(dataProvider = "scores")
    public void testScores(final ScoringStrategy strategy, final String cigar, final int quality, final boolean failsQc, final int expected) {
        final MarkDuplicatesEndToEndRunner input = header();
        final SAMRecord rec = input.fragment("r", "rg1", CHR1, 100, false, cigar, quality);
        rec.setReadFailsVendorQualityCheckFlag(failsQc);
        final MarkDuplicates markDuplicates = markDuplicates(input);
        markDuplicates.DUPLICATE_SCORING_STRATEGY = strategy;
        Assert.assertEquals(markDuplicates.getReadDuplicateScore(rec, new ReadEndsForMarkDuplicates()), (short) expected);
        Assert.assertEquals(markDuplicates.buildReadEnds(input.header, 0, rec, false).score, (short) expected);
    }

    @Test
    public void testDefaultScoringStrategy() {
        Assert.assertEquals(new MarkDuplicates().DUPLICATE_SCORING_STRATEGY, ScoringStrategy.SUM_OF_BASE_QUALITIES);
    }

    private static ReadEndsForMarkDuplicates base() {
        final ReadEndsForMarkDuplicates ends = new ReadEndsForMarkDuplicates();
        ends.libraryId = 2;
        ends.read1ReferenceIndex = 1;
        ends.read1Coordinate = 500;
        ends.orientation = ReadEnds.FR;
        ends.read2ReferenceIndex = 1;
        ends.read2Coordinate = 700;
        ends.setTile((short) 5);
        ends.setX(50);
        ends.setY(60);
        ends.read1IndexInFile = 10;
        ends.read2IndexInFile = 20;
        ends.score = 7;
        ends.readGroup = 3;
        return ends;
    }

    /** One setter per sort key, in the order of precedence of the comparator. */
    @SuppressWarnings("unchecked")
    private static final Consumer<ReadEndsForMarkDuplicates>[] INCREASE_SORT_KEY = new Consumer[]{
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.libraryId++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read1ReferenceIndex++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read1Coordinate++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.orientation++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read2ReferenceIndex++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read2Coordinate++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.setTile((short) (e.getTile() + 1)),
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.setX(e.getX() + 1),
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.setY(e.getY() + 1),
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read1IndexInFile++,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read2IndexInFile++,
    };

    @SuppressWarnings("unchecked")
    private static final Consumer<ReadEndsForMarkDuplicates>[] DECREASE_SORT_KEY = new Consumer[]{
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.libraryId -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read1ReferenceIndex -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read1Coordinate -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.orientation -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read2ReferenceIndex -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read2Coordinate -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.setTile((short) (e.getTile() - 2)),
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.setX(e.getX() - 2),
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.setY(e.getY() - 2),
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read1IndexInFile -= 2,
            (Consumer<ReadEndsForMarkDuplicates>) e -> e.read2IndexInFile -= 2,
    };

    @DataProvider(name = "sortKeys")
    public Object[][] sortKeys() {
        final Object[][] rows = new Object[INCREASE_SORT_KEY.length][];
        for (int i = 0; i < rows.length; i++) {
            rows[i] = new Object[]{i};
        }
        return rows;
    }

    /** The sort keys that order read ends within a duplicate set: tile, x, y and the index in the file of read 1. */
    @DataProvider(name = "keysWithinADuplicateSet")
    public Object[][] keysWithinADuplicateSet() {
        return new Object[][]{{6}, {7}, {8}, {9}};
    }

    /** Within a duplicate set, a larger value of a sort key orders later whatever the values of the keys after it. */
    @Test(dataProvider = "keysWithinADuplicateSet")
    public void testComparatorPrecedence(final int key) {
        final MarkDuplicates.ReadEndsMDComparator comparator = new MarkDuplicates.ReadEndsMDComparator(false);
        final ReadEndsForMarkDuplicates lower = base();
        final ReadEndsForMarkDuplicates higher = base();
        INCREASE_SORT_KEY[key].accept(higher);
        for (int later = key + 1; later < DECREASE_SORT_KEY.length; later++) {
            DECREASE_SORT_KEY[later].accept(higher);
        }
        Assert.assertTrue(comparator.compare(lower, higher) < 0);
        Assert.assertTrue(comparator.compare(higher, lower) > 0);
    }

    /**
     * Sorted read ends of one duplicate set, of pairs (compared on both ends) or of fragments (compared on read 1
     * only), are next to each other whatever their tile, x, y and indices in the file.
     */
    @Test
    public void testComparatorKeepsEachDuplicateSetContiguous() {
        final Random random = new Random(42);
        final List<ReadEndsForMarkDuplicates> ends = new ArrayList<>();
        // two values of each of the six keys compared for pairs, three read ends with each combination
        for (int combination = 0; combination < 64; combination++) {
            for (int copy = 0; copy < 3; copy++) {
                final ReadEndsForMarkDuplicates e = base();
                for (int key = 0; key < 6; key++) {
                    if ((combination >> key & 1) == 1) {
                        INCREASE_SORT_KEY[key].accept(e);
                    }
                }
                e.setTile((short) random.nextInt(3));
                e.setX(random.nextInt(3));
                e.setY(random.nextInt(3));
                e.read1IndexInFile = random.nextInt(1000);
                e.read2IndexInFile = random.nextInt(1000);
                ends.add(e);
            }
        }
        Collections.shuffle(ends, random);
        ends.sort(new MarkDuplicates.ReadEndsMDComparator(false));

        final MarkDuplicates markDuplicates = new MarkDuplicates();
        for (final boolean compareRead2 : new boolean[]{true, false}) {
            final List<ReadEndsForMarkDuplicates> earlierSets = new ArrayList<>();
            ReadEndsForMarkDuplicates firstOfSet = ends.get(0);
            for (final ReadEndsForMarkDuplicates e : ends) {
                if (!markDuplicates.areComparableForDuplicates(firstOfSet, e, compareRead2, false)) {
                    earlierSets.add(firstOfSet);
                    firstOfSet = e;
                }
                for (final ReadEndsForMarkDuplicates earlier : earlierSets) {
                    Assert.assertFalse(markDuplicates.areComparableForDuplicates(earlier, e, compareRead2, false));
                }
            }
            Assert.assertEquals(earlierSets.size() + 1, compareRead2 ? 64 : 16);
        }
    }

    /**
     * Score, read group and optical orientation, higher or lower, do not reorder read ends of one duplicate set, whose
     * order decides which read of an optical cluster without the keeper is flagged: the earlier read in the file
     * still sorts first.
     */
    @Test
    public void testComparatorIgnoresScoreAndReadGroup() {
        final MarkDuplicates.ReadEndsMDComparator comparator = new MarkDuplicates.ReadEndsMDComparator(false);
        for (final int sign : new int[]{-1, 1}) {
            final ReadEndsForMarkDuplicates first = base();
            first.orientationForOpticalDuplicates = ReadEnds.FR;
            final ReadEndsForMarkDuplicates second = base();
            second.read1IndexInFile++;
            second.score += sign * 50;
            second.readGroup += sign;
            second.orientationForOpticalDuplicates = sign < 0 ? ReadEnds.F : ReadEnds.RF;
            Assert.assertTrue(comparator.compare(first, second) < 0);
            Assert.assertTrue(comparator.compare(second, first) > 0);
        }
    }

    /** Only the library, the 5' ends and the orientation decide whether read ends are in one duplicate set. */
    @Test(dataProvider = "sortKeys")
    public void testAreComparableForDuplicates(final int key) {
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        final ReadEndsForMarkDuplicates other = base();
        INCREASE_SORT_KEY[key].accept(other);
        // keys 0 to 3 are compared for fragments and pairs, 4 and 5 for pairs only, the rest never
        Assert.assertEquals(markDuplicates.areComparableForDuplicates(base(), other, true, false), key > 5);
        Assert.assertEquals(markDuplicates.areComparableForDuplicates(base(), other, false, false), key > 3);
        Assert.assertEquals(markDuplicates.areComparableForDuplicates(other, base(), true, false), key > 5);
    }

    @Test
    public void testAreComparableIgnoresScoreAndReadGroup() {
        final ReadEndsForMarkDuplicates other = base();
        other.score = 100;
        other.readGroup = 0;
        Assert.assertTrue(new MarkDuplicates().areComparableForDuplicates(base(), other, true, false));
    }

    @Test
    public void testDuplicateTypeCodes() {
        Assert.assertEquals(MarkDuplicates.DuplicateType.LIBRARY.code(), "LB");
        Assert.assertEquals(MarkDuplicates.DuplicateType.SEQUENCING.code(), "SQ");
        Assert.assertEquals(MarkDuplicates.DUPLICATE_TYPE_LIBRARY, "LB");
        Assert.assertEquals(MarkDuplicates.DUPLICATE_TYPE_SEQUENCING, "SQ");
        Assert.assertEquals(MarkDuplicates.DUPLICATE_TYPE_TAG, "DT");
    }
}
