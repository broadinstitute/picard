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

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SamRecordTrackingBuffer;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Tests for {@link MarkQueue}, driven the way MarkDuplicatesWithMateCigarIterator drives it: each record goes into a
 * {@link SamRecordTrackingBuffer} and then into the queue, and the duplicate flags are read from the records the
 * buffer emits. Covers keeping the best end among comparable ends, ties, the two ends of one pair at one location,
 * duplicate pairs in either order, fragments at a pair's 5' position, ends that are not comparable, the poll order,
 * the duplicate counts and metrics, the optical-duplicate location sets, and the comparator that picks the
 * representative end.
 */
public class MarkQueueTest {

    private static final class Fixture implements AutoCloseable {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final SAMFileHeader header = records.header;
        final LibraryIdGenerator libraryIdGenerator = new LibraryIdGenerator(header);
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder();
        final MarkQueue queue = new MarkQueue(ScoringStrategy.SUM_OF_BASE_QUALITIES);
        final File tmpDir = IOUtil.createTempDir("markQueue").toFile();
        final SamRecordTrackingBuffer<SamRecordWithOrdinalAndSetDuplicateReadFlag> buffer = new SamRecordTrackingBuffer<>(
                1000, 100, Collections.singletonList(tmpDir), header, SamRecordWithOrdinalAndSetDuplicateReadFlag.class);
        private long ordinal = 0;

        Fixture() {
            final short libA = libraryIdGenerator.getLibraryId(records.fragment("x", "rgA", 0, 1, false, "10M", 30));
            final short libB = libraryIdGenerator.getLibraryId(records.fragment("x", "rgB", 0, 1, false, "10M", 30));
            Assert.assertNotEquals(libA, libB);
        }

        DuplicationMetrics metrics(final String library) {
            return libraryIdGenerator.getMetricsByLibrary(library);
        }

        ReadEndsForMateCigar add(final SAMRecord rec) {
            final SamRecordWithOrdinalAndSetDuplicateReadFlag wrapped = new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, ordinal++);
            buffer.add(wrapped);
            final ReadEndsForMateCigar ends = new ReadEndsForMateCigar(header, wrapped, finder, libraryIdGenerator.getLibraryId(rec));
            queue.add(ends, buffer, metrics(LibraryIdGenerator.getLibraryName(header, rec)));
            return ends;
        }

        ReadEndsForMateCigar poll() {
            return queue.poll(buffer, header, finder, libraryIdGenerator);
        }

        /** Polls until the queue is empty and returns the names of the polled ends in order. */
        List<String> pollAll() {
            final List<String> names = new ArrayList<>();
            while (!queue.isEmpty()) {
                names.add(key(poll().getRecord()));
            }
            Assert.assertEquals(queue.size(), 0);
            return names;
        }

        /** Emits every record whose state is known and returns its duplicate flag by name. */
        Map<String, Boolean> emit() {
            final Map<String, Boolean> flags = new HashMap<>();
            while (!buffer.isEmpty() && buffer.canEmit()) {
                final SAMRecord rec = buffer.next().getRecord();
                flags.put(key(rec), rec.getDuplicateReadFlag());
            }
            return flags;
        }

        @Override
        public void close() {
            buffer.close();
            IOUtil.deleteDirectoryTree(tmpDir);
        }
    }

    private static String key(final SAMRecord rec) {
        final String name = rec.getReadName().split(":")[0];
        if (!rec.getReadPairedFlag()) return name;
        return name + (rec.getFirstOfPairFlag() ? "/1" : "/2");
    }

    /** A read name whose last three fields give tile 1101 and x = y = {@code xy}. */
    private static String name(final String prefix, final int xy) {
        return prefix + ":1:FC:1:1101:" + xy + ":" + xy;
    }

    private static Map<String, Boolean> flags(final Object... keysAndFlags) {
        final Map<String, Boolean> flags = new HashMap<>();
        for (int i = 0; i < keysAndFlags.length; i += 2) {
            flags.put((String) keysAndFlags[i], (Boolean) keysAndFlags[i + 1]);
        }
        return flags;
    }

    private static void assertMetrics(final DuplicationMetrics metrics, final long unpairedDuplicates, final long pairDuplicates) {
        Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, unpairedDuplicates, "UNPAIRED_READ_DUPLICATES");
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES, pairDuplicates, "READ_PAIR_DUPLICATES");
    }

    @Test
    public void testNewQueue() {
        final MarkQueue queue = new MarkQueue(ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
        Assert.assertEquals(queue.size(), 0);
        Assert.assertTrue(queue.isEmpty());
        Assert.assertEquals(queue.getNumDuplicates(), 0);
        Assert.assertEquals(queue.getToMarkQueueMinimumDistance(), -1);
        queue.setToMarkQueueMinimumDistance(250);
        Assert.assertEquals(queue.getToMarkQueueMinimumDistance(), 250);
    }

    /** Sums of base qualities: a = 300, b = 200, c = 400; c has a 3 base soft clip that puts its 5' end at 100. */
    @Test
    public void testBestFragmentIsKept() {
        try (final Fixture f = new Fixture()) {
            f.add(f.records.fragment(name("a", 1), "rgA", 0, 100, false, "10M", 30));
            Assert.assertEquals(f.queue.size(), 1);
            Assert.assertEquals(f.queue.getNumDuplicates(), 0);

            f.add(f.records.fragment(name("b", 2), "rgA", 0, 100, false, "10M", 20));
            Assert.assertEquals(f.queue.size(), 1);
            Assert.assertEquals(f.queue.getNumDuplicates(), 1);
            assertMetrics(f.metrics("libA"), 1, 0);

            final ReadEndsForMateCigar c = f.add(f.records.fragment(name("c", 3), "rgA", 0, 103, false, "3S7M", 40));
            Assert.assertEquals(f.queue.size(), 1);
            Assert.assertEquals(f.queue.getNumDuplicates(), 2);
            Assert.assertSame(f.queue.peek(), c);

            Assert.assertSame(f.poll(), c);
            Assert.assertTrue(f.queue.isEmpty());
            Assert.assertEquals(f.emit(), flags("a", true, "b", true, "c", false));
            Assert.assertEquals(f.queue.getNumDuplicates(), 2);
            assertMetrics(f.metrics("libA"), 2, 0);
            assertMetrics(f.metrics("libB"), 0, 0);
        }
    }

    @DataProvider(name = "orders")
    public Object[][] orders() {
        return new Object[][]{{true}, {false}};
    }

    /** With equal scores the end with the lexicographically smaller name is kept. */
    @Test(dataProvider = "orders")
    public void testTiedScoresKeepSmallerName(final boolean smallerFirst) {
        try (final Fixture f = new Fixture()) {
            final SAMRecord d = f.records.fragment(name("d", 1), "rgA", 0, 100, true, "10M", 30);
            final SAMRecord n = f.records.fragment(name("n", 2), "rgA", 0, 100, true, "10M", 30);
            f.add(smallerFirst ? d : n);
            f.add(smallerFirst ? n : d);
            Assert.assertEquals(f.queue.getNumDuplicates(), 1);
            Assert.assertEquals(f.pollAll(), List.of("d"));
            Assert.assertEquals(f.emit(), flags("d", false, "n", true));
        }
    }

    /** Reverse reads are compared at their unclipped 3' alignment end. */
    @Test
    public void testReverseFragmentsUseUnclippedEnd() {
        try (final Fixture f = new Fixture()) {
            f.add(f.records.fragment(name("c", 3), "rgA", 0, 191, false, "10M", 40)); // forward, 5' end 191
            f.add(f.records.fragment(name("a", 1), "rgA", 0, 191, true, "10M", 20));  // reverse, 5' end 200
            f.add(f.records.fragment(name("b", 2), "rgA", 0, 191, true, "7M3S", 30)); // reverse, 5' end 200
            f.add(f.records.fragment(name("e", 4), "rgA", 0, 192, true, "10M", 10));  // reverse, 5' end 201
            Assert.assertEquals(f.queue.getNumDuplicates(), 1);
            Assert.assertEquals(f.queue.size(), 3);
            Assert.assertEquals(f.pollAll(), List.of("c", "b", "e"));
            Assert.assertEquals(f.emit(), flags("a", true, "b", false, "c", false, "e", false));
            assertMetrics(f.metrics("libA"), 1, 0);
        }
    }

    /**
     * Each row is one end of a pair at chr1:100-109 and a fragment of higher quality. The fragment is a duplicate when
     * its strand and 5' end match the pair end's, and it is found when the pair end is polled.
     */
    @DataProvider(name = "fragmentAtPairPosition")
    public Object[][] fragmentAtPairPosition() {
        return new Object[][]{
                // pair end strand, mate strand, fragment strand, fragment start, fragment is a duplicate
                {false, true, false, 100, true},   // FR, forward fragment at 100
                {false, false, false, 100, true},  // FF
                {true, false, true, 100, true},    // RF, reverse fragment ending at 109
                {true, true, true, 100, true},     // RR
                {false, true, true, 91, false},    // FR at 100, reverse fragment ending at 100
                {true, false, false, 109, false},  // RF ending at 109, forward fragment starting at 109
                {false, true, false, 101, false}   // FR at 100, forward fragment at 101
        };
    }

    @Test(dataProvider = "fragmentAtPairPosition")
    public void testFragmentAtPairPosition(final boolean negativeStrand, final boolean mateNegativeStrand,
                                           final boolean fragmentNegativeStrand, final int fragmentStart, final boolean duplicate) {
        try (final Fixture f = new Fixture()) {
            final ReadEndsForMateCigar pairEnd = f.add(f.records.pairEnd(name("p", 1), "rgA", 0, 100, negativeStrand, "10M", true,
                    0, 1000, mateNegativeStrand, "10M", 20));
            f.add(f.records.fragment(name("f", 2), "rgA", 0, fragmentStart, fragmentNegativeStrand, "10M", 40));
            Assert.assertEquals(f.queue.size(), 2);
            Assert.assertEquals(f.queue.getNumDuplicates(), 0);

            Assert.assertSame(f.poll(), pairEnd);
            Assert.assertEquals(f.queue.getNumDuplicates(), duplicate ? 1 : 0);
            Assert.assertEquals(f.queue.size(), duplicate ? 0 : 1);
            Assert.assertEquals(f.pollAll(), duplicate ? List.of() : List.of("f"));
            Assert.assertEquals(f.emit(), flags("p/1", false, "f", duplicate));
            assertMetrics(f.metrics("libA"), duplicate ? 1 : 0, 0);
        }
    }

    /** The two ends of one pair at the same 5' position and strand are not duplicates of each other. */
    @Test(dataProvider = "orders")
    public void testBothEndsOfOnePairAtOneLocation(final boolean firstEndFirst) {
        try (final Fixture f = new Fixture()) {
            final SAMRecord[] pair = f.records.pair(name("p", 5), "rgA", 0, 100, false, "10M", 100, false, "10M", 30);
            final ReadEndsForMateCigar end1;
            if (firstEndFirst) {
                end1 = f.add(pair[0]);
                f.add(pair[1]);
            } else {
                f.add(pair[1]);
                end1 = f.add(pair[0]);
            }
            Assert.assertEquals(f.queue.size(), 2);
            Assert.assertEquals(f.queue.getNumDuplicates(), 0);
            Assert.assertSame(f.queue.peek(), end1);

            Assert.assertSame(f.poll(), end1);
            Assert.assertTrue(f.queue.isEmpty());
            Assert.assertEquals(f.queue.size(), 0);
            Assert.assertEquals(f.emit(), flags("p/1", false, "p/2", false));

            Assert.assertTrue(f.queue.shouldBeInLocations(end1));
            final Set<ReadEnds> locations = f.queue.getLocations(end1);
            Assert.assertEquals(locations.size(), 1);
            Assert.assertTrue(locations.contains(end1));
            assertMetrics(f.metrics("libA"), 0, 0);
        }
    }

    @DataProvider(name = "duplicatePairs")
    public Object[][] duplicatePairs() {
        return new Object[][]{
                {true, "rgA", "libA", "libB"},
                {false, "rgA", "libA", "libB"},
                {true, "rgB", "libB", "libA"},
                {false, "rgB", "libB", "libA"}
        };
    }

    /**
     * Pairs q (sum of qualities 300 per end) and r (200 per end), all four ends forward at chr1:100. Both ends of r
     * are duplicates whichever pair comes first. The location set of the kept end holds one end per pair.
     */
    @Test(dataProvider = "duplicatePairs")
    public void testDuplicatePairsAtOneLocation(final boolean betterFirst, final String readGroup, final String library,
                                                final String otherLibrary) {
        try (final Fixture f = new Fixture()) {
            final SAMRecord[] q = f.records.pair(name("q", 100), readGroup, 0, 100, false, "10M", 100, false, "10M", 30);
            final SAMRecord[] r = f.records.pair(name("r", 900), readGroup, 0, 100, false, "10M", 100, false, "10M", 20);
            final ReadEndsForMateCigar q1;
            final ReadEndsForMateCigar r1;
            if (betterFirst) {
                q1 = f.add(q[0]);
                f.add(q[1]);
                r1 = f.add(r[0]);
                f.add(r[1]);
            } else {
                r1 = f.add(r[0]);
                f.add(r[1]);
                q1 = f.add(q[0]);
                f.add(q[1]);
            }
            Assert.assertEquals(f.queue.getNumDuplicates(), 2);
            Assert.assertEquals(f.queue.size(), 2);
            assertMetrics(f.metrics(library), 0, 2);

            Assert.assertSame(f.poll(), q1);
            Assert.assertTrue(f.queue.isEmpty());
            Assert.assertEquals(f.emit(), flags("q/1", false, "q/2", false, "r/1", true, "r/2", true));

            final Set<ReadEnds> locations = f.queue.getLocations(q1);
            Assert.assertEquals(locations.size(), 2);
            Assert.assertTrue(locations.contains(q1));
            Assert.assertTrue(locations.contains(r1));
            final PicardException e = Assert.expectThrows(PicardException.class, () -> f.queue.getLocations(q1));
            Assert.assertTrue(e.getMessage().contains("Already called"), e.getMessage());

            assertMetrics(f.metrics(library), 0, 2);
            assertMetrics(f.metrics(otherLibrary), 0, 0);
        }
    }

    /** Pairs of FR ends at the same 5' position whose mates map elsewhere are not comparable. */
    @Test
    public void testPairsWithDifferentMatesOrLibrariesAreNotDuplicates() {
        try (final Fixture f = new Fixture()) {
            f.add(f.records.pairEnd(name("s", 3), "rgB", 0, 100, false, "10M", true, 0, 300, true, "10M", 40));
            f.add(f.records.pairEnd(name("r", 2), "rgA", 0, 100, false, "10M", true, 0, 310, true, "10M", 30));
            f.add(f.records.pairEnd(name("q", 1), "rgA", 0, 100, false, "10M", true, 0, 300, true, "10M", 20));
            f.add(f.records.pairEnd(name("t", 4), "rgA", 0, 100, false, "10M", true, 1, 300, true, "10M", 30));
            Assert.assertEquals(f.queue.size(), 4);
            Assert.assertEquals(f.queue.getNumDuplicates(), 0);
            // mates at chr1:309, chr1:319 and chr2:309; library 2 last
            Assert.assertEquals(f.pollAll(), List.of("q/1", "r/1", "t/1", "s/1"));
            Assert.assertEquals(f.emit(), flags("q/1", false, "r/1", false, "s/1", false, "t/1", false));
        }
    }

    /** Reads whose mates are unmapped are compared like fragments, by score, and count as unpaired duplicates. */
    @Test(dataProvider = "orders")
    public void testMateUnmappedReads(final boolean betterFirst) {
        try (final Fixture f = new Fixture()) {
            final SAMRecord better = f.records.mateUnmapped(name("h", 1), "rgA", 0, 100, false, "10M", true, 30);
            final SAMRecord worse = f.records.mateUnmapped(name("g", 2), "rgA", 0, 100, false, "10M", false, 20);
            f.add(betterFirst ? better : worse);
            f.add(betterFirst ? worse : better);
            Assert.assertEquals(f.queue.getNumDuplicates(), 1);
            Assert.assertEquals(f.pollAll(), List.of("h/1"));
            Assert.assertEquals(f.emit(), flags("h/1", false, "g/2", true));
            assertMetrics(f.metrics("libA"), 1, 0);
        }
    }

    /**
     * Pair ends at one position are polled in the order RF, RR, FR, FF; polling an R-type end removes the reverse
     * fragment at that position and polling an F-type end the forward one.
     */
    @Test
    public void testPollOrderByOrientation() {
        try (final Fixture f = new Fixture()) {
            // each pair end has its 5' end at chr1:200
            f.add(f.records.pairEnd(name("ff", 1), "rgA", 0, 200, false, "10M", true, 0, 2000, false, "10M", 30));
            f.add(f.records.fragment(name("fwd", 5), "rgA", 0, 200, false, "10M", 40));
            f.add(f.records.pairEnd(name("fr", 2), "rgA", 0, 200, false, "10M", true, 0, 3000, true, "10M", 30));
            f.add(f.records.pairEnd(name("rr", 3), "rgA", 0, 191, true, "10M", true, 0, 4000, true, "10M", 30));
            f.add(f.records.fragment(name("rev", 6), "rgA", 0, 191, true, "10M", 40));
            f.add(f.records.pairEnd(name("rf", 4), "rgA", 0, 191, true, "10M", true, 0, 5000, false, "10M", 30));
            Assert.assertEquals(f.queue.size(), 6);
            Assert.assertEquals(f.queue.getNumDuplicates(), 0);

            final List<String> polled = new ArrayList<>();
            polled.add(key(f.poll().getRecord()));
            Assert.assertEquals(f.queue.getNumDuplicates(), 1);
            Assert.assertEquals(f.queue.size(), 4);
            polled.addAll(f.pollAll());
            Assert.assertEquals(polled, List.of("rf/1", "rr/1", "fr/1", "ff/1"));
            Assert.assertEquals(f.queue.getNumDuplicates(), 2);
            Assert.assertEquals(f.emit(), flags("ff/1", false, "fwd", true, "fr/1", false, "rr/1", false, "rev", true, "rf/1", false));
            assertMetrics(f.metrics("libA"), 2, 0);
        }
    }

    /** The queue is ordered by library id, reference, 5' position, orientation (pairs first) and the mate's position. */
    @Test
    public void testPollOrder() {
        try (final Fixture f = new Fixture()) {
            f.add(f.records.fragment(name("lb", 1), "rgB", 0, 100, false, "10M", 30));
            f.add(f.records.fragment(name("c2", 2), "rgA", 1, 50, false, "10M", 30));
            f.add(f.records.fragment(name("f300", 3), "rgA", 0, 300, false, "10M", 30));
            f.add(f.records.pairEnd(name("m900", 4), "rgA", 0, 150, false, "10M", true, 0, 891, true, "10M", 30));
            f.add(f.records.fragment(name("f100", 5), "rgA", 0, 100, false, "10M", 30));
            f.add(f.records.pairEnd(name("mchr2", 6), "rgA", 0, 150, false, "10M", true, 1, 1, true, "10M", 30));
            f.add(f.records.pairEnd(name("m800", 7), "rgA", 0, 150, false, "10M", true, 0, 791, true, "10M", 30));
            f.add(f.records.fragment(name("r100", 8), "rgC", 0, 91, true, "10M", 30));
            Assert.assertEquals(f.queue.size(), 8);
            Assert.assertEquals(f.queue.peek().getRecordReadName(), name("r100", 8));
            Assert.assertEquals(f.pollAll(), List.of("r100", "f100", "m800/1", "m900/1", "mchr2/1", "f300", "c2", "lb"));
            Assert.assertEquals(f.queue.getNumDuplicates(), 0);
            Assert.assertEquals(f.emit().size(), 8);
        }
    }

    @Test
    public void testShouldBeInLocations() {
        try (final Fixture f = new Fixture()) {
            final ReadEndsForMateCigar fragment = f.add(f.records.fragment(name("a", 1), "rgA", 0, 100, false, "10M", 30));
            final ReadEndsForMateCigar halfMapped = f.add(f.records.mateUnmapped(name("b", 2), "rgA", 0, 200, false, "10M", true, 30));
            final ReadEndsForMateCigar pairEnd = f.add(f.records.pairEnd(name("c", 3), "rgA", 0, 300, false, "10M", false, 0, 400, true, "10M", 30));
            Assert.assertFalse(f.queue.shouldBeInLocations(fragment));
            Assert.assertFalse(f.queue.shouldBeInLocations(halfMapped));
            Assert.assertTrue(f.queue.shouldBeInLocations(pairEnd));

            Assert.assertNull(fragment.getLocationSet());
            Assert.assertNull(halfMapped.getLocationSet());
            Assert.assertEquals(pairEnd.getLocationSet().size(), 1);

            for (final ReadEndsForMateCigar ends : List.of(fragment, halfMapped)) {
                final PicardException e = Assert.expectThrows(PicardException.class, () -> f.queue.getLocations(ends));
                Assert.assertTrue(e.getMessage().contains("Not implemented"), e.getMessage());
            }
            final Set<ReadEnds> locations = f.queue.getLocations(pairEnd);
            Assert.assertEquals(locations.size(), 1);
            Assert.assertTrue(locations.contains(pairEnd));
        }
    }

    /** Three pairs at one location, the best in the middle; every end of the other two pairs is a duplicate. */
    @Test
    public void testThreePairsAtOneLocation() {
        try (final Fixture f = new Fixture()) {
            final SAMRecord[] a = f.records.pair(name("a", 10), "rgA", 0, 100, true, "10M", 100, true, "10M", 20);
            final SAMRecord[] b = f.records.pair(name("b", 20), "rgA", 0, 100, true, "10M", 100, true, "10M", 40);
            final SAMRecord[] c = f.records.pair(name("c", 30), "rgA", 0, 100, true, "10M", 100, true, "10M", 30);
            f.add(a[0]);
            f.add(b[1]);
            f.add(c[0]);
            f.add(a[1]);
            final ReadEndsForMateCigar b1 = f.add(b[0]);
            f.add(c[1]);
            Assert.assertEquals(f.queue.getNumDuplicates(), 4);
            assertMetrics(f.metrics("libA"), 0, 4);
            Assert.assertSame(f.queue.peek(), b1);
            Assert.assertEquals(f.pollAll(), List.of("b/1"));
            Assert.assertEquals(f.emit(), flags("a/1", true, "a/2", true, "b/1", false, "b/2", false, "c/1", true, "c/2", true));
            Assert.assertEquals(f.queue.getLocations(b1).size(), 3);
        }
    }

    private static MarkQueue.ReadEndsMCComparator comparator() {
        final MarkQueue queue = new MarkQueue(ScoringStrategy.SUM_OF_BASE_QUALITIES);
        return queue.new ReadEndsMCComparator(ScoringStrategy.SUM_OF_BASE_QUALITIES);
    }

    private static ReadEndsForMateCigar ends(final MateCigarTestRecords records, final SAMRecord rec) {
        return new ReadEndsForMateCigar(records.header, new SamRecordWithOrdinalAndSetDuplicateReadFlag(rec, 0),
                new OpticalDuplicateFinder(), (short) 1);
    }

    /** Asserts that {@code first} sorts strictly before {@code second}, in both argument orders. */
    private static void assertBefore(final MarkQueue.ReadEndsMCComparator comparator, final ReadEndsForMateCigar first, final ReadEndsForMateCigar second) {
        Assert.assertTrue(comparator.compare(first, second) < 0, "compare(first, second) should be negative");
        Assert.assertTrue(comparator.compare(second, first) > 0, "compare(second, first) should be positive");
    }

    @Test
    public void testComparatorPositionFields() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final MarkQueue.ReadEndsMCComparator comparator = comparator();
        final SAMRecord rec = records.pairEnd("p", "rgA", 0, 100, false, "10M", true, 0, 300, true, "10M", 30);

        final ReadEndsForMateCigar base = ends(records, rec);
        final ReadEndsForMateCigar other = ends(records, rec);
        other.libraryId = 2;
        assertBefore(comparator, base, other);

        final ReadEndsForMateCigar laterReference = ends(records, rec);
        laterReference.read1ReferenceIndex = 1;
        laterReference.read1Coordinate = 1;
        assertBefore(comparator, base, laterReference);

        final ReadEndsForMateCigar laterCoordinate = ends(records, rec);
        laterCoordinate.read1Coordinate = 101;
        laterCoordinate.orientation = ReadEnds.RF;
        assertBefore(comparator, base, laterCoordinate);

        final ReadEndsForMateCigar lowerOrientation = ends(records, rec);
        lowerOrientation.orientation = ReadEnds.FF;
        assertBefore(comparator, base, lowerOrientation);

        final ReadEndsForMateCigar hasUnmapped = ends(records, rec);
        hasUnmapped.hasUnmapped = 1;
        hasUnmapped.read2ReferenceIndex = 0;
        hasUnmapped.read2Coordinate = 0;
        assertBefore(comparator, base, hasUnmapped);

        final ReadEndsForMateCigar laterMateReference = ends(records, rec);
        laterMateReference.read2ReferenceIndex = 1;
        laterMateReference.read2Coordinate = 0;
        assertBefore(comparator, base, laterMateReference);

        final ReadEndsForMateCigar laterMateCoordinate = ends(records, rec);
        laterMateCoordinate.read2Coordinate = 310;
        assertBefore(comparator, base, laterMateCoordinate);
    }

    @Test
    public void testComparatorOrientationOrder() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final MarkQueue.ReadEndsMCComparator comparator = comparator();
        final SAMRecord rec = records.pairEnd("p", "rgA", 0, 100, false, "10M", true, 0, 300, true, "10M", 30);
        final byte[] order = {ReadEnds.RF, ReadEnds.RR, ReadEnds.FR, ReadEnds.FF, ReadEnds.R, ReadEnds.F};
        for (int i = 0; i + 1 < order.length; i++) {
            final ReadEndsForMateCigar first = ends(records, rec);
            first.orientation = order[i];
            final ReadEndsForMateCigar second = ends(records, rec);
            second.orientation = order[i + 1];
            assertBefore(comparator, first, second);
        }
    }

    /** At the same position the higher sum of base qualities comes first, then the smaller read name. */
    @Test
    public void testComparatorScoreThenName() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final MarkQueue.ReadEndsMCComparator comparator = comparator();
        final ReadEndsForMateCigar high = ends(records, records.fragment("z", "rgA", 0, 100, false, "10M", 31));
        final ReadEndsForMateCigar low = ends(records, records.fragment("a", "rgA", 0, 100, false, "10M", 30));
        assertBefore(comparator, high, low);

        final ReadEndsForMateCigar a = ends(records, records.fragment("a", "rgA", 0, 100, false, "10M", 30));
        final ReadEndsForMateCigar b = ends(records, records.fragment("b", "rgA", 0, 100, false, "10M", 30));
        assertBefore(comparator, a, b);

        final ReadEndsForMateCigar same = ends(records, records.fragment("a", "rgA", 0, 100, false, "10M", 30));
        Assert.assertEquals(comparator.compare(a, same), 0);
    }

    /** Quality scores below 15 do not count toward the sum of base qualities. */
    @Test
    public void testComparatorIgnoresLowQualities() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final MarkQueue.ReadEndsMCComparator comparator = comparator();
        final ReadEndsForMateCigar lowQualities = ends(records, records.fragment("a", "rgA", 0, 100, false, "10M", 14));
        final ReadEndsForMateCigar threshold = ends(records, records.fragment("b", "rgA", 0, 100, false, "10M", 15));
        assertBefore(comparator, threshold, lowQualities);
        final ReadEndsForMateCigar zero = ends(records, records.fragment("c", "rgA", 0, 100, false, "10M", 2));
        // 14 and 2 both score zero, so the name decides
        assertBefore(comparator, lowQualities, zero);
    }

    /** When everything else ties for the two ends of one pair, the first end comes first. */
    @Test
    public void testComparatorPrefersFirstEndOfPair() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final MarkQueue.ReadEndsMCComparator comparator = comparator();
        final SAMRecord[] pair = records.pair("p", "rgA", 0, 100, false, "10M", 100, false, "10M", 30);
        final ReadEndsForMateCigar end1 = ends(records, pair[0]);
        final ReadEndsForMateCigar end2 = ends(records, pair[1]);
        assertBefore(comparator, end1, end2);
    }

    @Test
    public void testComparatorTotalMappedReferenceLength() {
        final MateCigarTestRecords records = new MateCigarTestRecords();
        final MarkQueue queue = new MarkQueue(ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
        final MarkQueue.ReadEndsMCComparator comparator = queue.new ReadEndsMCComparator(ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH);
        // same 5' end at 100; the longer alignment is kept even with lower base qualities
        final ReadEndsForMateCigar longer = ends(records, records.fragment("z", "rgA", 0, 100, false, "12M", 20));
        final ReadEndsForMateCigar shorter = ends(records, records.fragment("a", "rgA", 0, 100, false, "10M2S", 40));
        assertBefore(comparator, longer, shorter);
    }

    @Test
    public void testLibrariesAreKeptApartInAdd() {
        try (final Fixture f = new Fixture()) {
            f.add(f.records.fragment(name("a", 1), "rgA", 0, 100, false, "10M", 30));
            f.add(f.records.fragment(name("b", 2), "rgB", 0, 100, false, "10M", 20));
            f.add(f.records.fragment(name("c", 3), "rgC", 0, 100, false, "10M", 40));
            // rgC is in libA with rgA
            Assert.assertEquals(f.queue.getNumDuplicates(), 1);
            Assert.assertEquals(f.pollAll(), List.of("c", "b"));
            Assert.assertEquals(f.emit(), flags("a", true, "b", false, "c", false));
            assertMetrics(f.metrics("libA"), 1, 0);
            assertMetrics(f.metrics("libB"), 0, 0);
        }
    }
}
