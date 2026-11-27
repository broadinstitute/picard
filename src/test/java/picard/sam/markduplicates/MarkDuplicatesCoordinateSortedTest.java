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

import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.Result;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertMetrics;
import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertSameRecords;

/**
 * Tests which records {@link MarkDuplicates} flags as duplicates in coordinate-sorted input, and the per-library
 * metrics it reports, on inputs small enough to work out by hand: pairs, fragments, fragments at pair ends,
 * clipped reads, the four pair orientations, mates on other contigs, unmapped reads and mates, secondary and
 * supplementary records, missing mates, libraries and read groups, several inputs, and each duplicate scoring
 * strategy with its tie-breaking.
 *
 * Read names here have no tile and coordinates (except where a test says so), so no duplicate is optical. Every
 * test also checks that the output holds the input records in input order with only the duplicate flag and the PG
 * tag changed.
 */
public class MarkDuplicatesCoordinateSortedTest {
    private static final int CHR1 = 0, CHR2 = 1, CHR3 = 2;
    private static final boolean F = false, R = true;

    private static MarkDuplicatesEndToEndRunner oneLibrary() {
        return new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
    }

    private static Result runAndCheckRecords(final MarkDuplicatesEndToEndRunner input, final String... args) {
        final Result result = input.run(args);
        assertSameRecords(input.inputRecords(), result.records, SAMTag.PG.name());
        return result;
    }

    @Test
    public void testPairsAreDuplicatesOnlyWhenBothFivePrimeEndsMatch() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // 5' ends at chr1:100 forward and chr1:300 reverse
        input.pair("best", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("dup", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.pair("otherSecondEnd", "rg", CHR1, 100, F, "10M", CHR1, 292, R, "10M", 20);
        input.pair("otherContig", "rg", CHR2, 100, F, "10M", CHR2, 291, R, "10M", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("dup/1", "dup/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 4, 0, 0, 0, 1, 0);
        assertCoverageMultHistogram(result, 4, 1, 0);
    }

    @Test
    public void testFragmentsAreDuplicatesWhenFivePrimeEndAndStrandMatch() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.fragment("forwardBest", "rg", CHR1, 100, F, "10M", 30);
        input.fragment("forwardDup", "rg", CHR1, 100, F, "10M", 20);
        input.fragment("forwardOneBaseLater", "rg", CHR1, 101, F, "10M", 20);
        // reverse reads: the 5' end is the alignment end, chr1:109
        input.fragment("reverseBest", "rg", CHR1, 100, R, "10M", 30);
        input.fragment("reverseDupShorter", "rg", CHR1, 105, R, "5M", 30);
        input.fragment("reverseSameStartOtherEnd", "rg", CHR1, 100, R, "9M", 40);
        input.fragment("forwardAtReverseEnd", "rg", CHR1, 109, F, "10M", 10);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("forwardDup", "reverseDupShorter");
        assertMetrics(result.onlyMetrics(), "lib", 7, 0, 0, 0, 2, 0, 0);
    }

    @Test
    public void testFragmentsAtAPairEndAreDuplicatesWhateverTheirScore() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // quality 10 is below the threshold of 15 so the pair scores 0
        input.pair("pair", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 10);
        input.fragment("atFirstEnd", "rg", CHR1, 100, F, "10M", 40);
        input.fragment("atSecondEnd", "rg", CHR1, 291, R, "10M", 40);
        input.fragment("oppositeStrand", "rg", CHR1, 91, R, "10M", 40);
        input.fragment("oppositeStrandDup", "rg", CHR1, 91, R, "10M", 30);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("atFirstEnd", "atSecondEnd", "oppositeStrandDup");
        assertMetrics(result.onlyMetrics(), "lib", 4, 1, 0, 0, 3, 0, 0);
    }

    @Test
    public void testPairsAndFragmentsAtTheSamePosition() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.pair("pairBest", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("pairDup", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.fragment("fragment", "rg", CHR1, 100, F, "10M", 40);
        input.fragment("fragmentAtMateEnd", "rg", CHR1, 291, R, "10M", 40);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("pairDup/1", "pairDup/2", "fragment", "fragmentAtMateEnd");
        assertMetrics(result.onlyMetrics(), "lib", 2, 2, 0, 0, 2, 1, 0);
    }

    @Test
    public void testClippedBasesCountTowardTheFivePrimeEnd() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // forward fragments with an unclipped start of chr1:100
        input.fragment("forward", "rg", CHR1, 100, F, "10M", 30);
        input.fragment("softClippedStart", "rg", CHR1, 102, F, "2S8M", 20);
        input.fragment("hardClippedStart", "rg", CHR1, 103, F, "3H7M", 20);
        input.fragment("softClippedEnd", "rg", CHR1, 100, F, "8M2S", 20);
        input.fragment("alignedStartOnly", "rg", CHR1, 102, F, "8M", 40);
        // reverse fragments with an unclipped end of chr1:209
        input.fragment("reverse", "rg", CHR1, 200, R, "10M", 30);
        input.fragment("softClippedReverseEnd", "rg", CHR1, 200, R, "8M2S", 20);
        input.fragment("hardClippedReverseEnd", "rg", CHR1, 200, R, "5M5H", 20);
        input.fragment("alignedEndOnly", "rg", CHR1, 202, R, "8M2S", 40);
        // pairs with 5' ends chr1:1000 forward and chr1:1200 reverse
        input.pair("pair", "rg", CHR1, 1000, F, "10M", CHR1, 1191, R, "10M", 30);
        input.pair("clippedPair", "rg", CHR1, 1005, F, "5S5M", CHR1, 1191, R, "6M4S", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("softClippedStart", "hardClippedStart", "softClippedEnd",
                "softClippedReverseEnd", "hardClippedReverseEnd", "clippedPair/1", "clippedPair/2");
        assertMetrics(result.onlyMetrics(), "lib", 9, 2, 0, 0, 5, 1, 0);
    }

    @Test
    public void testPairWhoseFirstRecordHasTheLaterFivePrimeEnd() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // the reverse read comes first in the file but its 5' end (2009) is after its mate's (2005)
        input.pair("reverseFirst", "rg", CHR1, 2000, R, "10M", CHR1, 2005, F, "10M", 30);
        input.pair("forwardFirst", "rg", CHR1, 2005, F, "10M", CHR1, 2006, R, "4M", 20);
        // 5' ends 2005 reverse and 2009 forward
        input.pair("otherOrientation", "rg", CHR1, 1996, R, "10M", CHR1, 2009, F, "10M", 10);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("forwardFirst/1", "forwardFirst/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 3, 0, 0, 0, 1, 0);
    }

    @Test
    public void testEachPairOrientationIsItsOwnDuplicateSet() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // every pair has 5' ends at chr1:100 and chr1:300; in the second pair of each orientation the first of
        // pair is the right-hand read
        input.pair("FR1", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("FR2", "rg", CHR1, 291, R, "10M", CHR1, 100, F, "10M", 20);
        input.pair("RF1", "rg", CHR1, 91, R, "10M", CHR1, 300, F, "10M", 30);
        input.pair("RF2", "rg", CHR1, 300, F, "10M", CHR1, 91, R, "10M", 20);
        input.pair("FF1", "rg", CHR1, 100, F, "10M", CHR1, 300, F, "10M", 30);
        input.pair("FF2", "rg", CHR1, 300, F, "10M", CHR1, 100, F, "10M", 20);
        input.pair("RR1", "rg", CHR1, 91, R, "10M", CHR1, 291, R, "10M", 30);
        input.pair("RR2", "rg", CHR1, 291, R, "10M", CHR1, 91, R, "10M", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("FR2/1", "FR2/2", "RF2/1", "RF2/2", "FF2/1", "FF2/2", "RR2/1", "RR2/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 8, 0, 0, 0, 4, 0);
    }

    @Test
    public void testPairWithBothFivePrimeEndsAtOnePositionDoesNotDependOnRecordOrder() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // both ends of each pair have their 5' end at chr1:100, on opposite strands; the reverse read is first
        // in the file for one pair and second for the other
        input.pair("reverseFirst", "rg", CHR1, 100, R, "1M", CHR1, 100, F, "10M", 30);
        input.pair("forwardFirst", "rg", CHR1, 100, F, "10M", CHR1, 100, R, "1M", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("forwardFirst/1", "forwardFirst/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 2, 0, 0, 0, 1, 0);
    }

    @Test
    public void testMatesOnDifferentContigs() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.pair("best", "rg", CHR1, 100, F, "10M", CHR2, 500, R, "10M", 30);
        input.pair("dup", "rg", CHR1, 100, F, "10M", CHR2, 500, R, "10M", 20);
        input.pair("firstOfPairOnLaterContig", "rg", CHR2, 500, R, "10M", CHR1, 100, F, "10M", 20);
        input.pair("mateOnThirdContig", "rg", CHR1, 100, F, "10M", CHR3, 500, R, "10M", 20);
        input.pair("mateOneBaseLater", "rg", CHR1, 100, F, "10M", CHR2, 501, R, "10M", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("dup/1", "dup/2", "firstOfPairOnLaterContig/1", "firstOfPairOnLaterContig/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 5, 0, 0, 0, 2, 0);
    }

    @Test
    public void testUnmappedReadsAndUnmappedMatesAreNeverDuplicates() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // reads whose mate is unmapped are compared as fragments
        input.pairWithUnmappedMate("mateUnmapped1", "rg", CHR1, 100, F, "10M", 20);
        input.pairWithUnmappedMate("mateUnmapped2", "rg", CHR1, 100, F, "10M", 30);
        input.fragment("fragment", "rg", CHR1, 100, F, "10M", 10);
        input.pair("pair", "rg", CHR1, 1000, F, "10M", CHR1, 1191, R, "10M", 10);
        input.pairWithUnmappedMate("mateUnmappedAtPair", "rg", CHR1, 1000, F, "10M", 40);
        input.unplacedFragment("unplaced", "rg", 30);
        input.unplacedPair("unplacedPair", "rg", 30);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("mateUnmapped1/1", "fragment", "mateUnmappedAtPair/1");
        assertMetrics(result.onlyMetrics(), "lib", 4, 1, 0, 6, 3, 0, 0);
    }

    @Test
    public void testSecondaryAndSupplementaryRecordsAreNeverDuplicatesAndAreIgnored() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.pair("best", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        final SAMRecord[] dup = input.pair("dup", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.nonPrimary(dup[0], false, CHR2, 700, F, "10M", 40);
        input.nonPrimary(dup[1], true, CHR1, 5000, R, "10M", 40);
        input.fragment("fragment", "rg", CHR2, 2000, F, "10M", 20);
        input.fragment("fragmentDup", "rg", CHR2, 2000, F, "10M", 10);
        final SAMRecord other = input.fragment("other", "rg", CHR3, 100, F, "10M", 20);
        // high-scoring non-primary records at the positions of lower-scoring fragments
        input.nonPrimary(other, true, CHR2, 2000, F, "10M", 40);
        input.nonPrimary(other, false, CHR3, 3000, F, "10M", 40);
        input.fragment("aloneWithSecondary", "rg", CHR3, 3000, F, "10M", 10);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("dup/1", "dup/2", "fragmentDup");
        assertMetrics(result.onlyMetrics(), "lib", 4, 2, 4, 0, 1, 1, 0);
    }

    @Test
    public void testDuplicateFlagsInTheInputAreRecomputed() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        final List<SAMRecord> all = new ArrayList<>();
        all.addAll(Arrays.asList(input.pair("best", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30)));
        all.addAll(Arrays.asList(input.pair("dup", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20)));
        all.add(input.nonPrimary(all.get(0), false, CHR2, 100, F, "10M", 30));
        all.add(input.fragment("unique", "rg", CHR2, 500, F, "10M", 30));
        all.addAll(Arrays.asList(input.pairWithUnmappedMate("mateUnmapped", "rg", CHR3, 100, F, "10M", 30)));
        all.add(input.unplacedFragment("unplaced", "rg", 30));
        all.forEach(r -> r.setDuplicateReadFlag(true));

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("dup/1", "dup/2");
        assertMetrics(result.onlyMetrics(), "lib", 2, 2, 1, 2, 0, 1, 0);
    }

    /**
     * A read whose mate is flagged as mapped but is missing from the file counts as a pair end: a fragment at its 5'
     * end is a duplicate, and two such reads make one pair in READ_PAIRS_EXAMINED. The documentation does not specify
     * this; the expected value follows the current implementation.
     */
    @Test
    public void testReadsWhoseMateIsMissingCountAsPairEnds() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // each orphan says its mate is mapped, but the mate is not in the file
        final SAMRecord[] firstOrphan = input.pair("firstOrphan", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 10);
        input.remove(firstOrphan[1]);
        final SAMRecord[] secondOrphan = input.pair("secondOrphan", "rg", CHR2, 100, F, "10M", CHR2, 291, R, "10M", 10);
        input.remove(secondOrphan[0]);
        input.fragment("atFirstOrphan", "rg", CHR1, 100, F, "10M", 40);
        input.fragment("atSecondOrphan", "rg", CHR2, 291, R, "10M", 40);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("atFirstOrphan", "atSecondOrphan");
        assertMetrics(result.onlyMetrics(), "lib", 2, 1, 0, 0, 2, 0, 0);
    }

    @Test
    public void testLibrariesAndReadGroups() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner()
                .readGroup("rgA", "libX")
                .readGroup("rgB", "libX")
                .readGroup("rgC", "libY")
                .readGroup("rgD", null)
                .readGroup("rgE", "libZ")
                .readGroup("rgF", "libW");
        input.pair("a", "rgA", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("b", "rgB", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.pair("c", "rgC", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 10);
        input.pair("d", "rgD", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 10);
        input.fragment("fa", "rgA", CHR2, 100, F, "10M", 30);
        input.fragment("fb", "rgB", CHR2, 100, F, "10M", 20);
        input.fragment("fc", "rgC", CHR2, 100, F, "10M", 40);
        // at the first end of the pairs, but the only read of its library there
        input.fragment("fe", "rgE", CHR1, 100, F, "10M", 10);
        input.unplacedFragment("unplaced", "rgD", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("b/1", "b/2", "fb");
        Assert.assertEquals(result.metricsFile.getMetrics().stream().map(m -> m.LIBRARY).collect(Collectors.toList()),
                Arrays.asList("Unknown Library", "libW", "libX", "libY", "libZ"));
        assertMetrics(result.metrics("Unknown Library"), "Unknown Library", 0, 1, 0, 1, 0, 0, 0);
        assertMetrics(result.metrics("libX"), "libX", 2, 2, 0, 0, 1, 1, 0);
        assertMetrics(result.metrics("libY"), "libY", 1, 1, 0, 0, 0, 0, 0);
        assertMetrics(result.metrics("libZ"), "libZ", 1, 0, 0, 0, 0, 0, 0);
        assertMetrics(result.metrics("libW"), "libW", 0, 0, 0, 0, 0, 0, 0);
        Assert.assertNull(result.histogram("CoverageMult"));
    }

    @Test
    public void testPairsWithTheSameNameInDifferentReadGroups() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner()
                .readGroup("rgA", "lib")
                .readGroup("rgB", "lib");
        // the two pairs named "same" interleave; each must be matched to its mate in its own read group
        input.pair("same", "rgA", CHR1, 1000, F, "10M", CHR1, 1191, R, "10M", 30);
        input.pair("same", "rgB", CHR1, 1100, F, "10M", CHR1, 1291, R, "10M", 30);
        input.pair("dupOfB", "rgB", CHR1, 1100, F, "10M", CHR1, 1291, R, "10M", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("dupOfB/1", "dupOfB/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 3, 0, 0, 0, 1, 0);
    }

    @Test
    public void testSeveralInputsAreMergedAndMarkedTogether() {
        final MarkDuplicatesEndToEndRunner first = new MarkDuplicatesEndToEndRunner().readGroup("rg1", "lib");
        first.pair("p1", "rg1", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        first.fragment("f1", "rg1", CHR2, 100, F, "10M", 30);
        final MarkDuplicatesEndToEndRunner second = new MarkDuplicatesEndToEndRunner().readGroup("rg2", "lib");
        second.pair("p2", "rg2", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        second.fragment("f2", "rg2", CHR2, 100, F, "10M", 20);
        second.pair("q2", "rg2", CHR1, 50, F, "10M", CHR3, 100, R, "10M", 20);

        final File dir = IOUtil.createTempDir("MarkDuplicatesSeveralInputs").toFile();
        final Result result;
        try {
            final File firstFile = first.write(new File(dir, "first.sam"));
            final File secondFile = second.write(new File(dir, "second.sam"));
            result = MarkDuplicatesEndToEndRunner.runOn(dir, Arrays.asList(firstFile, secondFile));
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
        Assert.assertEquals(result.exitCode, 0);
        result.assertDuplicates("p1/1", "p1/2", "f2");
        assertMetrics(result.onlyMetrics(), "lib", 2, 3, 0, 0, 1, 1, 0);
        Assert.assertNotNull(result.header.getReadGroup("rg1"));
        Assert.assertNotNull(result.header.getReadGroup("rg2"));
        assertMergedRecords(result, first, second);
    }

    /** Asserts that every input record is written once, unchanged apart from the flag and PG tag, in coordinate order. */
    private static void assertMergedRecords(final Result result, final MarkDuplicatesEndToEndRunner... inputs) {
        final Map<String, SAMRecord> inputByKey = new HashMap<>();
        for (final MarkDuplicatesEndToEndRunner input : inputs) {
            for (final SAMRecord rec : input.inputRecords()) {
                inputByKey.put(MarkDuplicatesEndToEndRunner.key(rec), rec);
            }
        }
        Assert.assertEquals(result.records.size(), inputByKey.size());
        int previousContig = -1, previousStart = -1;
        for (final SAMRecord rec : result.records) {
            final SAMRecord in = inputByKey.get(MarkDuplicatesEndToEndRunner.key(rec));
            Assert.assertNotNull(in, MarkDuplicatesEndToEndRunner.key(rec));
            Assert.assertEquals(MarkDuplicatesEndToEndRunner.withoutDuplicateFlagAndTags(rec, SAMTag.PG.name()),
                    MarkDuplicatesEndToEndRunner.withoutDuplicateFlagAndTags(in, SAMTag.PG.name()));
            Assert.assertTrue(rec.getReferenceIndex() > previousContig ||
                    rec.getReferenceIndex() == previousContig && rec.getAlignmentStart() >= previousStart);
            previousContig = rec.getReferenceIndex();
            previousStart = rec.getAlignmentStart();
        }
    }

    @Test
    public void testSeveralInputsAssumedCoordinateSorted() {
        final MarkDuplicatesEndToEndRunner first = new MarkDuplicatesEndToEndRunner(SAMFileHeader.SortOrder.unsorted).readGroup("rg1", "lib");
        first.pair("p1", "rg1", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        first.fragment("f1", "rg1", CHR2, 100, F, "10M", 30);
        final MarkDuplicatesEndToEndRunner second = new MarkDuplicatesEndToEndRunner(SAMFileHeader.SortOrder.unsorted).readGroup("rg2", "lib");
        second.pair("p2", "rg2", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        second.fragment("f2", "rg2", CHR2, 100, F, "10M", 20);

        final File dir = IOUtil.createTempDir("MarkDuplicatesSeveralInputs").toFile();
        final Result result;
        try {
            final File firstFile = first.write(new File(dir, "first.sam"));
            final File secondFile = second.write(new File(dir, "second.sam"));
            result = MarkDuplicatesEndToEndRunner.runOn(dir, Arrays.asList(firstFile, secondFile), "ASSUME_SORTED=true");
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
        Assert.assertEquals(result.exitCode, 0);
        Assert.assertEquals(result.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        result.assertDuplicates("p1/1", "p1/2", "f2");
        // merged in coordinate order
        Assert.assertEquals(result.records.stream().map(SAMRecord::getAlignmentStart).collect(Collectors.toList()),
                Arrays.asList(100, 100, 291, 291, 100, 100));
        Assert.assertEquals(result.records.stream().map(SAMRecord::getReferenceIndex).collect(Collectors.toList()),
                Arrays.asList(CHR1, CHR1, CHR1, CHR1, CHR2, CHR2));
        assertMergedRecords(result, first, second);
        assertMetrics(result.onlyMetrics(), "lib", 2, 2, 0, 0, 1, 1, 0);
    }

    /**
     * With one library the metrics file holds the CoverageMult histogram: the unique pairs expected from x times
     * the sequencing, X * (1 - exp(-x * pairs / X)) for library size X, relative to the unique pairs observed.
     */
    @Test
    public void testCoverageMultHistogram() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // sets of three, two and one pairs; fragments do not enter the histograms
        for (int i = 0; i < 3; i++) {
            input.pair("a" + i, "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30 - i);
        }
        for (int i = 0; i < 2; i++) {
            input.pair("b" + i, "rg", CHR1, 1000, F, "10M", CHR1, 1191, R, "10M", 30 - i);
        }
        input.pair("c", "rg", CHR2, 100, F, "10M", CHR2, 291, R, "10M", 30);
        input.fragment("f", "rg", CHR3, 100, F, "10M", 30);
        input.fragment("g", "rg", CHR3, 100, F, "10M", 20);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("a1/1", "a1/2", "a2/1", "a2/2", "b1/1", "b1/2", "g");
        assertMetrics(result.onlyMetrics(), "lib", 2, 6, 0, 0, 1, 3, 0);
        assertCoverageMultHistogram(result, 6, 3, 0);
        result.assertSetSizeHistograms(new double[]{1, 1, 2, 1, 3, 1}, new double[0], new double[]{1, 1, 2, 1, 3, 1});
    }

    static void assertCoverageMultHistogram(final Result result, final long pairs, final long duplicatePairs, final long opticalPairs) {
        final long librarySize = MarkDuplicatesEndToEndRunner.expectedLibrarySize(pairs - opticalPairs, pairs - duplicatePairs);
        final Histogram<Double> histogram = result.histogram("CoverageMult");
        Assert.assertNotNull(histogram);
        Assert.assertEquals(histogram.values().size(), 100);
        for (int x = 1; x <= 100; x++) {
            final double expected = librarySize * (1 - Math.exp(-x * (double) pairs / librarySize)) / (pairs - duplicatePairs);
            Assert.assertEquals(histogram.get((double) x).getValue(), expected, 1e-5, "x=" + x);
        }
    }

    @Test
    public void testEmptyInput() {
        final Result result = runAndCheckRecords(oneLibrary());
        Assert.assertTrue(result.records.isEmpty());
        assertMetrics(result.onlyMetrics(), "lib", 0, 0, 0, 0, 0, 0, 0);
        Assert.assertNull(result.histogram("CoverageMult"));
        result.assertSetSizeHistograms(new double[0], new double[0], new double[0]);
    }

    @Test
    public void testEmptyInputWithoutReadGroupsHasNoMetricsRows() {
        final Result result = runAndCheckRecords(new MarkDuplicatesEndToEndRunner());
        Assert.assertTrue(result.records.isEmpty());
        Assert.assertTrue(result.metricsFile.getMetrics().isEmpty());
    }

    @Test
    public void testSingleFragment() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.fragment("only", "rg", CHR1, 100, F, "10M", 30);
        final Result result = runAndCheckRecords(input);
        result.assertDuplicates();
        assertMetrics(result.onlyMetrics(), "lib", 1, 0, 0, 0, 0, 0, 0);
    }

    @Test
    public void testSinglePair() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.pair("only", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        final Result result = runAndCheckRecords(input);
        result.assertDuplicates();
        assertMetrics(result.onlyMetrics(), "lib", 0, 1, 0, 0, 0, 0, 0);
        Assert.assertNull(result.histogram("CoverageMult"));
        result.assertSetSizeHistograms(new double[]{1, 1}, new double[0], new double[]{1, 1});
    }

    @Test
    public void testOnlyUnmappedReads() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        input.unplacedFragment("a", "rg", 30);
        input.unplacedPair("b", "rg", 30);
        final Result result = runAndCheckRecords(input);
        result.assertDuplicates();
        assertMetrics(result.onlyMetrics(), "lib", 0, 0, 0, 3, 0, 0, 0);
    }

    @DataProvider(name = "scoringStrategies")
    public Object[][] scoringStrategies() {
        return new Object[][]{
                // no argument: the MarkDuplicates default is SUM_OF_BASE_QUALITIES
                {null, "longReference"},
                {ScoringStrategy.SUM_OF_BASE_QUALITIES, "longReference"},
                {ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, "highQuality"},
        };
    }

    @Test(dataProvider = "scoringStrategies")
    public void testScoringStrategyChoosesTheKeeper(final ScoringStrategy strategy, final String expectedDuplicate) {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // base quality sums 600 and 400; reference lengths 20 and 40; same 5' ends
        input.pair("highQuality", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("longReference", "rg", CHR1, 100, F, "5M10D5M", CHR1, 281, R, "5M10D5M", 20);

        final Result result = strategy == null ? runAndCheckRecords(input) :
                runAndCheckRecords(input, "DUPLICATE_SCORING_STRATEGY=" + strategy.name());
        result.assertDuplicates(expectedDuplicate + "/1", expectedDuplicate + "/2");
        assertMetrics(result.onlyMetrics(), "lib", 0, 2, 0, 0, 0, 1, 0);
    }

    @Test
    public void testBaseQualitiesBelowFifteenDoNotCount() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // 20 bases of quality 14 score 0; 4 bases of quality 15 score 60
        input.pair("manyBasesQuality14", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 14);
        input.pair("fewBasesQuality15", "rg", CHR1, 100, F, "2M", CHR1, 299, R, "2M", 15);
        input.fragment("fragmentQuality14", "rg", CHR2, 100, F, "10M", 14);
        input.fragment("fragmentQuality15", "rg", CHR2, 100, F, "1M", 15);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("manyBasesQuality14/1", "manyBasesQuality14/2", "fragmentQuality14");
    }

    @DataProvider(name = "allScoringStrategies")
    public Object[][] allScoringStrategies() {
        return Arrays.stream(ScoringStrategy.values()).map(s -> new Object[]{s}).toArray(Object[][]::new);
    }

    @Test(dataProvider = "allScoringStrategies")
    public void testReadsFailingVendorQualityLose(final ScoringStrategy strategy) {
        // names whose RANDOM scores are at most 16383, so that the two ends of a pair sum without overflow
        final List<String> names = namesWithRandomScoreAtMost(16383, 4);
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        final SAMRecord[] failing = input.pair(names.get(0), "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 40);
        failing[1].setReadFailsVendorQualityCheckFlag(true);
        input.pair(names.get(1), "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.fragment(names.get(2), "rg", CHR2, 100, F, "10M", 40).setReadFailsVendorQualityCheckFlag(true);
        input.fragment(names.get(3), "rg", CHR2, 100, F, "10M", 20);

        final Result result = runAndCheckRecords(input, "DUPLICATE_SCORING_STRATEGY=" + strategy.name());
        result.assertDuplicates(names.get(0) + "/1", names.get(0) + "/2", names.get(2));
    }

    private static List<String> namesWithRandomScoreAtMost(final int maxScore, final int count) {
        final List<String> names = new ArrayList<>();
        for (int i = 0; names.size() < count; i++) {
            final String name = "read" + i;
            if (randomScore(name) <= maxScore) {
                names.add(name);
            }
        }
        return names;
    }

    private static short randomScore(final String name) {
        final SAMRecord rec = new SAMRecord(null);
        rec.setReadName(name);
        return DuplicateScoringStrategy.computeDuplicateScore(rec, ScoringStrategy.RANDOM);
    }

    @Test
    public void testRandomScoringKeepsTheHigherScoringPairAndFragment() {
        final List<String> names = namesWithRandomScoreAtMost(16383, 4);
        Assert.assertNotEquals(randomScore(names.get(0)), randomScore(names.get(1)));
        Assert.assertNotEquals(randomScore(names.get(2)), randomScore(names.get(3)));
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        for (int i = 0; i < 2; i++) {
            input.pair(names.get(i), "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
            input.fragment(names.get(i + 2), "rg", CHR2, 100, F, "10M", 30);
        }
        // the score of a pair is the sum of the scores of its ends, which share a name
        final String pairDup = randomScore(names.get(0)) < randomScore(names.get(1)) ? names.get(0) : names.get(1);
        final String fragmentDup = randomScore(names.get(2)) < randomScore(names.get(3)) ? names.get(2) : names.get(3);

        final Result result = runAndCheckRecords(input, "DUPLICATE_SCORING_STRATEGY=RANDOM");
        result.assertDuplicates(pairDup + "/1", pairDup + "/2", fragmentDup);
    }

    /**
     * Under every scoring strategy, of two pairs or two fragments with equal scores and no location in their names,
     * the one first in the file is kept. The documentation does not specify this; the expected value follows the
     * current implementation.
     */
    @Test(dataProvider = "allScoringStrategies")
    public void testTiedScoresKeepTheRecordFirstInTheFile(final ScoringStrategy strategy) {
        // RANDOM scores depend only on the name, so the tied records share a name in different read groups
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner()
                .readGroup("rgA", "lib")
                .readGroup("rgB", "lib");
        input.pair("tied", "rgB", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("tied", "rgA", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.fragment("tiedFragment", "rgB", CHR2, 100, F, "10M", 30);
        input.fragment("tiedFragment", "rgA", CHR2, 100, F, "10M", 30);

        final Result result = runAndCheckRecords(input, "DUPLICATE_SCORING_STRATEGY=" + strategy.name());
        final List<String> duplicates = result.records.stream().filter(SAMRecord::getDuplicateReadFlag)
                .map(r -> MarkDuplicatesEndToEndRunner.key(r) + "@" + r.getStringAttribute(SAMTag.RG.name()))
                .collect(Collectors.toList());
        Assert.assertEquals(duplicates, Arrays.asList("tied/1@rgA", "tied/2@rgA", "tiedFragment@rgA"));
    }

    /**
     * Of pairs or fragments with equal scores, the one with the lowest tile, then x, then y is kept, wherever it is in
     * the file. The documentation does not specify this; the expected value follows the current implementation.
     */
    @Test
    public void testTiedScoresKeepTheLowestTileThenXThenY() {
        final MarkDuplicatesEndToEndRunner input = oneLibrary();
        // read names give tile, x and y; the records later in the file have the lower values
        input.pair("RUN:1:2:5:5", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("RUN:1:1:9000:9000", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.fragment("RUN:1:3:700:900", "rg", CHR2, 100, F, "10M", 30);
        input.fragment("RUN:1:3:500:900", "rg", CHR2, 100, F, "10M", 30);
        input.fragment("RUN:1:4:50:60", "rg", CHR3, 100, F, "10M", 30);
        input.fragment("RUN:1:4:50:40", "rg", CHR3, 100, F, "10M", 30);

        final Result result = runAndCheckRecords(input);
        result.assertDuplicates("RUN:1:2:5:5/1", "RUN:1:2:5:5/2", "RUN:1:3:700:900", "RUN:1:4:50:60");
        // fragments are never optical duplicates and the pairs are on different tiles
        assertMetrics(result.onlyMetrics(), "lib", 4, 2, 0, 0, 2, 1, 0);
    }
}
