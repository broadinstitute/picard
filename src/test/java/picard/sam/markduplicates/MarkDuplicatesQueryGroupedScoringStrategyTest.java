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
import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.List;

import static picard.sam.markduplicates.MarkDuplicatesQueryNameTestSupport.*;

/**
 * Tests DUPLICATE_SCORING_STRATEGY TOTAL_MAPPED_REFERENCE_LENGTH and RANDOM in {@link MarkDuplicates} on query-grouped
 * input (SO:unsorted, GO:query, ASSUME_SORT_ORDER=queryname). A pair scores the sum of the scores of its two primary
 * reads and a fragment the score of its read; the highest score is kept. Each input is built so that the default
 * SUM_OF_BASE_QUALITIES would keep another pair or fragment, and so that the pair kept is not the first in the file.
 *
 * Pairs have read one forward with its 5' end at chr1:100 and read two reverse with its 5' end at chr1:309; fragments
 * are forward with their 5' end at chr1:2000. Every read has ten bases.
 */
public class MarkDuplicatesQueryGroupedScoringStrategyTest {
    private static final int QC_FAIL = 0x200;
    private static final String[] QUERY_GROUPED = {"READ_NAME_REGEX=null", "ASSUME_SORT_ORDER=queryname"};

    /** Adds read one at chr1:pos1 forward and read two at chr1:pos2 reverse, each with the other's CIGAR in MC. */
    private static void pair(final SamText input, final String name, final String cigar1, final int pos2, final String cigar2,
                             final char qual, final int read2ExtraFlags) {
        input.record(name, PAIRED | PROPER | FIRST | MATE_REVERSE, "chr1", 100, cigar1, "chr1", pos2, qual,
                "MC:Z:" + cigar2, "RG:Z:rgA");
        input.record(name, PAIRED | PROPER | SECOND | REVERSE | read2ExtraFlags, "chr1", pos2, cigar2, "chr1", 100, qual,
                "MC:Z:" + cigar1, "RG:Z:rgA");
    }

    private static void fragment(final SamText input, final String name, final String cigar, final char qual) {
        input.record(name, 0, "chr1", 2000, cigar, "*", 0, qual, "RG:Z:rgA");
    }

    @DataProvider(name = "referenceLengthStrategies")
    public Object[][] referenceLengthStrategies() {
        return new Object[][]{
                {"TOTAL_MAPPED_REFERENCE_LENGTH", "DDD..DDD."},
                {"SUM_OF_BASE_QUALITIES", "...DDDD.D"},
        };
    }

    /**
     * Mapped reference lengths of read one + read two: q1 10 + 10, q2 15 + 10, q3 10 + 14, so TOTAL_MAPPED_REFERENCE_LENGTH
     * keeps q2 although q1 has the highest base qualities ('I', against '5' for q2 and '?' for q3). q1's supplementary
     * alignment (110 bases of reference) does not count towards its score, and is marked with q1. Fragment r2 (16
     * bases of reference, quality '5') is kept over r1 (10 bases, quality 'I').
     */
    @Test(dataProvider = "referenceLengthStrategies")
    public void testTotalMappedReferenceLength(final String strategy, final String expectedFlags) throws IOException {
        final SamText input = new SamText("SO:unsorted\tGO:query").readGroup("rgA", "libA");
        pair(input, "q1", "10M", 300, "10M", 'I', 0);
        input.record("q1", PAIRED | FIRST | SUPPLEMENTARY | MATE_REVERSE, "chr2", 5000, "5M100D5M", "chr1", 300, 'I',
                "MC:Z:10M", "RG:Z:rgA");
        pair(input, "q2", "5M5D5M", 300, "10M", '5', 0);
        pair(input, "q3", "10M", 296, "3M4D7M", '?', 0);
        fragment(input, "r1", "10M", 'I');
        fragment(input, "r2", "4M6D6M", '5');

        final Result result = run(input.write(), args(QUERY_GROUPED, "DUPLICATE_SCORING_STRATEGY=" + strategy));

        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        assertDuplicateFlags(result, expectedFlags);
        assertMetrics(result.metricsFor("libA"), 2, 3, 1, 0, 1, 2, 0);
    }

    /** The RANDOM score htsjdk gives one read, which depends only on its name and its vendor quality flag. */
    private static int randomScore(final String name, final boolean failsVendorQuality) {
        final SAMRecord rec = new SAMRecord(null);
        rec.setReadName(name);
        rec.setReadFailsVendorQualityCheckFlag(failsVendorQuality);
        return DuplicateScoringStrategy.computeDuplicateScore(rec, DuplicateScoringStrategy.ScoringStrategy.RANDOM);
    }

    @DataProvider(name = "randomInputs")
    public Object[][] randomInputs() {
        return new Object[][]{{false, "q3"}, {true, "r5"}};
    }

    /**
     * Names are chosen so that no read scores above 16383 and the two reads of a pair cannot overflow a short when
     * summed. Of pairs r2 ('I', 5M5D5M), q3 ('5') and r5 ('?') the one with the highest sum of its reads' RANDOM scores
     * is kept, which is q3, or r5 when read two of q3 fails vendor quality checks; of fragments pA ('I') and pC ('5'),
     * pC.
     */
    @Test(dataProvider = "randomInputs")
    public void testRandom(final boolean q3FailsVendorQuality, final String expectedKeeper) throws IOException {
        final List<String> pairs = List.of("r2", "q3", "r5");
        final List<String> fragments = List.of("pA", "pC");
        String keeper = null;
        int best = Integer.MIN_VALUE;
        for (final String name : pairs) {
            Assert.assertTrue(randomScore(name, false) <= 16383, name);
            final int score = randomScore(name, false) + randomScore(name, q3FailsVendorQuality && name.equals("q3"));
            if (score > best) {
                best = score;
                keeper = name;
            }
        }
        Assert.assertEquals(keeper, expectedKeeper, "RANDOM scores of the read names");
        Assert.assertTrue(randomScore("pC", false) > randomScore("pA", false), "RANDOM scores of the fragment names");

        final SamText input = new SamText("SO:unsorted\tGO:query").readGroup("rgA", "libA");
        pair(input, "r2", "5M5D5M", 300, "10M", 'I', 0);
        pair(input, "q3", "10M", 300, "10M", '5', q3FailsVendorQuality ? QC_FAIL : 0);
        pair(input, "r5", "10M", 300, "10M", '?', 0);
        fragment(input, fragments.get(0), "10M", 'I');
        fragment(input, fragments.get(1), "10M", '5');

        final Result result = run(input.write(), args(QUERY_GROUPED, "DUPLICATE_SCORING_STRATEGY=RANDOM"));

        final StringBuilder expectedFlags = new StringBuilder();
        for (final String name : pairs) {
            expectedFlags.append(name.equals(expectedKeeper) ? ".." : "DD");
        }
        expectedFlags.append("D.");
        assertOnlyDuplicateFlagAndTagsChanged(input.recordLines(), result.recordLines);
        assertDuplicateFlags(result, expectedFlags.toString());
        assertMetrics(result.metricsFor("libA"), 2, 3, 0, 0, 1, 2, 0);
    }
}
