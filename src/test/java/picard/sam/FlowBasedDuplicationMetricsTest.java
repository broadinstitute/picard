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

package picard.sam;

import htsjdk.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.MarkDuplicatesForFlowHelper;

import java.util.Collections;

import static picard.sam.DuplicationMetricsBehaviorTest.MATE_UNMAPPED;
import static picard.sam.DuplicationMetricsBehaviorTest.PAIRED;
import static picard.sam.DuplicationMetricsBehaviorTest.PROPER;
import static picard.sam.DuplicationMetricsBehaviorTest.FIRST;
import static picard.sam.DuplicationMetricsBehaviorTest.REVERSE;
import static picard.sam.DuplicationMetricsBehaviorTest.UNMAPPED;
import static picard.sam.DuplicationMetricsBehaviorTest.DUPLICATE;

/**
 * Tests {@link FlowBasedDuplicationMetrics}: which single-end reads count as having a known fragment length (clipped
 * at the read's 3' end, or adapter clipped according to the tm tag), the counts kept for duplicates with and without a
 * known length, the two derived duplication rates including empty denominators, and merging.
 */
public class FlowBasedDuplicationMetricsTest {

    private static SAMRecord record(final int flags, final String cigar, final String clippingTag) {
        final SAMRecord rec = DuplicationMetricsBehaviorTest.record(flags, cigar);
        if (clippingTag != null) {
            rec.setAttribute(MarkDuplicatesForFlowHelper.CLIPPING_TAG_NAME, clippingTag);
        }
        return rec;
    }

    @DataProvider(name = "knownFragmentData")
    public Object[][] knownFragmentData() {
        return new Object[][]{
                // flags, cigar, tm tag, counted as a single-end read with known fragment length
                {0, "10M", null, false},
                {0, "8M2S", null, true},
                {0, "2S8M", null, false},
                {0, "2S6M2S", null, true},
                {0, "8M2H", null, true},
                {0, "2H8M", null, false},
                {REVERSE, "2S8M", null, true},
                {REVERSE, "8M2S", null, false},
                {REVERSE, "2H8M", null, true},
                {REVERSE, "10M", null, false},
                {0, "10M", "A", true},
                {REVERSE, "10M", "A", true},
                {0, "10M", "QA", true},
                {0, "10M", "Q", false},
                {0, "10M", "Z", false},
                {0, "10M", "QZ", false},
                {0, "10M", "", false},
                {DUPLICATE, "8M2S", null, true},
                {PAIRED | MATE_UNMAPPED | FIRST, "8M2S", null, false},
                {PAIRED | MATE_UNMAPPED | FIRST, "10M", "A", false},
                {PAIRED | PROPER | FIRST, "8M2S", "A", false},
                {UNMAPPED, "8M2S", "A", false},
        };
    }

    @Test(dataProvider = "knownFragmentData")
    public void testAddReadToLibraryMetrics(final int flags, final String cigar, final String clippingTag, final boolean knownFragment) {
        final FlowBasedDuplicationMetrics metrics = new FlowBasedDuplicationMetrics();
        final SAMRecord rec = record(flags, cigar, clippingTag);

        metrics.addReadToLibraryMetrics(rec);
        Assert.assertEquals(metrics.UNPAIRED_WITH_TLEN, knownFragment ? 1 : 0);
        Assert.assertEquals(metrics.UNPAIRED_DUPS_WITH_TLEN, 0);
        Assert.assertEquals(metrics.UNPAIRED_DUPS_WITHOUT_TLEN, 0);

        final DuplicationMetrics plain = new DuplicationMetrics();
        plain.addReadToLibraryMetrics(rec);
        Assert.assertEquals(DuplicationMetricsBehaviorTest.counts(metrics), DuplicationMetricsBehaviorTest.counts(plain));

        metrics.addReadToLibraryMetrics(rec);
        Assert.assertEquals(metrics.UNPAIRED_WITH_TLEN, knownFragment ? 2 : 0);
        Assert.assertNull(metrics.UNPAIRED_DUP_RATE_WITH_TLEN);
        Assert.assertNull(metrics.UNPAIRED_DUP_RATE_WITHOUT_TLEN);
    }

    @Test(dataProvider = "knownFragmentData")
    public void testAddDuplicateReadToMetrics(final int flags, final String cigar, final String clippingTag, final boolean knownFragment) {
        final FlowBasedDuplicationMetrics metrics = new FlowBasedDuplicationMetrics();
        final SAMRecord rec = record(flags, cigar, clippingTag);

        metrics.addDuplicateReadToMetrics(rec);
        final boolean unpairedDuplicate = (flags & UNMAPPED) == 0 && ((flags & PAIRED) == 0 || (flags & MATE_UNMAPPED) != 0);
        Assert.assertEquals(metrics.UNPAIRED_DUPS_WITH_TLEN, unpairedDuplicate && knownFragment ? 1 : 0);
        Assert.assertEquals(metrics.UNPAIRED_DUPS_WITHOUT_TLEN, unpairedDuplicate && !knownFragment ? 1 : 0);
        Assert.assertEquals(metrics.UNPAIRED_WITH_TLEN, 0);

        final DuplicationMetrics plain = new DuplicationMetrics();
        plain.addDuplicateReadToMetrics(rec);
        Assert.assertEquals(DuplicationMetricsBehaviorTest.counts(metrics), DuplicationMetricsBehaviorTest.counts(plain));
        Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, unpairedDuplicate ? 1 : 0);
    }

    @Test
    public void testDuplicatesOfSecondaryAndSupplementaryReadsAreNotCounted() {
        for (final int flags : new int[]{0x100, 0x800, 0x100 | 0x800}) {
            final FlowBasedDuplicationMetrics metrics = new FlowBasedDuplicationMetrics();
            metrics.addDuplicateReadToMetrics(record(flags, "8M2S", "A"));
            metrics.addDuplicateReadToMetrics(record(flags, "10M", null));
            Assert.assertEquals(metrics.UNPAIRED_DUPS_WITH_TLEN, 0);
            Assert.assertEquals(metrics.UNPAIRED_DUPS_WITHOUT_TLEN, 0);
            Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, 0);
        }
    }

    private static FlowBasedDuplicationMetrics flowMetrics(final long unpairedExamined, final long withLength,
                                                           final long duplicatesWithLength, final long duplicatesWithoutLength) {
        final FlowBasedDuplicationMetrics metrics = new FlowBasedDuplicationMetrics();
        metrics.LIBRARY = "lib";
        metrics.UNPAIRED_READS_EXAMINED = unpairedExamined;
        metrics.UNPAIRED_WITH_TLEN = withLength;
        metrics.UNPAIRED_DUPS_WITH_TLEN = duplicatesWithLength;
        metrics.UNPAIRED_DUPS_WITHOUT_TLEN = duplicatesWithoutLength;
        metrics.UNPAIRED_READ_DUPLICATES = duplicatesWithLength + duplicatesWithoutLength;
        return metrics;
    }

    @DataProvider(name = "derivedRatesData")
    public Object[][] derivedRatesData() {
        return new Object[][]{
                // unpaired examined, with known length, duplicates with, duplicates without, rate without, rate with
                {10L, 4L, 1L, 3L, 0.5, 0.25},
                {0L, 0L, 0L, 0L, 0.0, 0.0},
                {4L, 4L, 2L, 0L, 0.0, 0.5},
                {5L, 0L, 0L, 5L, 1.0, 0.0},
                {7L, 4L, 3L, 1L, 1.0 / 3.0, 0.75},
                {1000000L, 1L, 1L, 999L, 999.0 / 999999.0, 1.0},
        };
    }

    @Test(dataProvider = "derivedRatesData")
    public void testCalculateDerivedFields(final long unpairedExamined, final long withLength, final long duplicatesWithLength,
                                           final long duplicatesWithoutLength, final double rateWithout, final double rateWith) {
        final FlowBasedDuplicationMetrics metrics = flowMetrics(unpairedExamined, withLength, duplicatesWithLength, duplicatesWithoutLength);
        metrics.calculateDerivedFields();
        Assert.assertEquals(metrics.UNPAIRED_DUP_RATE_WITHOUT_TLEN.doubleValue(), rateWithout, 0.0);
        Assert.assertEquals(metrics.UNPAIRED_DUP_RATE_WITH_TLEN.doubleValue(), rateWith, 0.0);

        // the fields of the parent class are derived as before
        final double percent = unpairedExamined == 0 ? 0.0 : (duplicatesWithLength + duplicatesWithoutLength) / (double) unpairedExamined;
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), percent, 0.0);
        Assert.assertNull(metrics.ESTIMATED_LIBRARY_SIZE);
    }

    @Test
    public void testDerivedFieldsIncludeLibrarySize() {
        final FlowBasedDuplicationMetrics metrics = flowMetrics(10, 4, 1, 3);
        metrics.READ_PAIRS_EXAMINED = 100;
        metrics.READ_PAIR_DUPLICATES = 10;
        metrics.calculateDerivedFields();
        Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE.longValue(), DuplicationMetricsBehaviorTest.expectedLibrarySize(100, 90));
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), 24.0 / 210.0, 0.0);
        Assert.assertEquals(metrics.UNPAIRED_DUP_RATE_WITHOUT_TLEN.doubleValue(), 0.5, 0.0);
        Assert.assertEquals(metrics.UNPAIRED_DUP_RATE_WITH_TLEN.doubleValue(), 0.25, 0.0);
    }

    /** Counts gathered read by read lead to the rates computed by hand. */
    @Test
    public void testAddReadsThenDerive() {
        final FlowBasedDuplicationMetrics metrics = new FlowBasedDuplicationMetrics();
        final SAMRecord[] withLength = {record(0, "8M2S", null), record(REVERSE, "2S8M", null), record(0, "10M", "A")};
        final SAMRecord[] withoutLength = {record(0, "10M", null), record(REVERSE, "10M", null),
                record(0, "2S8M", null), record(PAIRED | MATE_UNMAPPED | FIRST, "8M2S", null)};
        for (final SAMRecord rec : withLength) metrics.addReadToLibraryMetrics(rec);
        for (final SAMRecord rec : withoutLength) metrics.addReadToLibraryMetrics(rec);
        metrics.addDuplicateReadToMetrics(withLength[0]);
        metrics.addDuplicateReadToMetrics(withoutLength[0]);
        metrics.addDuplicateReadToMetrics(withoutLength[3]);

        Assert.assertEquals(metrics.UNPAIRED_READS_EXAMINED, 7);
        Assert.assertEquals(metrics.UNPAIRED_WITH_TLEN, 3);
        Assert.assertEquals(metrics.UNPAIRED_DUPS_WITH_TLEN, 1);
        Assert.assertEquals(metrics.UNPAIRED_DUPS_WITHOUT_TLEN, 2);
        Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, 3);

        metrics.calculateDerivedFields();
        Assert.assertEquals(metrics.UNPAIRED_DUP_RATE_WITH_TLEN.doubleValue(), 1.0 / 3.0, 0.0);
        Assert.assertEquals(metrics.UNPAIRED_DUP_RATE_WITHOUT_TLEN.doubleValue(), 0.5, 0.0);
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), 3.0 / 7.0, 0.0);
    }

    @Test
    public void testMerge() {
        final FlowBasedDuplicationMetrics left = flowMetrics(10, 4, 1, 3);
        final FlowBasedDuplicationMetrics right = flowMetrics(6, 2, 2, 0);
        right.READ_PAIRS_EXAMINED = 5;
        right.UNMAPPED_READS = 9;

        left.merge(Collections.singletonList(right));
        Assert.assertEquals(left.LIBRARY, "lib");
        Assert.assertEquals(left.UNPAIRED_READS_EXAMINED, 16);
        Assert.assertEquals(left.UNPAIRED_WITH_TLEN, 6);
        Assert.assertEquals(left.UNPAIRED_DUPS_WITH_TLEN, 3);
        Assert.assertEquals(left.UNPAIRED_DUPS_WITHOUT_TLEN, 3);
        Assert.assertEquals(left.UNPAIRED_READ_DUPLICATES, 6);
        Assert.assertEquals(left.READ_PAIRS_EXAMINED, 5);
        Assert.assertEquals(left.UNMAPPED_READS, 9);
        // merging a collection derives the rates again
        Assert.assertEquals(left.UNPAIRED_DUP_RATE_WITHOUT_TLEN.doubleValue(), 0.3, 0.0);
        Assert.assertEquals(left.UNPAIRED_DUP_RATE_WITH_TLEN.doubleValue(), 0.5, 0.0);
        Assert.assertEquals(left.PERCENT_DUPLICATION.doubleValue(), 6.0 / 26.0, 0.0);

        // the other side is left alone
        Assert.assertEquals(right.UNPAIRED_WITH_TLEN, 2);
        Assert.assertNull(right.UNPAIRED_DUP_RATE_WITH_TLEN);

        final FlowBasedDuplicationMetrics otherLibrary = flowMetrics(1, 1, 0, 0);
        otherLibrary.LIBRARY = "other";
        Assert.assertThrows(IllegalStateException.class, () -> flowMetrics(1, 1, 0, 0).merge(otherLibrary));
        Assert.assertTrue(left.canMerge(flowMetrics(1, 1, 0, 0)));
        Assert.assertFalse(left.canMerge(new DuplicationMetrics()));
    }
}
