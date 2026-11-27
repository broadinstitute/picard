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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Tests {@link DuplicationMetrics}: the derived PERCENT_DUPLICATION and ESTIMATED_LIBRARY_SIZE, the library size
 * estimate against an independent numerical solution of the Lander-Waterman equation, the return-on-investment
 * estimate and histogram, the per-read add methods for every combination of the flags they read, and the debugging
 * main method.
 */
public class DuplicationMetricsBehaviorTest {

    static final SAMFileHeader HEADER = new SAMFileHeader(
            new SAMSequenceDictionary(Collections.singletonList(new SAMSequenceRecord("chr1", 10000))));

    // SAM flag bits
    static final int PAIRED = 0x1;
    static final int PROPER = 0x2;
    static final int UNMAPPED = 0x4;
    static final int MATE_UNMAPPED = 0x8;
    static final int REVERSE = 0x10;
    static final int MATE_REVERSE = 0x20;
    static final int FIRST = 0x40;
    static final int SECOND = 0x80;
    static final int SECONDARY = 0x100;
    static final int QC_FAIL = 0x200;
    static final int DUPLICATE = 0x400;
    static final int SUPPLEMENTARY = 0x800;

    /** Builds a ten-base record with the given flags and cigar at chr1:100. */
    static SAMRecord record(final int flags, final String cigar) {
        final SAMRecord rec = new SAMRecord(HEADER);
        rec.setReadName("read");
        rec.setFlags(flags);
        rec.setReferenceIndex(0);
        rec.setAlignmentStart(100);
        rec.setCigarString(cigar);
        rec.setReadString("ACGTACGTAC");
        rec.setBaseQualityString("IIIIIIIIII");
        if ((flags & PAIRED) != 0) {
            rec.setMateReferenceIndex(0);
            rec.setMateAlignmentStart(300);
        }
        return rec;
    }

    /**
     * Solves C / X = 1 - exp(-N / X) for X by bisection, written independently of the class under test.
     * The bracket starts at X = C, where the left side exceeds the right.
     */
    static double solveLibrarySize(final long pairs, final long uniquePairs) {
        final double n = pairs;
        final double c = uniquePairs;
        double lo = c;
        double hi = c;
        while (excess(hi, n, c) > 0) {
            hi *= 2;
        }
        for (int i = 0; i < 5000; i++) {
            final double mid = lo + (hi - lo) / 2;
            if (mid <= lo || mid >= hi) break;
            if (excess(mid, n, c) > 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return lo + (hi - lo) / 2;
    }

    /** C / X - (1 - exp(-N / X)): positive below the root and negative above it. */
    private static double excess(final double x, final double n, final double c) {
        return c / x + Math.expm1(-n / x);
    }

    /**
     * The integer part of the root, after checking the root is at least 1e-3 from an integer, so that the precision of
     * the estimate cannot decide it. The root always exceeds C, so a root just above C is unambiguous.
     */
    static long expectedLibrarySize(final long pairs, final long uniquePairs) {
        final double root = solveLibrarySize(pairs, uniquePairs);
        final double floor = Math.floor(root);
        final double fraction = root - floor;
        Assert.assertTrue((fraction > 1e-3 || floor == uniquePairs) && fraction < 1 - 1e-3,
                "root " + root + " for " + pairs + ", " + uniquePairs + " is too close to an integer for this test");
        return (long) floor;
    }

    /** Unique pairs expected from x times the sequencing, over the unique pairs observed, by Lander-Waterman. */
    private static double expectedRoi(final long librarySize, final double x, final long pairs, final long uniquePairs) {
        return -librarySize * Math.expm1(-x * pairs / (double) librarySize) / uniquePairs;
    }

    private static DuplicationMetrics metrics(final long unpairedExamined, final long pairsExamined, final long unpairedDuplicates,
                                              final long pairDuplicates, final long opticalPairDuplicates) {
        final DuplicationMetrics metrics = new DuplicationMetrics();
        metrics.LIBRARY = "lib";
        metrics.UNPAIRED_READS_EXAMINED = unpairedExamined;
        metrics.READ_PAIRS_EXAMINED = pairsExamined;
        metrics.UNPAIRED_READ_DUPLICATES = unpairedDuplicates;
        metrics.READ_PAIR_DUPLICATES = pairDuplicates;
        metrics.READ_PAIR_OPTICAL_DUPLICATES = opticalPairDuplicates;
        return metrics;
    }

    @DataProvider(name = "derivedFieldsData")
    public Object[][] derivedFieldsData() {
        return new Object[][]{
                // unpaired examined, pairs examined, unpaired dups, pair dups, optical pair dups, percent duplication,
                // and the (pairs, unique pairs) the library size is estimated from, or null when there is no estimate
                {0L, 0L, 0L, 0L, 0L, 0.0, null},
                {10L, 0L, 3L, 0L, 0L, 0.3, null},
                {7L, 0L, 7L, 0L, 0L, 1.0, null},
                {0L, 100L, 0L, 10L, 0L, 0.1, new long[]{100, 90}},
                {50L, 100L, 5L, 10L, 0L, 0.1, new long[]{100, 90}},
                {1L, 3L, 1L, 1L, 0L, 3.0 / 7.0, new long[]{3, 2}},
                {0L, 1000L, 0L, 200L, 100L, 0.2, new long[]{900, 800}},
                {0L, 1000L, 0L, 200L, 199L, 0.2, new long[]{801, 800}},
                {0L, 1000L, 0L, 100L, 100L, 0.1, null},
                {0L, 1000L, 0L, 0L, 0L, 0.0, null},
                {20L, 1000L, 20L, 999L, 0L, 2018.0 / 2020.0, new long[]{1000, 1}},
        };
    }

    @Test(dataProvider = "derivedFieldsData")
    public void testCalculateDerivedFields(final long unpairedExamined, final long pairsExamined, final long unpairedDuplicates,
                                           final long pairDuplicates, final long opticalPairDuplicates,
                                           final double percentDuplication, final long[] libraryInputs) {
        final DuplicationMetrics metrics = metrics(unpairedExamined, pairsExamined, unpairedDuplicates, pairDuplicates, opticalPairDuplicates);
        Assert.assertNull(metrics.PERCENT_DUPLICATION);
        Assert.assertNull(metrics.ESTIMATED_LIBRARY_SIZE);

        metrics.calculateDerivedFields();
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), percentDuplication, 0.0);
        if (libraryInputs == null) {
            Assert.assertNull(metrics.ESTIMATED_LIBRARY_SIZE);
        } else {
            Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE.longValue(), expectedLibrarySize(libraryInputs[0], libraryInputs[1]));
        }

        // counts are not changed
        Assert.assertEquals(metrics.UNPAIRED_READS_EXAMINED, unpairedExamined);
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, pairsExamined);
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES, pairDuplicates);
        Assert.assertEquals(metrics.READ_PAIR_OPTICAL_DUPLICATES, opticalPairDuplicates);
        Assert.assertEquals(metrics.LIBRARY, "lib");

        // recomputing replaces an earlier value
        metrics.READ_PAIRS_EXAMINED += 1000;
        metrics.UNPAIRED_READ_DUPLICATES = 0;
        metrics.READ_PAIR_DUPLICATES = 0;
        metrics.READ_PAIR_OPTICAL_DUPLICATES = 0;
        metrics.calculateDerivedFields();
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), 0.0, 0.0);
        Assert.assertNull(metrics.ESTIMATED_LIBRARY_SIZE);
    }

    @Test
    @SuppressWarnings("deprecation")
    public void testCalculateDerivedMetricsMatchesDerivedFields() {
        final DuplicationMetrics deprecated = metrics(50, 1000, 5, 200, 100);
        final DuplicationMetrics current = metrics(50, 1000, 5, 200, 100);
        deprecated.calculateDerivedMetrics();
        current.calculateDerivedFields();
        Assert.assertEquals(deprecated.PERCENT_DUPLICATION, current.PERCENT_DUPLICATION);
        Assert.assertEquals(deprecated.ESTIMATED_LIBRARY_SIZE, current.ESTIMATED_LIBRARY_SIZE);
        Assert.assertEquals(deprecated.PERCENT_DUPLICATION.doubleValue(), 405.0 / 2050.0, 0.0);
        Assert.assertEquals(deprecated.ESTIMATED_LIBRARY_SIZE.longValue(), expectedLibrarySize(900, 800));
    }

    @DataProvider(name = "librarySizeData")
    public Object[][] librarySizeData() {
        return new Object[][]{
                {2L, 1L},
                {3L, 2L},
                {10L, 1L},
                {10L, 9L},
                {50L, 49L},
                {100L, 90L},
                {801L, 800L},
                {900L, 800L},
                {1000L, 1L},
                {1000L, 500L},
                {1000L, 900L},
                {1000L, 990L},
                {1000L, 998L},
                {1000L, 999L},
                {10000L, 9000L},
                {10000L, 9990L},
                {12345L, 12295L},
                {100000L, 100L},
                {100000L, 99000L},
                {123456L, 100000L},
                {1000000L, 999000L},
                {5000000L, 4000000L},
        };
    }

    /** The estimate is the integer part of the root of C / X = 1 - exp(-N / X). */
    @Test(dataProvider = "librarySizeData")
    public void testEstimateLibrarySize(final long pairs, final long uniquePairs) {
        Assert.assertEquals(DuplicationMetrics.estimateLibrarySize(pairs, uniquePairs).longValue(), expectedLibrarySize(pairs, uniquePairs));
    }

    @DataProvider(name = "lowDuplicationLibrarySizeData")
    public Object[][] lowDuplicationLibrarySizeData() {
        return new Object[][]{
                {100000L, 99990L},
                {1000000L, 999990L},
                {1000000L, 999900L},
                {10000000L, 9999900L},
                {10000000L, 9999000L},
                {100000000L, 99999000L},
        };
    }

    /** With few duplicates the root is large and only its leading digits are well determined; agree to 1e-6. */
    @Test(dataProvider = "lowDuplicationLibrarySizeData")
    public void testEstimateLibrarySizeLowDuplication(final long pairs, final long uniquePairs) {
        final double expected = solveLibrarySize(pairs, uniquePairs);
        final long actual = DuplicationMetrics.estimateLibrarySize(pairs, uniquePairs);
        Assert.assertEquals(actual, expected, expected * 1e-6);
        // the library is larger than the number of distinct pairs seen, and grows as duplicates get rarer
        Assert.assertTrue(actual > pairs);
        Assert.assertTrue(DuplicationMetrics.estimateLibrarySize(pairs, uniquePairs + 1) > actual);
    }

    @DataProvider(name = "pinnedLibrarySizeData")
    public Object[][] pinnedLibrarySizeData() {
        return new Object[][]{
                {10000000L, 9800671L, 247497047L},
                {1000000000L, 900000098L, 4660798385L},
                {1000000000L, 799999304L, 2154175923L},
        };
    }

    /**
     * Pins the current output for inputs whose root lies close enough to an integer that the precision of the
     * bisection decides the result; each value is within one of the integer part of the root.
     */
    @Test(dataProvider = "pinnedLibrarySizeData")
    public void testEstimateLibrarySizePinned(final long pairs, final long uniquePairs, final long pinned) {
        Assert.assertEquals(DuplicationMetrics.estimateLibrarySize(pairs, uniquePairs).longValue(), pinned);
        Assert.assertEquals(pinned, Math.floor(solveLibrarySize(pairs, uniquePairs)), 1.0);
    }

    @DataProvider(name = "noEstimateData")
    public Object[][] noEstimateData() {
        return new Object[][]{
                {0L, 0L},
                {1L, 1L},
                {100L, 100L},
                {Long.MAX_VALUE, Long.MAX_VALUE},
        };
    }

    /** No duplicates, or no pairs, give no estimate. */
    @Test(dataProvider = "noEstimateData")
    public void testEstimateLibrarySizeWithoutDuplicates(final long pairs, final long uniquePairs) {
        Assert.assertNull(DuplicationMetrics.estimateLibrarySize(pairs, uniquePairs));
    }

    /** When every pair after the first is a duplicate the library holds a single molecule. */
    @Test
    public void testEstimateLibrarySizeAllDuplicates() {
        for (final long pairs : new long[]{2, 3, 10, 1000, 1000000, Integer.MAX_VALUE, 1L << 40}) {
            Assert.assertEquals(DuplicationMetrics.estimateLibrarySize(pairs, 1).longValue(), 1L, "pairs = " + pairs);
        }
    }

    @Test
    public void testEstimateLibrarySizeGrowsWithUniquePairs() {
        long previous = 0;
        for (long unique = 1; unique < 1000; unique++) {
            final long estimate = DuplicationMetrics.estimateLibrarySize(1000, unique);
            Assert.assertTrue(estimate >= unique, "unique = " + unique);
            Assert.assertTrue(estimate >= previous, "unique = " + unique);
            previous = estimate;
        }
    }

    @Test
    public void testEstimateRoi() {
        Assert.assertEquals(DuplicationMetrics.estimateRoi(1000, 0, 10, 5), 0.0, 0.0);
        // sequencing far deeper sees every molecule: library size over unique pairs
        Assert.assertEquals(DuplicationMetrics.estimateRoi(1000, 1e7, 10, 5), 200.0, 0.0);
        Assert.assertEquals(DuplicationMetrics.estimateRoi(466, 1e6, 100, 90), 466.0 / 90.0, 1e-12);

        for (final double x : new double[]{0.5, 1, 2, 3.5, 10, 100}) {
            Assert.assertEquals(DuplicationMetrics.estimateRoi(466, x, 100, 90), expectedRoi(466, x, 100, 90), 1e-12, "x = " + x);
            Assert.assertEquals(DuplicationMetrics.estimateRoi(5000, x, 1000, 800), expectedRoi(5000, x, 1000, 800), 1e-12, "x = " + x);
            Assert.assertEquals(DuplicationMetrics.estimateRoi(1L << 30, x, 1000, 1000), expectedRoi(1L << 30, x, 1000, 1000), 1e-8, "x = " + x);
        }

        // at x = 1 with the library size the counts imply, the yield is what was observed
        Assert.assertEquals(DuplicationMetrics.estimateRoi(expectedLibrarySize(1000000, 999000), 1, 1000000, 999000), 1.0, 1e-8);
    }

    @Test
    public void testRoiHistogram() {
        final DuplicationMetrics metrics = metrics(0, 100, 0, 10, 0);
        final Histogram<Double> histogram = metrics.calculateRoiHistogram();

        // the derived fields are filled in on the way
        Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE.longValue(), expectedLibrarySize(100, 90));
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), 0.1, 0.0);

        assertRoiHistogram(histogram, metrics.ESTIMATED_LIBRARY_SIZE, 100, 90);
        Assert.assertEquals(histogram.get(1.0).getValue(), 1.0, 1e-3);
    }

    /** The histogram uses all pairs examined, while the library size excludes optical duplicates. */
    @Test
    public void testRoiHistogramWithOpticalDuplicates() {
        final DuplicationMetrics metrics = metrics(30, 1000, 3, 200, 100);
        final Histogram<Double> histogram = metrics.calculateRoiHistogram();
        Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE.longValue(), expectedLibrarySize(900, 800));
        assertRoiHistogram(histogram, metrics.ESTIMATED_LIBRARY_SIZE, 1000, 800);
    }

    /** A library size already present is used as is, and the derived fields are not recomputed. */
    @Test
    public void testRoiHistogramUsesExistingLibrarySize() {
        final DuplicationMetrics metrics = metrics(0, 100, 0, 10, 0);
        metrics.ESTIMATED_LIBRARY_SIZE = 1000L;
        final Histogram<Double> histogram = metrics.calculateRoiHistogram();
        Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE.longValue(), 1000L);
        Assert.assertNull(metrics.PERCENT_DUPLICATION);
        assertRoiHistogram(histogram, 1000, 100, 90);
    }

    @Test
    public void testRoiHistogramWithoutEstimate() {
        final DuplicationMetrics noDuplicates = metrics(10, 100, 2, 0, 0);
        Assert.assertNull(noDuplicates.calculateRoiHistogram());
        Assert.assertNull(noDuplicates.ESTIMATED_LIBRARY_SIZE);
        Assert.assertEquals(noDuplicates.PERCENT_DUPLICATION.doubleValue(), 2.0 / 210.0, 0.0);

        Assert.assertNull(metrics(10, 0, 5, 0, 0).calculateRoiHistogram());
        Assert.assertNull(metrics(0, 0, 0, 0, 0).calculateRoiHistogram());
        Assert.assertNull(metrics(0, 100, 0, 10, 10).calculateRoiHistogram());
    }

    private static void assertRoiHistogram(final Histogram<Double> histogram, final long librarySize, final long pairs, final long uniquePairs) {
        Assert.assertNotNull(histogram);
        Assert.assertEquals(histogram.getValueLabel(), "CoverageMult");
        Assert.assertEquals(histogram.size(), 100);
        final List<Double> keys = new ArrayList<>(histogram.keySet());
        double previous = 0;
        for (int i = 0; i < 100; i++) {
            final double x = i + 1;
            Assert.assertEquals(keys.get(i), x);
            final double value = histogram.get(x).getValue();
            Assert.assertEquals(value, expectedRoi(librarySize, x, pairs, uniquePairs), 1e-12, "x = " + x);
            Assert.assertTrue(value > previous, "x = " + x);
            previous = value;
        }
    }

    // indexes into counts()
    private static final int UNPAIRED_EXAMINED = 0;
    private static final int PAIRS_EXAMINED = 1;
    private static final int SECONDARY_OR_SUPPLEMENTARY = 2;
    private static final int UNMAPPED_READS = 3;
    private static final int UNPAIRED_DUPLICATES = 4;
    private static final int PAIR_DUPLICATES = 5;
    private static final int OPTICAL_DUPLICATES = 6;
    private static final int NONE = -1;

    static long[] counts(final DuplicationMetrics metrics) {
        return new long[]{
                metrics.UNPAIRED_READS_EXAMINED,
                metrics.READ_PAIRS_EXAMINED,
                metrics.SECONDARY_OR_SUPPLEMENTARY_RDS,
                metrics.UNMAPPED_READS,
                metrics.UNPAIRED_READ_DUPLICATES,
                metrics.READ_PAIR_DUPLICATES,
                metrics.READ_PAIR_OPTICAL_DUPLICATES,
        };
    }

    private static long[] oneAt(final int... indexes) {
        final long[] expected = new long[7];
        for (final int index : indexes) {
            if (index != NONE) expected[index]++;
        }
        return expected;
    }

    @DataProvider(name = "readFlagsData")
    public Object[][] readFlagsData() {
        return new Object[][]{
                // flags, counter for addReadToLibraryMetrics, counter for addDuplicateReadToMetrics
                {0, UNPAIRED_EXAMINED, UNPAIRED_DUPLICATES},
                {REVERSE, UNPAIRED_EXAMINED, UNPAIRED_DUPLICATES},
                {DUPLICATE, UNPAIRED_EXAMINED, UNPAIRED_DUPLICATES},
                {QC_FAIL, UNPAIRED_EXAMINED, UNPAIRED_DUPLICATES},
                {PAIRED | PROPER | FIRST | MATE_REVERSE, PAIRS_EXAMINED, PAIR_DUPLICATES},
                {PAIRED | SECOND | REVERSE, PAIRS_EXAMINED, PAIR_DUPLICATES},
                {PAIRED | FIRST | DUPLICATE, PAIRS_EXAMINED, PAIR_DUPLICATES},
                {PAIRED | MATE_UNMAPPED | FIRST, UNPAIRED_EXAMINED, UNPAIRED_DUPLICATES},
                {PAIRED | MATE_UNMAPPED | SECOND | REVERSE, UNPAIRED_EXAMINED, UNPAIRED_DUPLICATES},
                {UNMAPPED, UNMAPPED_READS, NONE},
                {UNMAPPED | PAIRED | MATE_UNMAPPED | FIRST, UNMAPPED_READS, NONE},
                {UNMAPPED | PAIRED | SECOND, UNMAPPED_READS, NONE},
                {UNMAPPED | SECONDARY, UNMAPPED_READS, NONE},
                {UNMAPPED | DUPLICATE, UNMAPPED_READS, NONE},
                {SECONDARY, SECONDARY_OR_SUPPLEMENTARY, NONE},
                {SUPPLEMENTARY, SECONDARY_OR_SUPPLEMENTARY, NONE},
                {SECONDARY | SUPPLEMENTARY, SECONDARY_OR_SUPPLEMENTARY, NONE},
                {PAIRED | FIRST | SECONDARY, SECONDARY_OR_SUPPLEMENTARY, NONE},
                {PAIRED | SECOND | SUPPLEMENTARY | DUPLICATE, SECONDARY_OR_SUPPLEMENTARY, NONE},
                {PAIRED | MATE_UNMAPPED | FIRST | SUPPLEMENTARY, SECONDARY_OR_SUPPLEMENTARY, NONE},
        };
    }

    @Test(dataProvider = "readFlagsData")
    public void testAddReadMethods(final int flags, final int libraryCounter, final int duplicateCounter) {
        for (final DuplicationMetrics metrics : new DuplicationMetrics[]{new DuplicationMetrics(), new FlowBasedDuplicationMetrics()}) {
            metrics.LIBRARY = "lib";
            final SAMRecord rec = record(flags, "10M");

            metrics.addReadToLibraryMetrics(rec);
            Assert.assertEquals(counts(metrics), oneAt(libraryCounter));

            metrics.addDuplicateReadToMetrics(rec);
            Assert.assertEquals(counts(metrics), oneAt(libraryCounter, duplicateCounter));

            metrics.addReadToLibraryMetrics(rec);
            metrics.addDuplicateReadToMetrics(rec);
            final long[] twice = oneAt(libraryCounter, duplicateCounter);
            for (int i = 0; i < twice.length; i++) twice[i] *= 2;
            Assert.assertEquals(counts(metrics), twice);

            Assert.assertEquals(metrics.LIBRARY, "lib");
            Assert.assertNull(metrics.PERCENT_DUPLICATION);
            Assert.assertNull(metrics.ESTIMATED_LIBRARY_SIZE);
            // the record is not modified
            Assert.assertEquals(rec.getFlags(), flags);
        }
    }

    /** Reads are counted one at a time: both ends of a pair add to READ_PAIRS_EXAMINED and READ_PAIR_DUPLICATES. */
    @Test
    public void testAddReadsThenDerive() {
        final DuplicationMetrics metrics = new DuplicationMetrics();
        final List<SAMRecord> examined = new ArrayList<>();
        final List<SAMRecord> duplicates = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            final SAMRecord first = record(PAIRED | FIRST, "10M");
            final SAMRecord second = record(PAIRED | SECOND | REVERSE, "10M");
            examined.add(first);
            examined.add(second);
            if (i < 3) {
                duplicates.add(first);
                duplicates.add(second);
            }
        }
        for (int i = 0; i < 4; i++) {
            final SAMRecord fragment = record(0, "10M");
            examined.add(fragment);
            if (i == 0) duplicates.add(fragment);
        }
        examined.add(record(UNMAPPED, "10M"));
        examined.add(record(SECONDARY, "10M"));
        duplicates.add(record(SECONDARY, "10M"));

        examined.forEach(metrics::addReadToLibraryMetrics);
        duplicates.forEach(metrics::addDuplicateReadToMetrics);
        Assert.assertEquals(counts(metrics), new long[]{4, 20, 1, 1, 1, 6, 0});

        // the callers halve the pair counts before deriving
        metrics.READ_PAIRS_EXAMINED /= 2;
        metrics.READ_PAIR_DUPLICATES /= 2;
        metrics.calculateDerivedFields();
        Assert.assertEquals(metrics.PERCENT_DUPLICATION.doubleValue(), 7.0 / 24.0, 0.0);
        Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE.longValue(), expectedLibrarySize(10, 7));
    }

    @Test
    public void testMain() {
        final PrintStream originalOut = System.out;
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        try (final PrintStream capture = new PrintStream(bytes, true, StandardCharsets.UTF_8)) {
            System.setOut(capture);
            DuplicationMetrics.main(new String[]{"100", "10"});
        } finally {
            System.setOut(originalOut);
        }
        final String[] lines = bytes.toString(StandardCharsets.UTF_8).split("\\R");
        final long librarySize = expectedLibrarySize(100, 90);

        Assert.assertEquals(lines.length, 104);
        Assert.assertEquals(lines[0], "Percent Duplication: 0.1");
        Assert.assertEquals(lines[1], "Est. Library Size  : " + librarySize);
        Assert.assertEquals(lines[2], "");
        Assert.assertEquals(lines[3], "X Seq\tX Unique");
        for (int i = 0; i < 100; i++) {
            final String[] fields = lines[4 + i].split("\t");
            Assert.assertEquals(fields.length, 2, lines[4 + i]);
            Assert.assertEquals(fields[0], (i + 1) + ".0");
            Assert.assertEquals(Double.parseDouble(fields[1]), expectedRoi(librarySize, i + 1, 100, 90), 1e-12);
        }
    }
}
