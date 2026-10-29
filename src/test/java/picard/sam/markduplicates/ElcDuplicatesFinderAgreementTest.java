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

import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;
import picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequenceWithBarcodes;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

/**
 * Runs {@link ElcIdenticalBasesDuplicatesFinder}, {@link ElcHashBasedDuplicatesFinder} and
 * {@link ElcDuplicatesFinderResolver} on the same groups of read pairs, where the duplicate sets are known by
 * construction, and checks that both finders produce those sets: exact and near copies, differences spread over
 * both reads, reads of different lengths, MAX_READ_LENGTH, a zero difference rate, order-dependent chains, optical
 * duplicates, barcodes, and a hash-based finder reused across groups. Also covers fillHistogram and minLength.
 */
public class ElcDuplicatesFinderAgreementTest {

    private static final int MIN_IDENTICAL_BASES = 5;
    private static final String PREFIX1 = "ACGTA";
    private static final String PREFIX2 = "TTGCA";

    private static String randomBases(final Random random, final int length) {
        final StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        return sb.toString();
    }

    private static String mutate(final String bases, final int... positions) {
        final char[] chars = bases.toCharArray();
        for (final int p : positions) {
            chars[p] = "CGTA".charAt("ACGT".indexOf(chars[p]));
        }
        return new String(chars);
    }

    private static PairedReadSequence pair(final String read1, final String read2) {
        final PairedReadSequence prs = new PairedReadSequence();
        prs.read1 = read1.getBytes();
        prs.read2 = read2.getBytes();
        return prs;
    }

    /** Pairs that share both 5-base prefixes, as a group handed to a finder would. */
    private static List<PairedReadSequence> templates(final long seed, final int n, final int length) {
        final Random random = new Random(seed);
        final List<PairedReadSequence> out = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            out.add(pair(PREFIX1 + randomBases(random, length - 5), PREFIX2 + randomBases(random, length - 5)));
        }
        return out;
    }

    private static PairedReadSequence copy(final PairedReadSequence p) {
        return pair(new String(p.read1), new String(p.read2));
    }

    private static Map<Integer, Double> bins(final Histogram<Integer> histogram) {
        final Map<Integer, Double> out = new TreeMap<>();
        for (final Integer bin : histogram.keySet()) {
            if (histogram.get(bin).getValue() != 0) {
                out.put(bin, histogram.get(bin).getValue());
            }
        }
        return out;
    }

    private static Map<Integer, Double> histogram(final int... sizeAndCount) {
        final Map<Integer, Double> out = new TreeMap<>();
        for (int i = 0; i < sizeAndCount.length; i += 2) {
            out.put(sizeAndCount[i], (double) sizeAndCount[i + 1]);
        }
        return out;
    }

    private static final class Search {
        final Map<Integer, Double> duplicates;
        final Map<Integer, Double> optical;

        Search(final ElcDuplicatesFinder finder, final List<PairedReadSequence> group) {
            final Histogram<Integer> duplicationHisto = new Histogram<>();
            final Histogram<Integer> opticalHisto = new Histogram<>();
            finder.searchDuplicates(new ArrayList<>(group), duplicationHisto, opticalHisto);
            duplicates = bins(duplicationHisto);
            optical = bins(opticalHisto);
        }
    }

    private static ElcDuplicatesFinder identical(final double rate, final int maxReadLength) {
        return new ElcIdenticalBasesDuplicatesFinder(rate, maxReadLength, MIN_IDENTICAL_BASES, false, new OpticalDuplicateFinder());
    }

    private static ElcDuplicatesFinder hashBased(final double rate, final int maxReadLength) {
        return new ElcHashBasedDuplicatesFinder(rate, maxReadLength, MIN_IDENTICAL_BASES, new OpticalDuplicateFinder());
    }

    @DataProvider(name = "groups")
    public Object[][] groups() {
        final List<Object[]> rows = new ArrayList<>();

        // exact copies of three templates
        {
            final List<PairedReadSequence> t = templates(1, 3, 50);
            rows.add(new Object[]{"exact copies", 0.03, 0,
                    Arrays.asList(t.get(0), t.get(1), copy(t.get(0)), copy(t.get(0)), t.get(2), copy(t.get(1))),
                    histogram(3, 1, 2, 1, 1, 1)});
        }
        // copies one base away from their template are at most two apart from each other
        {
            final List<PairedReadSequence> t = templates(2, 2, 50);
            final PairedReadSequence a = t.get(0);
            final PairedReadSequence b = t.get(1);
            rows.add(new Object[]{"near copies", 0.03, 0,
                    Arrays.asList(pair(mutate(new String(a.read1), 10), new String(a.read2)), b,
                            pair(new String(a.read1), mutate(new String(a.read2), 20)), a,
                            pair(mutate(new String(b.read1), 49), new String(b.read2)),
                            pair(mutate(new String(a.read1), 30), new String(a.read2))),
                    histogram(4, 1, 2, 1)});
        }
        // three differences over both reads are within floor(100 * 0.035); four are not
        {
            final PairedReadSequence a = templates(3, 1, 50).get(0);
            final String r1 = new String(a.read1);
            final String r2 = new String(a.read2);
            rows.add(new Object[]{"split differences", 0.035, 0,
                    Arrays.asList(a, pair(mutate(r1, 10, 12), mutate(r2, 9)), pair(mutate(r1, 30, 31, 32, 33), r2)),
                    histogram(2, 1, 1, 1)});
        }
        // 60-base pairs whose last ten bases differ in five places, and the same pairs cut to 50 bases
        {
            final PairedReadSequence x = templates(4, 1, 60).get(0);
            final PairedReadSequence y = pair(mutate(new String(x.read1), 51, 53, 55, 57, 59), new String(x.read2));
            final PairedReadSequence z = pair(new String(x.read1, 0, 50), new String(x.read2, 0, 50));
            rows.add(new Object[]{"tails, long first", 0.035, 0, Arrays.asList(x, y, z), histogram(2, 1, 1, 1)});
            rows.add(new Object[]{"tails, short first", 0.035, 0, Arrays.asList(copy(z), copy(x), copy(y)), histogram(3, 1)});
        }
        // differences past MAX_READ_LENGTH are not compared
        {
            final PairedReadSequence a = templates(5, 1, 50).get(0);
            final PairedReadSequence b = pair(mutate(new String(a.read1), 35, 40, 45, 49), new String(a.read2));
            rows.add(new Object[]{"past max read length", 0.03, 30, Arrays.asList(a, b), histogram(2, 1)});
            rows.add(new Object[]{"no max read length", 0.03, 0, Arrays.asList(copy(a), copy(b)), histogram(1, 2)});
            rows.add(new Object[]{"max read length at the difference", 0.03, 36, Arrays.asList(copy(a), copy(b)), histogram(2, 1)});
        }
        // a zero rate accepts exact copies only
        {
            final List<PairedReadSequence> t = templates(6, 1, 50);
            final PairedReadSequence a = t.get(0);
            rows.add(new Object[]{"zero rate", 0.0, 0,
                    Arrays.asList(a, pair(mutate(new String(a.read1), 49), new String(a.read2)), copy(a),
                            pair(new String(a.read1), mutate(new String(a.read2), 5))),
                    histogram(2, 1, 1, 2)});
        }
        // a~b and b~c within three differences, a and c four apart: the earliest pair anchors a set
        {
            final PairedReadSequence a = templates(7, 1, 50).get(0);
            final PairedReadSequence b = pair(mutate(new String(a.read1), 10, 20), new String(a.read2));
            final PairedReadSequence c = pair(mutate(new String(b.read1), 30, 40), new String(a.read2));
            rows.add(new Object[]{"chain abc", 0.035, 0, Arrays.asList(a, b, c), histogram(2, 1, 1, 1)});
            rows.add(new Object[]{"chain bac", 0.035, 0, Arrays.asList(copy(b), copy(a), copy(c)), histogram(3, 1)});
            rows.add(new Object[]{"chain cab", 0.035, 0, Arrays.asList(copy(c), copy(a), copy(b)), histogram(2, 1, 1, 1)});
        }
        // a group of 120: 30 templates, each with 3 copies a single base away
        {
            final List<PairedReadSequence> t = templates(8, 30, 50);
            final List<PairedReadSequence> group = new ArrayList<>(t);
            for (int c = 1; c <= 3; c++) {
                for (final PairedReadSequence template : t) {
                    group.add(pair(mutate(new String(template.read1), 5 + 7 * c), new String(template.read2)));
                }
            }
            rows.add(new Object[]{"large group", 0.03, 0, group, histogram(4, 30)});
        }
        // single pair
        rows.add(new Object[]{"single pair", 0.03, 0, templates(9, 1, 50), histogram(1, 1)});

        return rows.toArray(new Object[0][]);
    }

    /**
     * The chain and tails rows expect the first pair in input order to anchor its set. The documentation does not
     * specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "groups")
    public void testFindersFindTheConstructedSets(final String description, final double rate, final int maxReadLength,
                                                  final List<PairedReadSequence> group, final Map<Integer, Double> expected) {
        final Search identical = new Search(identical(rate, maxReadLength), group);
        final Search hashBased = new Search(hashBased(rate, maxReadLength), group);
        Assert.assertEquals(identical.duplicates, expected, "identical-bases finder");
        Assert.assertEquals(hashBased.duplicates, expected, "hash-based finder");
        Assert.assertTrue(identical.optical.isEmpty());
        Assert.assertTrue(hashBased.optical.isEmpty());

        final Histogram<Integer> duplicationHisto = new Histogram<>();
        final Histogram<Integer> opticalHisto = new Histogram<>();
        new ElcDuplicatesFinderResolver(rate, maxReadLength, MIN_IDENTICAL_BASES, false, new OpticalDuplicateFinder())
                .resolveAndSearch(new ArrayList<>(group), duplicationHisto, opticalHisto);
        Assert.assertEquals(bins(duplicationHisto), expected, "resolver");
    }

    @Test
    public void testHashBasedFinderReusedAcrossGroups() {
        // short reads first lower the shared part boundary; long reads after it must still compare their tails
        final PairedReadSequence s = templates(21, 1, 30).get(0);
        final PairedReadSequence x = templates(22, 1, 60).get(0);
        final PairedReadSequence y = pair(mutate(new String(x.read1), 51, 53, 55, 57, 59), new String(x.read2));
        final PairedReadSequence w = pair(mutate(new String(x.read1), 31, 41), new String(x.read2));

        final ElcDuplicatesFinder reused = hashBased(0.035, 0);
        Assert.assertEquals(new Search(reused, Arrays.asList(s, copy(s))).duplicates, histogram(2, 1));
        Assert.assertEquals(new Search(reused, Arrays.asList(x, y, w)).duplicates, histogram(2, 1, 1, 1));
        Assert.assertEquals(new Search(reused, Arrays.asList(copy(y), copy(x), copy(w))).duplicates, histogram(2, 1, 1, 1));
        Assert.assertEquals(new Search(reused, Arrays.asList(copy(s), copy(s), copy(s))).duplicates, histogram(3, 1));

        // the same groups on a fresh finder, long reads first
        final ElcDuplicatesFinder fresh = hashBased(0.035, 0);
        Assert.assertEquals(new Search(fresh, Arrays.asList(copy(x), copy(y), copy(w))).duplicates, histogram(2, 1, 1, 1));
        Assert.assertEquals(new Search(fresh, Arrays.asList(copy(s), copy(s))).duplicates, histogram(2, 1));
    }

    @Test
    public void testOpticalDuplicatesAgree() {
        final PairedReadSequence template = templates(31, 1, 50).get(0);
        final int[][] tileXy = {{1101, 1000, 1000}, {1101, 1100, 1000}, {1101, 1000, 1101}, {1101, 5000, 5000}, {1102, 1000, 1000}};
        final List<PairedReadSequence> group = new ArrayList<>();
        for (final int[] loc : tileXy) {
            final PairedReadSequence p = copy(template);
            p.setReadGroup((short) 0);
            p.setTile((short) loc[0]);
            p.setX(loc[1]);
            p.setY(loc[2]);
            group.add(p);
        }
        // the first pair is kept; only the second lies within 100 pixels of it, nothing else is that close
        for (final ElcDuplicatesFinder finder : Arrays.asList(identical(0.03, 0), hashBased(0.03, 0))) {
            final Search search = new Search(finder, group);
            Assert.assertEquals(search.duplicates, histogram(5, 1));
            Assert.assertEquals(search.optical, histogram(5, 1));
        }
    }

    private static PairedReadSequenceWithBarcodes withBarcodes(final PairedReadSequence p, final int barcode,
                                                               final int readOneBarcode, final int readTwoBarcode) {
        final PairedReadSequenceWithBarcodes out = new PairedReadSequenceWithBarcodes(p);
        out.barcode = barcode;
        out.readOneBarcode = readOneBarcode;
        out.readTwoBarcode = readTwoBarcode;
        return out;
    }

    @DataProvider(name = "barcodeGroups")
    public Object[][] barcodeGroups() {
        return new Object[][]{
                {150, new int[]{1, 2}, new int[]{0}, new int[]{0}, histogram(75, 2)},
                {150, new int[]{1}, new int[]{1, 2, 3}, new int[]{0}, histogram(50, 3)},
                {150, new int[]{1}, new int[]{0}, new int[]{4, 5}, histogram(75, 2)},
                {10, new int[]{1, 2}, new int[]{0}, new int[]{0}, histogram(5, 2)},
        };
    }

    @Test(dataProvider = "barcodeGroups")
    public void testResolverHonoursBarcodesAtAnyGroupSize(final int n, final int[] barcodes, final int[] readOneBarcodes,
                                                          final int[] readTwoBarcodes, final Map<Integer, Double> expected) {
        final PairedReadSequence template = templates(41, 1, 50).get(0);
        final List<PairedReadSequence> group = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            group.add(withBarcodes(template, barcodes[i % barcodes.length], readOneBarcodes[i % readOneBarcodes.length],
                    readTwoBarcodes[i % readTwoBarcodes.length]));
        }
        final Histogram<Integer> withBarcodes = new Histogram<>();
        new ElcDuplicatesFinderResolver(0.03, 0, MIN_IDENTICAL_BASES, true, new OpticalDuplicateFinder())
                .resolveAndSearch(new ArrayList<>(group), withBarcodes, new Histogram<>());
        Assert.assertEquals(bins(withBarcodes), expected);

        final Histogram<Integer> withoutBarcodes = new Histogram<>();
        new ElcDuplicatesFinderResolver(0.03, 0, MIN_IDENTICAL_BASES, false, new OpticalDuplicateFinder())
                .resolveAndSearch(new ArrayList<>(group), withoutBarcodes, new Histogram<>());
        Assert.assertEquals(bins(withoutBarcodes), histogram(n, 1));
    }

    @Test
    public void testFillHistogram() {
        final ElcDuplicatesFinder finder = identical(0.03, 0);
        final PairedReadSequence keeper = templates(51, 1, 50).get(0);
        keeper.setReadGroup((short) 0);
        keeper.setTile((short) 5);
        keeper.setX(100);
        keeper.setY(100);
        final PairedReadSequence near = copy(keeper);
        near.setReadGroup((short) 0);
        near.setTile((short) 5);
        near.setX(150);
        near.setY(150);
        final PairedReadSequence far = copy(keeper);
        far.setReadGroup((short) 0);
        far.setTile((short) 5);
        far.setX(1000);
        far.setY(100);

        final Histogram<Integer> duplicationHisto = new Histogram<>();
        final Histogram<Integer> opticalHisto = new Histogram<>();
        finder.fillHistogram(duplicationHisto, opticalHisto, keeper, new ArrayList<>());
        Assert.assertEquals(bins(duplicationHisto), histogram(1, 1));
        Assert.assertTrue(bins(opticalHisto).isEmpty());

        finder.fillHistogram(duplicationHisto, opticalHisto, keeper, new ArrayList<>(Arrays.asList(near, far)));
        Assert.assertEquals(bins(duplicationHisto), histogram(1, 1, 3, 1));
        Assert.assertEquals(bins(opticalHisto), histogram(3, 1));

        finder.fillHistogram(duplicationHisto, opticalHisto, keeper, new ArrayList<>(Collections.singletonList(far)));
        Assert.assertEquals(bins(duplicationHisto), histogram(1, 1, 2, 1, 3, 1));
        Assert.assertEquals(bins(opticalHisto), histogram(3, 1));
    }

    @DataProvider(name = "minLength")
    public Object[][] minLength() {
        return new Object[][]{
                {0, 50, 60, 50},
                {-1, 60, 50, 50},
                {40, 50, 60, 40},
                {70, 50, 60, 50},
                {50, 50, 60, 50},
        };
    }

    /**
     * A negative maximum read length, like 0, means no maximum (the second row). The documentation does not specify
     * this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "minLength")
    public void testMinLength(final int maxReadLength, final int length1, final int length2, final int expected) {
        for (final ElcDuplicatesFinder finder : Arrays.asList(identical(0.03, maxReadLength), hashBased(0.03, maxReadLength))) {
            Assert.assertEquals(finder.minLength(new byte[length1], new byte[length2]), expected);
        }
    }
}
