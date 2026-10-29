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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.SequenceUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

/**
 * End-to-end tests of {@link EstimateLibraryComplexity} on small paired inputs whose duplicate sets are known by
 * construction: the duplicate-set histogram and per-library metrics, MIN_IDENTICAL_BASES, MAX_DIFF_RATE,
 * MIN_MEAN_QUALITY, MAX_GROUP_RATIO, MIN_GROUP_COUNT, MAX_READ_LENGTH, the barcode tags, optical duplicates,
 * libraries and read groups, multiple inputs, aligned records, and spilling the sort to disk.
 */
public class EstimateLibraryComplexityBehaviorTest extends CommandLineProgramTest {

    private static final int LENGTH = 50;
    private static final byte QUAL = 30;
    private static final String RG_A = "rgA";

    @Override
    public String getCommandLineProgramName() {
        return EstimateLibraryComplexity.class.getSimpleName();
    }

    /** One read pair in sequencing orientation, with everything needed to write it as SAM records. */
    private static final class Pair {
        String read1;
        String read2;
        byte[] quals1;
        byte[] quals2;
        String readGroup = RG_A;
        int tile;
        int x = 1000;
        int y = 1000;
        final Map<String, String> tags1 = new HashMap<>();
        final Map<String, String> tags2 = new HashMap<>();

        Pair(final String read1, final String read2) {
            this.read1 = read1;
            this.read2 = read2;
            this.quals1 = quals(read1.length(), QUAL);
            this.quals2 = quals(read2.length(), QUAL);
        }

        Pair copy() {
            final Pair p = new Pair(read1, read2);
            p.quals1 = quals1.clone();
            p.quals2 = quals2.clone();
            p.readGroup = readGroup;
            p.tags1.putAll(tags1);
            p.tags2.putAll(tags2);
            return p;
        }

        Pair at(final int tile, final int x, final int y) {
            this.tile = tile;
            this.x = x;
            this.y = y;
            return this;
        }
    }

    private static byte[] quals(final int length, final byte q) {
        final byte[] quals = new byte[length];
        Arrays.fill(quals, q);
        return quals;
    }

    private static String randomBases(final Random random, final int length) {
        final StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append("ACGT".charAt(random.nextInt(4)));
        }
        return sb.toString();
    }

    /** Substitutes the base at each position with the next base in ACGT order. */
    private static String mutate(final String bases, final int... positions) {
        final char[] chars = bases.toCharArray();
        for (final int p : positions) {
            chars[p] = "CGTA".charAt("ACGT".indexOf(chars[p]));
        }
        return new String(chars);
    }

    private static List<Pair> randomPairs(final long seed, final int n) {
        final Random random = new Random(seed);
        final List<Pair> pairs = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            pairs.add(new Pair(randomBases(random, LENGTH), randomBases(random, LENGTH)));
        }
        return pairs;
    }

    private static SAMFileHeader header(final String... readGroupIdAndLibrary) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("chr1", 100_000));
        for (int i = 0; i < readGroupIdAndLibrary.length; i += 2) {
            final SAMReadGroupRecord rg = new SAMReadGroupRecord(readGroupIdAndLibrary[i]);
            rg.setSample("sample");
            if (readGroupIdAndLibrary[i + 1] != null) {
                rg.setLibrary(readGroupIdAndLibrary[i + 1]);
            }
            header.addReadGroup(rg);
        }
        return header;
    }

    private static SAMFileHeader defaultHeader() {
        return header(RG_A, "libA");
    }

    private static SAMRecord unalignedMate(final SAMFileHeader header, final String name, final boolean first,
                                           final String bases, final byte[] quals, final String readGroup,
                                           final Map<String, String> tags) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName(name);
        rec.setReadPairedFlag(true);
        rec.setFirstOfPairFlag(first);
        rec.setSecondOfPairFlag(!first);
        rec.setReadUnmappedFlag(true);
        rec.setMateUnmappedFlag(true);
        rec.setReadBases(bases.getBytes());
        rec.setBaseQualities(quals.clone());
        if (readGroup != null) {
            rec.setAttribute("RG", readGroup);
        }
        tags.forEach(rec::setAttribute);
        return rec;
    }

    /** Builds unaligned mate records in pair order; read names carry tile/x/y in the default Illumina layout. */
    private static List<SAMRecord> toRecords(final SAMFileHeader header, final List<Pair> pairs, final String namePrefix) {
        final List<SAMRecord> records = new ArrayList<>();
        int index = 0;
        for (final Pair p : pairs) {
            final int tile = p.tile == 0 ? 1 + index : p.tile;
            final String name = namePrefix + index + ":" + tile + ":" + p.x + ":" + p.y;
            records.add(unalignedMate(header, name, true, p.read1, p.quals1, p.readGroup, p.tags1));
            records.add(unalignedMate(header, name, false, p.read2, p.quals2, p.readGroup, p.tags2));
            index++;
        }
        return records;
    }

    private static File writeSam(final SAMFileHeader header, final List<SAMRecord> records) throws IOException {
        final File file = File.createTempFile("elcBehavior.", ".sam");
        file.deleteOnExit();
        try (SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(header, true, file)) {
            records.forEach(writer::addAlignment);
        }
        return file;
    }

    private static File writePairs(final SAMFileHeader header, final List<Pair> pairs) throws IOException {
        return writeSam(header, toRecords(header, pairs, "RUN:"));
    }

    private static final class ElcResult {
        final Map<String, DuplicationMetrics> metricsByLibrary = new LinkedHashMap<>();
        final Map<String, Map<Integer, Double>> histogramByLibrary = new HashMap<>();
    }

    private ElcResult run(final List<File> inputs, final String... extraArgs) throws IOException {
        final File output = File.createTempFile("elcBehavior.", ".elc_metrics");
        output.deleteOnExit();
        final List<String> args = new ArrayList<>();
        for (final File input : inputs) {
            args.add("INPUT=" + input.getAbsolutePath());
        }
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.addAll(Arrays.asList(extraArgs));
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final ElcResult result = new ElcResult();
        final MetricsFile<DuplicationMetrics, Integer> metricsFile = new MetricsFile<>();
        try (Reader reader = new FileReader(output)) {
            metricsFile.read(reader);
        }
        for (final DuplicationMetrics m : metricsFile.getMetrics()) {
            Assert.assertNull(result.metricsByLibrary.put(m.LIBRARY, m), "library reported twice: " + m.LIBRARY);
        }
        for (final Histogram<Integer> histogram : metricsFile.getAllHistograms()) {
            final Map<Integer, Double> bins = new TreeMap<>();
            for (final Integer bin : histogram.keySet()) {
                final double value = histogram.get(bin).getValue();
                if (value != 0) {
                    bins.put(bin, value);
                }
            }
            result.histogramByLibrary.put(histogram.getValueLabel(), bins);
        }
        return result;
    }

    private ElcResult run(final File input, final String... extraArgs) throws IOException {
        return run(Collections.singletonList(input), extraArgs);
    }

    /** Histogram expectation from alternating (set size, number of sets) values. */
    private static Map<Integer, Double> histogram(final int... sizeAndCount) {
        final Map<Integer, Double> bins = new TreeMap<>();
        for (int i = 0; i < sizeAndCount.length; i += 2) {
            bins.put(sizeAndCount[i], (double) sizeAndCount[i + 1]);
        }
        return bins;
    }

    /** Solves the Lander-Waterman relation c/x = 1 - exp(-n/x) for x by bisection; null when nothing is duplicated. */
    private static Long referenceLibrarySize(final long pairs, final long uniquePairs) {
        if (pairs <= 0 || pairs - uniquePairs <= 0) {
            return null;
        }
        double lo = uniquePairs;
        double hi = uniquePairs;
        while (uniquePairs / hi - 1 + Math.exp(-pairs / hi) > 0) {
            hi *= 2;
        }
        for (int i = 0; i < 200; i++) {
            final double mid = (lo + hi) / 2;
            if (uniquePairs / mid - 1 + Math.exp(-pairs / mid) > 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return (long) lo;
    }

    private static void assertMetrics(final DuplicationMetrics m, final String library, final long examined,
                                      final long duplicates, final long optical) {
        Assert.assertNotNull(m, "no metrics for library " + library);
        Assert.assertEquals(m.LIBRARY, library);
        Assert.assertEquals(m.READ_PAIRS_EXAMINED, examined, "READ_PAIRS_EXAMINED");
        Assert.assertEquals(m.READ_PAIR_DUPLICATES, duplicates, "READ_PAIR_DUPLICATES");
        Assert.assertEquals(m.READ_PAIR_OPTICAL_DUPLICATES, optical, "READ_PAIR_OPTICAL_DUPLICATES");
        Assert.assertEquals(m.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(m.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(m.UNMAPPED_READS, 0);
        Assert.assertEquals(m.SECONDARY_OR_SUPPLEMENTARY_RDS, 0);
        Assert.assertEquals(m.PERCENT_DUPLICATION, examined == 0 ? 0.0 : duplicates / (double) examined, 1e-6);
        final Long expectedSize = referenceLibrarySize(examined - optical, examined - duplicates);
        if (expectedSize == null) {
            Assert.assertNull(m.ESTIMATED_LIBRARY_SIZE);
        } else {
            Assert.assertNotNull(m.ESTIMATED_LIBRARY_SIZE);
            Assert.assertEquals((double) m.ESTIMATED_LIBRARY_SIZE, (double) expectedSize, 1.0);
        }
    }

    /** Sets of 3 and 2 exact copies and 4 unrelated pairs, interleaved. */
    private static List<Pair> mixedSets() {
        final List<Pair> t = randomPairs(11, 6);
        return Arrays.asList(t.get(0).copy(), t.get(1).copy(), t.get(2).copy(), t.get(0).copy(), t.get(3).copy(),
                t.get(1).copy(), t.get(4).copy(), t.get(0).copy(), t.get(5).copy());
    }

    @DataProvider(name = "minGroupCount")
    public Object[][] minGroupCount() {
        return new Object[][]{
                // bins: one set of 3, one set of 2, four singletons
                {"MIN_GROUP_COUNT=1", 9, 3},
                {null, 4, 0},
                {"MIN_GROUP_COUNT=4", 4, 0},
                {"MIN_GROUP_COUNT=5", 0, 0},
        };
    }

    @Test(dataProvider = "minGroupCount")
    public void testDuplicateSetsAndMinGroupCount(final String arg, final long examined, final long duplicates) throws IOException {
        final File input = writePairs(defaultHeader(), mixedSets());
        final ElcResult result = arg == null ? run(input) : run(input, arg);

        Assert.assertEquals(result.metricsByLibrary.keySet(), Collections.singleton("libA"));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", examined, duplicates, 0);
        // the histogram is written unfiltered
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(1, 4, 2, 1, 3, 1));
    }

    @Test
    public void testSpillingTheSortGivesTheSameResult() throws IOException {
        final File input = writePairs(defaultHeader(), mixedSets());
        final ElcResult result = run(input, "MIN_GROUP_COUNT=1", "MAX_RECORDS_IN_RAM=2");
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 9, 3, 0);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(1, 4, 2, 1, 3, 1));
    }

    @Test
    public void testNoPairsGivesNoMetrics() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final SAMRecord unpaired = new SAMRecord(header);
        unpaired.setReadName("RUN:0:1:1:1");
        unpaired.setReadUnmappedFlag(true);
        unpaired.setReadBases(randomBases(new Random(3), LENGTH).getBytes());
        unpaired.setBaseQualities(quals(LENGTH, QUAL));
        unpaired.setAttribute("RG", RG_A);

        for (final List<SAMRecord> records : Arrays.asList(Collections.<SAMRecord>emptyList(), Collections.singletonList(unpaired))) {
            final ElcResult result = run(writeSam(header, records), "MIN_GROUP_COUNT=1");
            Assert.assertTrue(result.metricsByLibrary.isEmpty());
            Assert.assertTrue(result.histogramByLibrary.isEmpty());
        }
    }

    @Test
    public void testSinglePair() throws IOException {
        final ElcResult result = run(writePairs(defaultHeader(), randomPairs(5, 1)), "MIN_GROUP_COUNT=1");
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 1, 0, 0);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(1, 1));
    }

    @DataProvider(name = "maxDiffRate")
    public Object[][] maxDiffRate() {
        // 100 bases per pair; mismatches are placed after the first 5 bases
        return new Object[][]{
                {0.035, new int[]{10, 20}, new int[]{15}, true},
                {0.035, new int[]{10, 20, 30}, new int[]{}, true},
                {0.035, new int[]{}, new int[]{15, 25, 35}, true},
                {0.035, new int[]{10, 20}, new int[]{15, 25}, false},
                {0.035, new int[]{10, 20, 30, 40}, new int[]{}, false},
                {0.045, new int[]{10, 20}, new int[]{15, 25}, true},
                {0.045, new int[]{10, 20, 30}, new int[]{15, 25}, false},
                {0.0, new int[]{}, new int[]{}, true},
                {0.0, new int[]{49}, new int[]{}, false},
                {0.0, new int[]{}, new int[]{5}, false},
        };
    }

    @Test(dataProvider = "maxDiffRate")
    public void testMaxDiffRate(final double rate, final int[] read1Mismatches, final int[] read2Mismatches,
                                final boolean duplicate) throws IOException {
        final Pair p = randomPairs(21, 1).get(0);
        final Pair q = new Pair(mutate(p.read1, read1Mismatches), mutate(p.read2, read2Mismatches));
        final ElcResult result = run(writePairs(defaultHeader(), Arrays.asList(p, q)),
                "MIN_GROUP_COUNT=1", "MAX_DIFF_RATE=" + rate);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), duplicate ? histogram(2, 1) : histogram(1, 2));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 2, duplicate ? 1 : 0, 0);
    }

    @Test
    public void testDiffRateUsesTheShorterReadOfEachEnd() throws IOException {
        // read 1 of q and r is 40 bases, so 90 bases are compared with p and floor(90 * 0.042) = 3 mismatches are
        // allowed; with 100 bases the limit would be 4. q differs from p in 3 places, r in 4 and from q in 7.
        final Pair p = randomPairs(23, 1).get(0);
        final Pair q = new Pair(mutate(p.read1, 10, 20, 30).substring(0, 40), p.read2);
        final Pair r = new Pair(mutate(p.read1, 11, 21, 31, 36).substring(0, 40), p.read2);
        final ElcResult result = run(writePairs(defaultHeader(), Arrays.asList(p, q, r)),
                "MIN_GROUP_COUNT=1", "MAX_DIFF_RATE=0.042");
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(2, 1, 1, 1));
    }

    @DataProvider(name = "minIdenticalBases")
    public Object[][] minIdenticalBases() {
        return new Object[][]{
                // a difference inside the prefix splits pairs into different groups, which are never compared
                {5, true, 2, false},
                {2, true, 2, true},
                {3, true, 2, false},
                {5, false, 4, false},
                {4, false, 4, true},
                {1, false, 0, false},
                {1, false, 1, true},
        };
    }

    @Test(dataProvider = "minIdenticalBases")
    public void testMinIdenticalBases(final int minIdenticalBases, final boolean inRead1, final int position,
                                      final boolean duplicate) throws IOException {
        final Pair p = randomPairs(31, 1).get(0);
        final Pair q = inRead1 ? new Pair(mutate(p.read1, position), p.read2) : new Pair(p.read1, mutate(p.read2, position));
        final ElcResult result = run(writePairs(defaultHeader(), Arrays.asList(p, q)),
                "MIN_GROUP_COUNT=1", "MIN_IDENTICAL_BASES=" + minIdenticalBases, "MAX_DIFF_RATE=0.035");
        Assert.assertEquals(result.histogramByLibrary.get("libA"), duplicate ? histogram(2, 1) : histogram(1, 2));
    }

    @DataProvider(name = "maxReadLength")
    public Object[][] maxReadLength() {
        // read 1 of q differs from p at every base from 30 on
        return new Object[][]{
                {0, false},
                {30, true},
                // 2 compared mismatches of at most floor(64 * 0.035) = 2
                {32, true},
                // 3 compared mismatches of at most floor(66 * 0.035) = 2
                {33, false},
        };
    }

    @Test(dataProvider = "maxReadLength")
    public void testMaxReadLengthLimitsComparedBases(final int maxReadLength, final boolean duplicate) throws IOException {
        final Pair p = randomPairs(41, 1).get(0);
        final int[] tail = new int[LENGTH - 30];
        for (int i = 0; i < tail.length; i++) {
            tail[i] = 30 + i;
        }
        final Pair q = new Pair(mutate(p.read1, tail), p.read2);
        final ElcResult result = run(writePairs(defaultHeader(), Arrays.asList(p, q)),
                "MIN_GROUP_COUNT=1", "MAX_DIFF_RATE=0.035", "MAX_READ_LENGTH=" + maxReadLength);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), duplicate ? histogram(2, 1) : histogram(1, 2));
    }

    /**
     * MAX_READ_LENGTH also limits the bases whose mean quality is compared with MIN_MEAN_QUALITY. The documentation
     * does not specify this; the expected value follows the current implementation.
     */
    @Test
    public void testMaxReadLengthLimitsTheQualityWindow() throws IOException {
        // mean quality over 50 bases is (30 * 30 + 20 * 2) / 50 = 18.8, over the first 30 bases it is 30
        final Pair p = randomPairs(43, 1).get(0);
        final Pair q = p.copy();
        Arrays.fill(q.quals1, 30, LENGTH, (byte) 2);

        final File input = writePairs(defaultHeader(), Arrays.asList(p, q));
        final ElcResult all = run(input, "MIN_GROUP_COUNT=1");
        Assert.assertEquals(all.histogramByLibrary.get("libA"), histogram(1, 1));
        assertMetrics(all.metricsByLibrary.get("libA"), "libA", 1, 0, 0);

        final ElcResult limited = run(input, "MIN_GROUP_COUNT=1", "MAX_READ_LENGTH=30");
        Assert.assertEquals(limited.histogramByLibrary.get("libA"), histogram(2, 1));
        assertMetrics(limited.metricsByLibrary.get("libA"), "libA", 2, 1, 0);
    }

    @DataProvider(name = "qualityFilter")
    public Object[][] qualityFilter() {
        return new Object[][]{
                // {which read of q is degraded, how, MIN_MEAN_QUALITY, q kept}
                {2, "all20", 20, true},
                {2, "mean19.98", 20, false},
                {1, "mean19.98", 20, false},
                {1, "mean19.98", 19, true},
                {1, "all0", 0, true},
                {1, "all0", 1, false},
                {1, "noCallInPrefix", 20, false},
                {2, "noCallInPrefix", 20, false},
                // a no-call after the prefix is kept and counts as one mismatch
                {1, "noCallAfterPrefix", 20, true},
                {1, "shorterThanPrefix", 20, false},
        };
    }

    @Test(dataProvider = "qualityFilter")
    public void testQualityFilter(final int read, final String how, final int minMeanQuality, final boolean kept) throws IOException {
        final Pair p = randomPairs(51, 1).get(0);
        final Pair q = p.copy();
        String bases = read == 1 ? q.read1 : q.read2;
        byte[] quals = read == 1 ? q.quals1 : q.quals2;
        switch (how) {
            case "all20":
                Arrays.fill(quals, (byte) 20);
                break;
            case "mean19.98":
                Arrays.fill(quals, (byte) 20);
                quals[LENGTH - 1] = 19;
                break;
            case "all0":
                Arrays.fill(quals, (byte) 0);
                break;
            case "noCallInPrefix":
                bases = bases.substring(0, 4) + "N" + bases.substring(5);
                break;
            case "noCallAfterPrefix":
                bases = bases.substring(0, 5) + "N" + bases.substring(6);
                break;
            case "shorterThanPrefix":
                bases = bases.substring(0, 4);
                quals = quals(4, QUAL);
                break;
            default:
                throw new IllegalArgumentException(how);
        }
        if (read == 1) {
            q.read1 = bases;
            q.quals1 = quals;
        } else {
            q.read2 = bases;
            q.quals2 = quals;
        }
        // p passes MIN_MEAN_QUALITY of up to 30
        final ElcResult result = run(writePairs(defaultHeader(), Arrays.asList(p, q)),
                "MIN_GROUP_COUNT=1", "MIN_MEAN_QUALITY=" + minMeanQuality);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), kept ? histogram(2, 1) : histogram(1, 1));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", kept ? 2 : 1, kept ? 1 : 0, 0);
    }

    /**
     * With MIN_IDENTICAL_BASES=1 groups are keyed on the first base of each read and the mean group size is
     * max(1, pairs / 4^2). A group of 5 copies on A/A, a group of 4 copies on C/C, and singletons spread over
     * the other 14 first-base combinations at most 2 per combination.
     */
    private static List<Pair> groupRatioInput(final int singletons, final int lowQualityPairs) {
        final Random random = new Random(61);
        final List<Pair> pairs = new ArrayList<>();
        final Pair aa = new Pair("A" + randomBases(random, LENGTH - 1), "A" + randomBases(random, LENGTH - 1));
        final Pair cc = new Pair("C" + randomBases(random, LENGTH - 1), "C" + randomBases(random, LENGTH - 1));
        for (int i = 0; i < 5; i++) {
            pairs.add(aa.copy());
        }
        for (int i = 0; i < 4; i++) {
            pairs.add(cc.copy());
        }
        final List<String> combos = new ArrayList<>();
        for (final char a : "ACGT".toCharArray()) {
            for (final char b : "ACGT".toCharArray()) {
                if (!(a == 'A' && b == 'A') && !(a == 'C' && b == 'C')) {
                    combos.add("" + a + b);
                }
            }
        }
        for (int i = 0; i < singletons; i++) {
            final String combo = combos.get(i % combos.size());
            pairs.add(new Pair(combo.charAt(0) + randomBases(random, LENGTH - 1), combo.charAt(1) + randomBases(random, LENGTH - 1)));
        }
        for (int i = 0; i < lowQualityPairs; i++) {
            final Pair low = new Pair("G" + randomBases(random, LENGTH - 1), "T" + randomBases(random, LENGTH - 1));
            low.quals1 = quals(LENGTH, (byte) 2);
            pairs.add(low);
        }
        return pairs;
    }

    @DataProvider(name = "maxGroupRatio")
    public Object[][] maxGroupRatio() {
        return new Object[][]{
                // 32 pairs: mean group size 2; groups larger than 2 * ratio are skipped
                {23, 0, 2, histogram(4, 1, 1, 23)},
                {23, 0, 3, histogram(5, 1, 4, 1, 1, 23)},
                // filtered pairs still count toward the mean group size
                {21, 2, 2, histogram(4, 1, 1, 21)},
                // 30 pairs: mean group size 1, groups larger than 2 are skipped
                {21, 0, 2, histogram(1, 21)},
        };
    }

    /**
     * A group of exactly MAX_GROUP_RATIO times the mean group size is kept (the first and third rows). The
     * documentation does not specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "maxGroupRatio")
    public void testMaxGroupRatio(final int singletons, final int lowQualityPairs, final int ratio,
                                  final Map<Integer, Double> expected) throws IOException {
        final ElcResult result = run(writePairs(defaultHeader(), groupRatioInput(singletons, lowQualityPairs)),
                "MIN_GROUP_COUNT=1", "MIN_IDENTICAL_BASES=1", "MAX_GROUP_RATIO=" + ratio);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), expected);
    }

    @DataProvider(name = "barcodes")
    public Object[][] barcodes() {
        return new Object[][]{
                // {argument, tag, read-1 values, read-2 values, expected histogram}
                {null, "BC", new String[]{"AAAA", "AAAA", "CCCC"}, new String[]{null, null, null}, histogram(3, 1)},
                {"BARCODE_TAG=BC", "BC", new String[]{"AAAA", "AAAA", "CCCC"}, new String[]{null, null, null}, histogram(2, 1, 1, 1)},
                // BARCODE_TAG is taken from the first read of the pair
                {"BARCODE_TAG=BC", "BC", new String[]{"AAAA", "AAAA", "AAAA"}, new String[]{"AAAA", "CCCC", "GGGG"}, histogram(3, 1)},
                // a missing tag only matches other missing tags
                {"BARCODE_TAG=BC", "BC", new String[]{"AAAA", null, null}, new String[]{null, null, null}, histogram(2, 1, 1, 1)},
                {"READ_ONE_BARCODE_TAG=RX", "RX", new String[]{"AAAA", "CCCC", "AAAA"}, new String[]{null, null, null}, histogram(2, 1, 1, 1)},
                {"READ_ONE_BARCODE_TAG=RX", "RX", new String[]{"AAAA", "AAAA", "AAAA"}, new String[]{"AAAA", "CCCC", "GGGG"}, histogram(3, 1)},
                {"READ_TWO_BARCODE_TAG=RX", "RX", new String[]{"AAAA", "CCCC", "GGGG"}, new String[]{"TTTT", "TTTT", "TTTT"}, histogram(3, 1)},
                {"READ_TWO_BARCODE_TAG=RX", "RX", new String[]{"AAAA", "AAAA", "AAAA"}, new String[]{"TTTT", "GGGG", "TTTT"}, histogram(2, 1, 1, 1)},
        };
    }

    /**
     * BARCODE_TAG is read from the first read of each pair only (the third row). The documentation does not specify
     * this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "barcodes")
    public void testBarcodes(final String arg, final String tag, final String[] read1Values, final String[] read2Values,
                             final Map<Integer, Double> expected) throws IOException {
        final Pair template = randomPairs(71, 1).get(0);
        final List<Pair> pairs = new ArrayList<>();
        for (int i = 0; i < read1Values.length; i++) {
            final Pair p = template.copy();
            if (read1Values[i] != null) {
                p.tags1.put(tag, read1Values[i]);
            }
            if (read2Values[i] != null) {
                p.tags2.put(tag, read2Values[i]);
            }
            pairs.add(p);
        }
        final File input = writePairs(defaultHeader(), pairs);
        for (final String ram : new String[]{null, "MAX_RECORDS_IN_RAM=2"}) {
            final List<String> args = new ArrayList<>();
            args.add("MIN_GROUP_COUNT=1");
            if (arg != null) {
                args.add(arg);
            }
            if (ram != null) {
                args.add(ram);
            }
            final ElcResult result = run(input, args.toArray(new String[0]));
            Assert.assertEquals(result.histogramByLibrary.get("libA"), expected, String.valueOf(ram));
        }
    }

    @Test
    public void testBarcodesAreHonouredInLargeGroups() throws IOException {
        // 120 identical pairs in one group, alternating between two barcodes
        final Pair template = randomPairs(73, 1).get(0);
        final List<Pair> pairs = new ArrayList<>();
        for (int i = 0; i < 120; i++) {
            final Pair p = template.copy();
            p.tags1.put("BC", i % 2 == 0 ? "AAAA" : "CCCC");
            pairs.add(p);
        }
        final ElcResult result = run(writePairs(defaultHeader(), pairs), "MIN_GROUP_COUNT=1", "BARCODE_TAG=BC");
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(60, 2));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 120, 118, 0);
    }

    /**
     * All pairs share their first five bases; copies of a template differ from it by at most one base after the
     * prefix, so copies of one template are mutual duplicates and pairs of different templates are not.
     */
    private static List<Pair> templatesInOneGroup(final int... copiesPerTemplate) {
        final Random random = new Random(81);
        final List<Pair> templates = new ArrayList<>();
        for (int t = 0; t < copiesPerTemplate.length; t++) {
            templates.add(new Pair("ACGTA" + randomBases(random, LENGTH - 5), "TTGCA" + randomBases(random, LENGTH - 5)));
        }
        final List<Pair> pairs = new ArrayList<>();
        final int max = Arrays.stream(copiesPerTemplate).max().orElse(0);
        for (int c = 0; c < max; c++) {
            for (int t = 0; t < templates.size(); t++) {
                if (c < copiesPerTemplate[t]) {
                    final Pair template = templates.get(t);
                    pairs.add(c == 0 ? template.copy() : new Pair(mutate(template.read1, 5 + c), template.read2));
                }
            }
        }
        return pairs;
    }

    private static int[] repeat(final int value, final int times) {
        final int[] out = new int[times];
        Arrays.fill(out, value);
        return out;
    }

    @DataProvider(name = "groupSizes")
    public Object[][] groupSizes() {
        final int[] mixed = new int[33];
        Arrays.fill(mixed, 3);
        mixed[0] = 5;
        return new Object[][]{
                {repeat(3, 33), histogram(3, 33)},
                {repeat(4, 25), histogram(4, 25)},
                {mixed, histogram(5, 1, 3, 32)},
        };
    }

    @Test(dataProvider = "groupSizes")
    public void testGroupsBelowAndAboveOneHundredPairs(final int[] copiesPerTemplate, final Map<Integer, Double> expected) throws IOException {
        final List<Pair> pairs = templatesInOneGroup(copiesPerTemplate);
        final ElcResult result = run(writePairs(defaultHeader(), pairs), "MIN_GROUP_COUNT=1");
        Assert.assertEquals(result.histogramByLibrary.get("libA"), expected);
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", pairs.size(), pairs.size() - copiesPerTemplate.length, 0);
    }

    @DataProvider(name = "greedyOrder")
    public Object[][] greedyOrder() {
        // a~b and b~c within 3 mismatches, a and c 4 apart: the first pair of a set collects its own duplicates
        return new Object[][]{
                {"abc", histogram(2, 1, 1, 1), null},
                {"bac", histogram(3, 1), null},
                {"cab", histogram(2, 1, 1, 1), null},
                {"bac", histogram(3, 1), "MAX_RECORDS_IN_RAM=1"},
        };
    }

    /**
     * The first pair in input order anchors a duplicate set and collects every pair within MAX_DIFF_RATE of it, so the
     * chain a~b~c, with a and c too far apart, splits differently in different orders. The documentation does not
     * specify this; the expected value follows the current implementation.
     */
    @Test(dataProvider = "greedyOrder")
    public void testFirstPairInInputOrderAnchorsEachSet(final String order, final Map<Integer, Double> expected,
                                                        final String extraArg) throws IOException {
        final Pair a = randomPairs(91, 1).get(0);
        final Pair b = new Pair(mutate(a.read1, 10, 20), a.read2);
        final Pair c = new Pair(mutate(b.read1, 30, 40), a.read2);
        final Map<Character, Pair> byName = new HashMap<>();
        byName.put('a', a);
        byName.put('b', b);
        byName.put('c', c);
        final List<Pair> pairs = new ArrayList<>();
        for (final char ch : order.toCharArray()) {
            pairs.add(byName.get(ch).copy());
        }
        final File input = writePairs(defaultHeader(), pairs);
        final ElcResult result = extraArg == null ?
                run(input, "MIN_GROUP_COUNT=1", "MAX_DIFF_RATE=0.035") :
                run(input, "MIN_GROUP_COUNT=1", "MAX_DIFF_RATE=0.035", extraArg);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), expected);
    }

    /** Four identical pairs on tile 1101; the first is the representative. */
    private static List<Pair> opticalSet(final String readGroup) {
        final Pair template = randomPairs(101, 1).get(0);
        final int[][] xy = {{1000, 1000}, {1100, 1000}, {1000, 1101}, {5000, 5000}};
        final List<Pair> pairs = new ArrayList<>();
        for (final int[] p : xy) {
            final Pair copy = template.copy().at(1101, p[0], p[1]);
            copy.readGroup = readGroup;
            pairs.add(copy);
        }
        return pairs;
    }

    @DataProvider(name = "opticalSetOfFour")
    public Object[][] opticalSetOfFour() {
        return new Object[][]{
                {new String[]{"MIN_GROUP_COUNT=1"}, 4, 1},
                {new String[]{"MIN_GROUP_COUNT=1", "OPTICAL_DUPLICATE_PIXEL_DISTANCE=101"}, 4, 2},
                {new String[]{"MIN_GROUP_COUNT=1", "OPTICAL_DUPLICATE_PIXEL_DISTANCE=99"}, 4, 0},
                {new String[]{"MIN_GROUP_COUNT=1", "OPTICAL_DUPLICATE_PIXEL_DISTANCE=5000"}, 4, 3},
                {new String[]{"MIN_GROUP_COUNT=1", "READ_NAME_REGEX=null"}, 4, 0},
                {new String[]{"MIN_GROUP_COUNT=1", "MAX_OPTICAL_DUPLICATE_SET_SIZE=3"}, 4, 0},
                {new String[]{"MIN_GROUP_COUNT=1", "MAX_OPTICAL_DUPLICATE_SET_SIZE=4"}, 4, 1},
                // the only bin holds one set, below MIN_GROUP_COUNT=2, so its optical duplicates are not counted
                {new String[]{}, 0, 0},
        };
    }

    @Test(dataProvider = "opticalSetOfFour")
    public void testOpticalDuplicates(final String[] args, final long examined, final long optical) throws IOException {
        // records without a read group, so the library is Unknown however the read names parse
        final File input = writePairs(defaultHeader(), opticalSet(null));
        final ElcResult result = run(input, args);
        Assert.assertEquals(result.histogramByLibrary.get("Unknown"), histogram(4, 1));
        assertMetrics(result.metricsByLibrary.get("Unknown"), "Unknown", examined, examined == 0 ? 0 : 3, optical);
    }

    @DataProvider(name = "opticalPairs")
    public Object[][] opticalPairs() {
        return new Object[][]{
                // {dx, dy, same tile, same read group, optical}
                {100, 0, true, true, 1},
                {101, 0, true, true, 0},
                {0, 100, true, true, 1},
                {0, 101, true, true, 0},
                {-100, -100, true, true, 1},
                {0, 0, false, true, 0},
                {0, 0, true, false, 0},
        };
    }

    @Test(dataProvider = "opticalPairs")
    public void testOpticalDuplicatePairs(final int dx, final int dy, final boolean sameTile, final boolean sameReadGroup,
                                          final long optical) throws IOException {
        final Pair p = randomPairs(103, 1).get(0).at(1101, 2000, 2000);
        final Pair q = p.copy().at(sameTile ? 1101 : 1102, 2000 + dx, 2000 + dy);
        q.readGroup = sameReadGroup ? RG_A : "rgB";
        final File input = writePairs(header(RG_A, "libA", "rgB", "libA"), Arrays.asList(p, q));
        final ElcResult result = run(input, "MIN_GROUP_COUNT=1");
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 2, 1, optical);
    }

    @Test
    public void testLibrariesAreCountedSeparately() throws IOException {
        final Pair template = randomPairs(111, 1).get(0);
        final String[] readGroups = {RG_A, "rgB", "rgC", RG_A, "rgB", "rgD"};
        final List<Pair> pairs = new ArrayList<>();
        for (final String rg : readGroups) {
            final Pair p = template.copy();
            p.readGroup = rg;
            pairs.add(p);
        }
        final File input = writePairs(header(RG_A, "libA", "rgB", "libB", "rgC", "libA", "rgD", null), pairs);
        final ElcResult result = run(input, "MIN_GROUP_COUNT=1");

        Assert.assertEquals(result.metricsByLibrary.keySet(), new java.util.HashSet<>(Arrays.asList("libA", "libB", "Unknown")));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 3, 2, 0);
        assertMetrics(result.metricsByLibrary.get("libB"), "libB", 2, 1, 0);
        assertMetrics(result.metricsByLibrary.get("Unknown"), "Unknown", 1, 0, 0);
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(3, 1));
        Assert.assertEquals(result.histogramByLibrary.get("libB"), histogram(2, 1));
        Assert.assertEquals(result.histogramByLibrary.get("Unknown"), histogram(1, 1));
    }

    @Test
    public void testMultipleInputs() throws IOException {
        final List<Pair> t = randomPairs(121, 2);
        final SAMFileHeader header = defaultHeader();
        final List<SAMRecord> first = toRecords(header, Arrays.asList(t.get(0).copy().at(11, 1000, 1000), t.get(1).copy().at(12, 1000, 1000)), "RUN:A");
        final List<SAMRecord> second = toRecords(header, Arrays.asList(t.get(0).copy().at(21, 1000, 1000), t.get(1).copy().at(22, 1000, 1000)), "RUN:B");
        // mates of one pair split across the two files are never joined
        final List<SAMRecord> split = toRecords(header, Collections.singletonList(t.get(0).copy().at(31, 1000, 1000)), "RUN:C");
        first.add(split.get(0));
        second.add(split.get(1));
        final ElcResult result = run(Arrays.asList(writeSam(header, first), writeSam(header, second)), "MIN_GROUP_COUNT=1");
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(2, 2));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 4, 2, 0);
    }

    private static SAMRecord alignedMate(final SAMFileHeader header, final String name, final boolean first,
                                         final String sequencedBases, final boolean negative, final int start) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName(name);
        rec.setReadPairedFlag(true);
        rec.setFirstOfPairFlag(first);
        rec.setSecondOfPairFlag(!first);
        rec.setReferenceIndex(0);
        rec.setAlignmentStart(start);
        rec.setCigarString(sequencedBases.length() + "M");
        rec.setMappingQuality(60);
        rec.setReadNegativeStrandFlag(negative);
        rec.setMateReferenceIndex(0);
        rec.setMateAlignmentStart(start + 100);
        rec.setMateNegativeStrandFlag(!negative);
        rec.setReadBases((negative ? SequenceUtil.reverseComplement(sequencedBases) : sequencedBases).getBytes());
        rec.setBaseQualities(quals(sequencedBases.length(), QUAL));
        rec.setAttribute("RG", RG_A);
        return rec;
    }

    @Test
    public void testAlignedRecordsAreRestoredToSequencingOrientation() throws IOException {
        final SAMFileHeader header = defaultHeader();
        final Random random = new Random(131);
        final String s1 = randomBases(random, LENGTH);
        final String s2 = randomBases(random, LENGTH);
        final String other = randomBases(random, LENGTH);
        final List<SAMRecord> records = new ArrayList<>();

        records.add(alignedMate(header, "RUN:0:1:1000:1000", true, s1, false, 100));
        // records that must not take part: secondary, supplementary, unpaired, and paired without a read number
        final SAMRecord secondary = alignedMate(header, "RUN:0:1:1000:1000", true, other, false, 5000);
        secondary.setSecondaryAlignment(true);
        records.add(secondary);
        final SAMRecord supplementary = alignedMate(header, "RUN:0:1:1000:1000", false, other, false, 6000);
        supplementary.setSupplementaryAlignmentFlag(true);
        records.add(supplementary);
        final SAMRecord unpaired = alignedMate(header, "RUN:0:1:1000:1000", true, other, false, 7000);
        unpaired.setReadPairedFlag(false);
        unpaired.setFirstOfPairFlag(false);
        unpaired.setSecondOfPairFlag(false);
        unpaired.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
        unpaired.setMateAlignmentStart(0);
        unpaired.setMateNegativeStrandFlag(false);
        records.add(unpaired);
        final SAMRecord noReadNumber = alignedMate(header, "RUN:0:1:1000:1000", true, other, false, 8000);
        noReadNumber.setFirstOfPairFlag(false);
        records.add(noReadNumber);
        records.add(alignedMate(header, "RUN:0:1:1000:1000", false, s2, true, 200));

        // the same fragment sequenced from the other strand
        records.add(alignedMate(header, "RUN:1:2:1000:1000", true, s1, true, 300));
        records.add(alignedMate(header, "RUN:1:2:1000:1000", false, s2, false, 400));

        // stored bases equal to s1 on the reverse strand are a different sequenced read
        records.add(alignedMate(header, "RUN:2:3:1000:1000", true, SequenceUtil.reverseComplement(s1), true, 500));
        records.add(alignedMate(header, "RUN:2:3:1000:1000", false, s2, false, 600));

        // a paired record with neither read number fails strict validation
        final ElcResult result = run(writeSam(header, records), "MIN_GROUP_COUNT=1", "VALIDATION_STRINGENCY=LENIENT");
        Assert.assertEquals(result.histogramByLibrary.get("libA"), histogram(2, 1, 1, 1));
        assertMetrics(result.metricsByLibrary.get("libA"), "libA", 3, 1, 0);
    }

    @Test
    public void testInvalidArgumentsFailTheCommandLine() throws IOException {
        final File input = writePairs(defaultHeader(), randomPairs(141, 2));
        final File output = File.createTempFile("elcBehavior.", ".elc_metrics");
        output.deleteOnExit();
        for (final String[] bad : new String[][]{{"MIN_IDENTICAL_BASES=0"}, {"MIN_IDENTICAL_BASES=6", "MAX_READ_LENGTH=5"}}) {
            final List<String> args = new ArrayList<>(Arrays.asList("INPUT=" + input.getAbsolutePath(), "OUTPUT=" + output.getAbsolutePath()));
            args.addAll(Arrays.asList(bad));
            Assert.assertEquals(runPicardCommandLine(args), 1, Arrays.toString(bad));
        }
    }
}
