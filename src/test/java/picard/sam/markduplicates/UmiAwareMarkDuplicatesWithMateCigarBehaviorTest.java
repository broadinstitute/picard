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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * End-to-end tests of {@link UmiAwareMarkDuplicatesWithMateCigar}: how UMIs split duplicate sets at, below and past
 * MAX_EDIT_DISTANCE_TO_JOIN, transitive joining, N bases, hyphenated and duplex UMIs, missing UMIs, the UMI tag name,
 * the molecular identifier written to each record, and the values of the UMI metrics.
 */
public class UmiAwareMarkDuplicatesWithMateCigarBehaviorTest {
    private static final String LIBRARY = "libA";
    private static final double TOLERANCE = 1e-5;

    private static class Result {
        final Map<String, SAMRecord> records;
        final Map<String, DuplicationMetrics> metrics;
        final List<UmiMetrics> umiMetrics;

        Result(final File output, final File metrics, final File umiMetrics) throws IOException {
            this.records = byKey(readRecords(output));
            this.metrics = readDuplicationMetrics(metrics);
            this.umiMetrics = readMetrics(umiMetrics);
        }
    }

    private static Result run(final File input, final String... extraArgs) throws IOException {
        final File output = newTempFile("UmiAwareMarkDuplicatesWithMateCigarBehaviorTest", ".sam");
        final File metrics = newTempFile("UmiAwareMarkDuplicatesWithMateCigarBehaviorTest", ".metrics");
        final File umiMetrics = newTempFile("UmiAwareMarkDuplicatesWithMateCigarBehaviorTest", ".umi_metrics");
        final List<String> args = new ArrayList<>(List.of(
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath(),
                "METRICS_FILE=" + metrics.getAbsolutePath(),
                "UMI_METRICS_FILE=" + umiMetrics.getAbsolutePath(),
                "DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH"));
        args.addAll(Arrays.asList(extraArgs));
        Assert.assertEquals(new UmiAwareMarkDuplicatesWithMateCigar().instanceMain(args.toArray(new String[0])), 0);
        return new Result(output, metrics, umiMetrics);
    }

    private static SAMRecordSetBuilder newUmiBuilder() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", LIBRARY);
        return builder;
    }

    /** Adds FR pairs u0, u1, ... with the first end at 1000 and the second at 1200, tagged with the given UMIs. */
    private static File pairsWithUmis(final String tag, final List<String> umis) throws IOException {
        final SAMRecordSetBuilder builder = newUmiBuilder();
        for (int i = 0; i < umis.size(); i++) {
            final String umi = umis.get(i);
            addPair(builder, "u" + i, 0, 1000, 1200).forEach(r -> r.setAttribute(tag, umi));
        }
        return writeCoordinateSorted(builder);
    }

    private static void assertPairsDuplicateAndMolecularIdentifiers(final Result result, final List<Boolean> duplicates,
                                                                    final List<String> molecularIdentifiers) {
        for (int i = 0; i < duplicates.size(); i++) {
            for (final String end : List.of("/1", "/2")) {
                final SAMRecord record = result.records.get("u" + i + end);
                Assert.assertEquals(record.getDuplicateReadFlag(), (boolean) duplicates.get(i), record.getSAMString());
                Assert.assertEquals(record.getAttribute("MI"), molecularIdentifiers.get(i), record.getSAMString());
            }
        }
    }

    @DataProvider(name = "umiSets")
    public Object[][] umiSets() {
        return new Object[][]{
                // UMIs of the pairs, MAX_EDIT_DISTANCE_TO_JOIN, duplicate flag of each pair, inferred UMI of each pair
                {List.of("AAAA", "AAAA", "AAAT"), 0, List.of(false, true, false), List.of("AAAA", "AAAA", "AAAT")},
                {List.of("AAAA", "AAAA", "AAAT"), 1, List.of(false, true, true), List.of("AAAA", "AAAA", "AAAA")},
                {List.of("AAAA", "AAAA", "AATT"), 1, List.of(false, true, false), List.of("AAAA", "AAAA", "AATT")},
                {List.of("AAAA", "AAAA", "AATT"), 2, List.of(false, true, true), List.of("AAAA", "AAAA", "AAAA")},
                {List.of("AAAA", "AAAA", "ATTT"), 2, List.of(false, true, false), List.of("AAAA", "AAAA", "ATTT")},
                // Joined through AAAT although AAAA and AATT are two apart.
                {List.of("AAAA", "AAAA", "AAAT", "AATT"), 1, List.of(false, true, true, true), List.of("AAAA", "AAAA", "AAAA", "AAAA")},
                // The minority UMI without N is inferred over one with N; the keeper is still chosen by name.
                {List.of("AANA", "AAAA", "AAAA"), 1, List.of(false, true, true), List.of("AAAA", "AAAA", "AAAA")},
                {List.of("AANA", "AANA", "AANA", "AAAA"), 1, List.of(false, true, true, true), List.of("AAAA", "AAAA", "AAAA", "AAAA")},
                // N matches N at the same position, but not a base or an N elsewhere.
                {List.of("ANAA", "ANAA"), 0, List.of(false, true), List.of("ANAA", "ANAA")},
                {List.of("ANAA", "AANA"), 1, List.of(false, false), List.of("ANAA", "AANA")},
                // When every UMI has an N, the one with the fewest is inferred.
                {List.of("NNAA", "NAAA", "NNAA"), 1, List.of(false, true, true), List.of("NAAA", "NAAA", "NAAA")},
                {List.of("AC-GT", "AC-GT", "AC-GA"), 1, List.of(false, true, true), List.of("AC-GT", "AC-GT", "AC-GT")},
                {List.of("AC-GT", "AC-GT", "AC-GA"), 0, List.of(false, true, false), List.of("AC-GT", "AC-GT", "AC-GA")},
        };
    }

    @Test(dataProvider = "umiSets")
    public void testUmisSplitDuplicateSets(final List<String> umis, final int maxEditDistance, final List<Boolean> duplicates,
                                           final List<String> inferredUmis) throws IOException {
        final Result result = run(pairsWithUmis("RX", umis), "MAX_EDIT_DISTANCE_TO_JOIN=" + maxEditDistance,
                "MOLECULAR_IDENTIFIER_TAG=MI");

        final List<String> molecularIdentifiers = new ArrayList<>();
        inferredUmis.forEach(umi -> molecularIdentifiers.add("chr1:1200/" + umi));
        assertPairsDuplicateAndMolecularIdentifiers(result, duplicates, molecularIdentifiers);
        for (int i = 0; i < umis.size(); i++) {
            Assert.assertEquals(result.records.get("u" + i + "/1").getAttribute("RX"), umis.get(i));
        }

        final long duplicatePairs = duplicates.stream().filter(d -> d).count();
        assertDuplicationMetrics(result.metrics.get(LIBRARY), LIBRARY, 0, umis.size(), 0, 0, 0, duplicatePairs, 0);

        // One duplicate set per end of the pairs, each split into one set per molecule.
        final UmiMetrics umiMetrics = result.umiMetrics.get(0);
        Assert.assertEquals(result.umiMetrics.size(), 1);
        Assert.assertEquals(umiMetrics.DUPLICATE_SETS_IGNORING_UMI, 2);
        Assert.assertEquals(umiMetrics.DUPLICATE_SETS_WITH_UMI, 2 * new HashSet<>(inferredUmis).size());
    }

    private static double entropyBase4(final double... counts) {
        final double total = Arrays.stream(counts).sum();
        double entropy = 0;
        for (final double count : counts) {
            entropy -= count / total * Math.log(count / total) / Math.log(4);
        }
        return entropy;
    }

    @DataProvider(name = "umiMetrics")
    public Object[][] umiMetrics() {
        return new Object[][]{
                // Six reads: AAAA x4 and AAAT x2; the two AAAT reads are corrected, so 2 errors in 24 bases.
                {List.of("AAAA", "AAAA", "AAAT"), 1,
                        new UmiMetrics(LIBRARY, 4, 2, 1, 2, 2, 2, 0, entropyBase4(4, 2), 11, 0)},
                // Reads whose UMI has an N only count toward PCT_UMI_WITH_N: 2 of 8.
                {List.of("AANA", "AAAA", "AAAA", "AAAT"), 1,
                        new UmiMetrics(LIBRARY, 4, 2, 1, 2, 2, 2, 0, entropyBase4(4, 2), 11, 0.25)},
                // Fourteen reads, five UMIs joined into two: AAAAAA x10 (six corrected, 1 error each), CCCCCC x4.
                {List.of("AAAAAA", "AAAAAT", "AAAAAA", "AAAAAG", "CCCCCC", "CCCCCC", "AAAATA"), 1,
                        new UmiMetrics(LIBRARY, 6, 5, 2, 6, 2, 4,
                                entropyBase4(10, 4), entropyBase4(4, 2, 2, 4, 2), 11, 0)},
        };
    }

    @Test(dataProvider = "umiMetrics")
    public void testUmiMetricsValues(final List<String> umis, final int maxEditDistance, final UmiMetrics expected) throws IOException {
        final Result result = run(pairsWithUmis("RX", umis), "MAX_EDIT_DISTANCE_TO_JOIN=" + maxEditDistance);
        Assert.assertEquals(result.umiMetrics.size(), 1);
        assertUmiMetrics(result.umiMetrics.get(0), expected);
    }

    private static void assertUmiMetrics(final UmiMetrics actual, final UmiMetrics expected) {
        Assert.assertEquals(actual.LIBRARY, expected.LIBRARY);
        Assert.assertEquals(actual.MEAN_UMI_LENGTH, expected.MEAN_UMI_LENGTH, TOLERANCE, "MEAN_UMI_LENGTH");
        Assert.assertEquals(actual.OBSERVED_UNIQUE_UMIS, expected.OBSERVED_UNIQUE_UMIS, "OBSERVED_UNIQUE_UMIS");
        Assert.assertEquals(actual.INFERRED_UNIQUE_UMIS, expected.INFERRED_UNIQUE_UMIS, "INFERRED_UNIQUE_UMIS");
        Assert.assertEquals(actual.OBSERVED_BASE_ERRORS, expected.OBSERVED_BASE_ERRORS, "OBSERVED_BASE_ERRORS");
        Assert.assertEquals(actual.DUPLICATE_SETS_IGNORING_UMI, expected.DUPLICATE_SETS_IGNORING_UMI, "DUPLICATE_SETS_IGNORING_UMI");
        Assert.assertEquals(actual.DUPLICATE_SETS_WITH_UMI, expected.DUPLICATE_SETS_WITH_UMI, "DUPLICATE_SETS_WITH_UMI");
        Assert.assertEquals(actual.OBSERVED_UMI_ENTROPY, expected.OBSERVED_UMI_ENTROPY, TOLERANCE, "OBSERVED_UMI_ENTROPY");
        Assert.assertEquals(actual.INFERRED_UMI_ENTROPY, expected.INFERRED_UMI_ENTROPY, TOLERANCE, "INFERRED_UMI_ENTROPY");
        Assert.assertEquals(actual.UMI_BASE_QUALITIES, expected.UMI_BASE_QUALITIES, TOLERANCE, "UMI_BASE_QUALITIES");
        Assert.assertEquals(actual.PCT_UMI_WITH_N, expected.PCT_UMI_WITH_N, TOLERANCE, "PCT_UMI_WITH_N");
    }

    /**
     * d1 and d3 come from the top strand and d2, with its halves swapped, from the bottom strand of the same molecule;
     * d3 has one error. d4 is another molecule.
     */
    private static File duplexInput() throws IOException {
        final SAMRecordSetBuilder builder = newUmiBuilder();
        addPair(builder, "d1", 0, 1000, false, "50M", 1200, true, "50M").forEach(r -> r.setAttribute("RX", "AAA-CCC"));
        addPair(builder, "d2", 0, 1200, true, "50M", 1000, false, "50M").forEach(r -> r.setAttribute("RX", "CCC-AAA"));
        addPair(builder, "d3", 0, 1000, false, "50M", 1200, true, "50M").forEach(r -> r.setAttribute("RX", "AAA-CCG"));
        addPair(builder, "d4", 0, 1000, false, "50M", 1200, true, "50M").forEach(r -> r.setAttribute("RX", "GGG-TTT"));
        return writeCoordinateSorted(builder);
    }

    private static void assertPairs(final Result result, final Map<String, Boolean> duplicates, final Map<String, String> molecularIdentifiers) {
        for (final String name : duplicates.keySet()) {
            for (final String end : List.of("/1", "/2")) {
                final SAMRecord record = result.records.get(name + end);
                Assert.assertEquals(record.getDuplicateReadFlag(), (boolean) duplicates.get(name), record.getSAMString());
                Assert.assertEquals(record.getAttribute("MI"), molecularIdentifiers.get(name), record.getSAMString());
            }
        }
    }

    @Test
    public void testDuplexUmis() throws IOException {
        final Result result = run(duplexInput(), "DUPLEX_UMI=true", "MAX_EDIT_DISTANCE_TO_JOIN=1", "MOLECULAR_IDENTIFIER_TAG=MI");
        assertPairs(result,
                Map.of("d1", false, "d2", true, "d3", true, "d4", false),
                Map.of("d1", "chr1:1200/AAA-CCC/A", "d2", "chr1:1200/AAA-CCC/B",
                        "d3", "chr1:1200/AAA-CCC/A", "d4", "chr1:1200/GGG-TTT/A"));

        // Observed top-strand UMIs AAA-CCC x4, AAA-CCG x2, GGG-TTT x2; inferred AAA-CCC x6, GGG-TTT x2.
        final UmiMetrics metrics = result.umiMetrics.get(0);
        Assert.assertEquals(metrics.MEAN_UMI_LENGTH, 6.0);
        Assert.assertEquals(metrics.OBSERVED_UNIQUE_UMIS, 3);
        Assert.assertEquals(metrics.INFERRED_UNIQUE_UMIS, 2);
        Assert.assertEquals(metrics.OBSERVED_BASE_ERRORS, 2);
        Assert.assertEquals(metrics.DUPLICATE_SETS_IGNORING_UMI, 2);
        Assert.assertEquals(metrics.DUPLICATE_SETS_WITH_UMI, 4);
        Assert.assertEquals(metrics.OBSERVED_UMI_ENTROPY, entropyBase4(4, 2, 2), TOLERANCE);
        Assert.assertEquals(metrics.INFERRED_UMI_ENTROPY, entropyBase4(6, 2), TOLERANCE);
        Assert.assertEquals(metrics.PCT_UMI_WITH_N, 0.0, TOLERANCE);
    }

    /**
     * Without DUPLEX_UMI, d1's AAA-CCC and d3's AAA-CCG, equally common, are joined and both get AAA-CCC, the UMI of
     * the first record. The documentation does not specify this; the expected value follows the current
     * implementation.
     */
    @Test
    public void testDuplexUmisReadAsSingleStranded() throws IOException {
        final Result result = run(duplexInput(), "MAX_EDIT_DISTANCE_TO_JOIN=1", "MOLECULAR_IDENTIFIER_TAG=MI");
        assertPairs(result,
                Map.of("d1", false, "d2", false, "d3", true, "d4", false),
                Map.of("d1", "chr1:1200/AAA-CCC", "d2", "chr1:1200/CCC-AAA",
                        "d3", "chr1:1200/AAA-CCC", "d4", "chr1:1200/GGG-TTT"));
        Assert.assertEquals(result.umiMetrics.get(0).DUPLICATE_SETS_WITH_UMI, 6);
    }

    @Test
    public void testUmiTagName() throws IOException {
        final Result result = run(pairsWithUmis("ZU", List.of("AAAA", "AAAA", "CCCC")), "UMI_TAG_NAME=ZU",
                "MOLECULAR_IDENTIFIER_TAG=MI");
        assertPairsDuplicateAndMolecularIdentifiers(result, List.of(false, true, false),
                List.of("chr1:1200/AAAA", "chr1:1200/AAAA", "chr1:1200/CCCC"));
        Assert.assertEquals(result.records.get("u2/2").getAttribute("ZU"), "CCCC");
        Assert.assertNull(result.records.get("u2/2").getAttribute("RX"));
    }

    @Test
    public void testRemoveDuplicates() throws IOException {
        final Result result = run(pairsWithUmis("RX", List.of("AAAA", "AAAA", "CCCC")), "REMOVE_DUPLICATES=true");
        Assert.assertEquals(result.records.keySet(), Set.of("u0/1", "u0/2", "u2/1", "u2/2"));
        result.records.values().forEach(r -> Assert.assertFalse(r.getDuplicateReadFlag()));
        assertDuplicationMetrics(result.metrics.get(LIBRARY), LIBRARY, 0, 3, 0, 0, 0, 1, 0);
    }

    @Test
    public void testNoMolecularIdentifierTagUnlessRequested() throws IOException {
        final Result result = run(pairsWithUmis("RX", List.of("AAAA", "AAAA", "CCCC")));
        final List<String> noTags = Arrays.asList(null, null, null);
        assertPairsDuplicateAndMolecularIdentifiers(result, List.of(false, true, false), noTags);
    }

    /**
     * With ALLOW_MISSING_UMIS, pairs without a UMI get no molecular identifier tag even when MOLECULAR_IDENTIFIER_TAG
     * is set. The documentation does not specify this; the expected value follows the current implementation.
     */
    @Test
    public void testMissingUmisAllowed() throws IOException {
        final SAMRecordSetBuilder builder = newUmiBuilder();
        addPair(builder, "u0", 0, 1000, 1200);
        addPair(builder, "u1", 0, 1000, 1200);
        addPair(builder, "u2", 0, 1000, 1200);
        final Result result = run(writeCoordinateSorted(builder), "ALLOW_MISSING_UMIS=true", "MOLECULAR_IDENTIFIER_TAG=MI");

        final List<String> noTags = Arrays.asList(null, null, null);
        assertPairsDuplicateAndMolecularIdentifiers(result, List.of(false, true, true), noTags);
        result.records.values().forEach(r -> Assert.assertNull(r.getAttribute("RX")));

        final UmiMetrics metrics = result.umiMetrics.get(0);
        Assert.assertEquals(metrics.DUPLICATE_SETS_IGNORING_UMI, 2);
        Assert.assertEquals(metrics.DUPLICATE_SETS_WITH_UMI, 2);
        Assert.assertEquals(metrics.OBSERVED_UNIQUE_UMIS, 0);
        Assert.assertEquals(metrics.INFERRED_UNIQUE_UMIS, 0);
        Assert.assertEquals(metrics.OBSERVED_BASE_ERRORS, 0);
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*u1.*does not contain a UMI.*")
    public void testMissingUmisRejected() throws IOException {
        final SAMRecordSetBuilder builder = newUmiBuilder();
        addPair(builder, "u0", 0, 1000, 1200).forEach(r -> r.setAttribute("RX", "AAAA"));
        addPair(builder, "u1", 0, 1000, 1200);
        run(writeCoordinateSorted(builder));
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*illegal characters.*")
    public void testIllegalUmiCharactersRejected() throws IOException {
        run(pairsWithUmis("RX", List.of("AAAA", "AAUA")));
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*differing lengths.*")
    public void testUmisOfDifferentLengthsRejected() throws IOException {
        final SAMRecordSetBuilder builder = newUmiBuilder();
        addPair(builder, "u0", 0, 1000, 1200).forEach(r -> r.setAttribute("RX", "AAAA"));
        addPair(builder, "u1", 0, 3000, 3200).forEach(r -> r.setAttribute("RX", "AAA"));
        run(writeCoordinateSorted(builder));
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*X-Y.*")
    public void testMalformedDuplexUmiRejected() throws IOException {
        run(pairsWithUmis("RX", List.of("AAACCC", "AAACCC")), "DUPLEX_UMI=true");
    }
}
