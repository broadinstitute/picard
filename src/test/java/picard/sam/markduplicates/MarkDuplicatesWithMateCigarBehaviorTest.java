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
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTag;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * End-to-end tests of {@link MarkDuplicatesWithMateCigar}: which records are marked or removed and in what order they are
 * written, the handling of pairs without a mate CIGAR, the MINIMUM_DISTANCE window, the 5' unclipped coordinates of both
 * ends, optical duplicates, per-library metrics, records spilled to disk by small MAX_RECORDS_IN_RAM and BLOCK_SIZE values,
 * program record chaining, multiple inputs, and the inputs the program rejects.
 */
public class MarkDuplicatesWithMateCigarBehaviorTest {

    private static class Result {
        final List<SAMRecord> records;
        final SAMFileHeader header;
        final Map<String, DuplicationMetrics> metrics;
        final Map<String, Map<Double, Double>> histograms;

        Result(final File output, final File metrics) throws IOException {
            this.records = readRecords(output);
            this.header = readHeader(output);
            this.metrics = readDuplicationMetrics(metrics);
            this.histograms = readHistograms(metrics);
        }
    }

    private static Result run(final List<File> inputs, final String... extraArgs) throws IOException {
        final File output = newTempFile("MarkDuplicatesWithMateCigarBehaviorTest", ".sam");
        final File metrics = newTempFile("MarkDuplicatesWithMateCigarBehaviorTest", ".metrics");
        final List<String> args = new ArrayList<>();
        inputs.forEach(input -> args.add("INPUT=" + input.getAbsolutePath()));
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());
        args.addAll(Arrays.asList(extraArgs));
        Assert.assertEquals(new MarkDuplicatesWithMateCigar().instanceMain(args.toArray(new String[0])), 0);
        return new Result(output, metrics);
    }

    private static Result run(final File input, final String... extraArgs) throws IOException {
        return run(List.of(input), extraArgs);
    }

    @DataProvider(name = "bufferSettings")
    public Object[][] bufferSettings() {
        return new Object[][]{
                {false, 500000, 100000},
                {true, 500000, 100000},
                {false, 2, 1},
                {true, 3, 2},
                {false, 1, 1},
        };
    }

    @Test(dataProvider = "bufferSettings")
    public void testMixedInputMarksDuplicatesInInputOrder(final boolean removeDuplicates, final int maxRecordsInRam,
                                                          final int blockSize) throws IOException {
        final File input = writeCoordinateSorted(mixedInput());
        final List<String> inputOrder = keys(readRecords(input));

        final Result result = run(input, "REMOVE_DUPLICATES=" + removeDuplicates,
                "MAX_RECORDS_IN_RAM=" + maxRecordsInRam, "BLOCK_SIZE=" + blockSize);

        final List<String> expectedOrder = inputOrder.stream()
                .filter(key -> !removeDuplicates || !MIXED_DUPLICATES.contains(key))
                .collect(Collectors.toList());
        Assert.assertEquals(keys(result.records), expectedOrder);
        assertDuplicates(result.records, removeDuplicates ? List.of() : MIXED_DUPLICATES);

        Assert.assertEquals(result.metrics.keySet(), Set.of(MIXED_LIBRARY));
        assertMixedMetrics(result.metrics.get(MIXED_LIBRARY));
    }

    @Test
    public void testInputDuplicateFlagsAreReplaced() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        addPair(builder, "p2", 0, 1000, 1200);
        addFragment(builder, "f1", 0, 3000, false, "50M");
        // The keeper arrives flagged and the duplicate arrives unflagged.
        builder.getRecords().stream().filter(r -> r.getReadName().equals("p1")).forEach(r -> r.setDuplicateReadFlag(true));
        builder.getRecords().stream().filter(r -> r.getReadName().equals("f1")).forEach(r -> r.setDuplicateReadFlag(true));

        final Result result = run(writeCoordinateSorted(builder));
        assertDuplicates(result.records, List.of("p2/1", "p2/2"));
    }

    /**
     * A pair without MC tags is not considered for duplicate marking but still counts in READ_PAIRS_EXAMINED. The
     * documentation does not specify this; the expected value follows the current implementation.
     */
    @Test
    public void testPairsWithoutMateCigarAreSkipped() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        // pX0 would be the keeper by name, but without MC tags it is passed through unexamined.
        final List<SAMRecord> noMateCigar = addPair(builder, "pX0", 0, 1000, 1200);
        addPair(builder, "pX1", 0, 1000, 1200);
        addPair(builder, "pX2", 0, 1000, 1200);
        for (final SAMRecord record : noMateCigar) {
            record.setAttribute(SAMTag.MC.name(), null);
            record.setDuplicateReadFlag(true);
        }

        final Result result = run(writeCoordinateSorted(builder));
        assertDuplicates(result.records, List.of("pX2/1", "pX2/2"));
        Assert.assertNull(byKey(result.records).get("pX0/1").getAttribute(SAMTag.MC.name()));
        Assert.assertEquals(result.records.size(), 6);
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 3, 0, 0, 0, 1, 0);
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*pX0.*no mate cigar.*")
    public void testPairsWithoutMateCigarFailWhenNotSkipped() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "pX1", 0, 1000, 1200);
        addPair(builder, "pX0", 0, 1000, 1200).forEach(r -> r.setAttribute(SAMTag.MC.name(), null));

        run(writeCoordinateSorted(builder), "SKIP_PAIRS_WITH_NO_MATE_CIGAR=false");
    }

    /**
     * Fragments A at 1000 and B at 1000 + d, then C aligned at 1001 + d whose soft clip puts its 5' end at 1000.
     * C can only be compared with A if A is still buffered when C is read, which needs a window of at least d.
     */
    private static File clippedFragmentInput(final int d) throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addFragment(builder, "fragA", 0, 1000, false, "50M");
        addFragment(builder, "fragB", 0, 1000 + d, false, "50M");
        addFragment(builder, "fragC", 0, 1001 + d, false, (d + 1) + "S50M");
        return writeCoordinateSorted(builder);
    }

    @DataProvider(name = "sufficientDistances")
    public Object[][] sufficientDistances() {
        return new Object[][]{
                {150, new String[]{"MINIMUM_DISTANCE=150"}},
                {150, new String[]{"MINIMUM_DISTANCE=400"}},
                // With MINIMUM_DISTANCE=-1 and 50-base reads the window is 100.
                {100, new String[]{}},
                {100, new String[]{"MINIMUM_DISTANCE=-1"}},
        };
    }

    @Test(dataProvider = "sufficientDistances")
    public void testClippedDuplicateWithinMinimumDistance(final int d, final String[] args) throws IOException {
        final Result result = run(clippedFragmentInput(d), args);
        assertDuplicates(result.records, List.of("fragC/0"));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 3, 0, 0, 0, 1, 0, 0);
    }

    @DataProvider(name = "insufficientDistances")
    public Object[][] insufficientDistances() {
        return new Object[][]{
                {150, new String[]{"MINIMUM_DISTANCE=149"}},
                {150, new String[]{"MINIMUM_DISTANCE=10"}},
                {101, new String[]{}},
        };
    }

    @Test(dataProvider = "insufficientDistances", expectedExceptions = PicardException.class,
            expectedExceptionsMessageRegExp = "(?s).*Please increase the minimum distance.*")
    public void testClippedDuplicateBeyondMinimumDistanceFails(final int d, final String[] args) throws IOException {
        run(clippedFragmentInput(d), args);
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*skipped reference bases.*")
    public void testSplicedReadBeyondMinimumDistanceFails() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addFragment(builder, "fragA", 0, 1000, false, "50M");
        addFragment(builder, "fragB", 0, 1150, false, "10M100N40M");
        addFragment(builder, "fragC", 0, 1151, false, "151S50M");
        run(writeCoordinateSorted(builder), "MINIMUM_DISTANCE=149");
    }

    @Test
    public void testFivePrimeEndsIncludeClipsOnBothEnds() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        // 5' ends of cl1: 1000 forward and 1349 reverse.
        addPair(builder, "cl1", 0, 1000, false, "50M", 1300, true, "50M");
        // Soft clips: 1005 - 5 = 1000 and 1300 + 45 - 1 + 5 = 1349.
        addPair(builder, "cl2", 0, 1005, false, "5S45M", 1300, true, "45M5S");
        // Hard clips: 1003 - 3 = 1000 and 1302 + 45 - 1 + 3 = 1349.
        addPair(builder, "cl3", 0, 1003, false, "3H47M", 1302, true, "2H45M3H");
        // The mate's 5' end is 1344.
        addPair(builder, "cl4", 0, 1000, false, "50M", 1300, true, "45M");
        // The first end's 5' end is 1005.
        addPair(builder, "cl5", 0, 1010, false, "5S45M", 1300, true, "50M");

        final Result result = run(writeCoordinateSorted(builder), "MINIMUM_DISTANCE=200");
        assertDuplicates(result.records, List.of("cl2/1", "cl2/2", "cl3/1", "cl3/2"));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 5, 0, 0, 0, 2, 0);
    }

    @DataProvider(name = "opticalSettings")
    public Object[][] opticalSettings() {
        // The histograms count sets by number of pairs: all_sets by size, optical_sets by optical duplicates plus one,
        // and non_optical_sets by size less optical duplicates.
        final Map<Double, Double> noOptical = Map.of();
        return new Object[][]{
                {new String[]{}, 1, Map.of(2.0, 1.0), Map.of(4.0, 1.0, 2.0, 1.0)},
                {new String[]{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=100"}, 1, Map.of(2.0, 1.0), Map.of(4.0, 1.0, 2.0, 1.0)},
                {new String[]{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=99"}, 0, noOptical, Map.of(5.0, 1.0, 2.0, 1.0)},
                // o1, o2 and o3 form one cluster, and q1 and q2 another.
                {new String[]{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=101"}, 3, Map.of(3.0, 1.0, 2.0, 1.0), Map.of(3.0, 1.0, 1.0, 1.0)},
        };
    }

    @Test(dataProvider = "opticalSettings")
    public void testOpticalDuplicates(final String[] args, final int expectedOpticalDuplicates,
                                      final Map<Double, Double> opticalSets, final Map<Double, Double> nonOpticalSets) throws IOException {
        final Result result = run(opticalInput(true), args);
        assertOpticalInputDuplicates(result, expectedOpticalDuplicates);
        Assert.assertEquals(result.histograms.get("all_sets"), Map.of(5.0, 1.0, 2.0, 1.0));
        // A histogram without counts is not written.
        Assert.assertEquals(result.histograms.getOrDefault("optical_sets", Map.of()), opticalSets);
        Assert.assertEquals(result.histograms.get("non_optical_sets"), nonOpticalSets);
    }

    @Test
    public void testOpticalDuplicatesWithoutReadNameRegex() throws IOException {
        final Result result = run(opticalInput(true), "READ_NAME_REGEX=null");
        assertOpticalInputDuplicates(result, 0);
    }

    private static void assertOpticalInputDuplicates(final Result result, final int expectedOpticalDuplicates) {
        assertDuplicates(result.records, opticalDuplicates(true));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 7, 0, 0, 0, 5, expectedOpticalDuplicates);
    }

    @Test
    public void testLibrariesAreMarkedSeparately() throws IOException {
        final Result result = run(writeCoordinateSorted(multiLibraryInput()));
        assertMultiLibraryResult(result.records, result.metrics);
    }

    @Test
    public void testPairsWithMatesOnOtherContigs() throws IOException {
        final Result result = run(writeCoordinateSorted(interContigInput()));
        assertDuplicates(result.records, List.of("x2/1", "x2/2"));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 3, 0, 0, 0, 1, 0);
    }

    @Test
    public void testCommentsAreAddedToHeader() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        final Result result = run(writeCoordinateSorted(builder), "COMMENT=first note", "COMMENT=second");
        Assert.assertEquals(result.header.getComments(), List.of("@CO\tfirst note", "@CO\tsecond"));
    }

    @Test
    public void testMultipleInputsAreMarkedTogether() throws IOException {
        final SAMRecordSetBuilder first = newBuilder();
        useReadGroup(first, "rgA", "libA");
        addPair(first, "p1", 0, 1000, 1200);
        addFragment(first, "f2", 0, 3000, false, "50M");
        final SAMRecordSetBuilder second = newBuilder();
        useReadGroup(second, "rgA", "libA");
        addPair(second, "p2", 0, 1000, 1200);
        addFragment(second, "f1", 0, 3000, false, "50M");

        final Result result = run(List.of(writeCoordinateSorted(first), writeCoordinateSorted(second)));
        Assert.assertEquals(result.records.size(), 6);
        assertDuplicates(result.records, List.of("p2/1", "p2/2", "f2/0"));
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 2, 2, 0, 0, 1, 1, 0);
    }

    @Test
    public void testEmptyInput() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");

        final Result result = run(writeCoordinateSorted(builder));
        Assert.assertTrue(result.records.isEmpty());
        Assert.assertEquals(result.metrics.size(), 1);
        assertDuplicationMetrics(result.metrics.get("libA"), "libA", 0, 0, 0, 0, 0, 0, 0);
    }

    /** Records tagged with an intermediate program, with the last program in the chain, and with no program. */
    private static File programRecordInput() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        builder.setProgramRecord(new SAMProgramRecord("aligner"));
        addPair(builder, "pgInter", 0, 1000, 1200);
        final SAMProgramRecord sorter = new SAMProgramRecord("sorter");
        sorter.setPreviousProgramGroupId("aligner");
        builder.setProgramRecord(sorter);
        addPair(builder, "pgLeaf", 0, 1000, 1200);
        builder.setProgramRecord(null);
        addPair(builder, "pgNone", 0, 1000, 1200);
        return writeCoordinateSorted(builder);
    }

    @DataProvider(name = "programRecordSettings")
    public Object[][] programRecordSettings() {
        return new Object[][]{
                // args, the expected new program ID (or a prefix for a collision), its name, whether reads are retagged
                {new String[]{}, "MarkDuplicates", "MarkDuplicatesWithMateCigar", true},
                {new String[]{"PROGRAM_RECORD_ID=MDMC", "PROGRAM_GROUP_NAME=marker"}, "MDMC", "marker", true},
                {new String[]{"PROGRAM_RECORD_ID=sorter"}, "sorter.", "MarkDuplicatesWithMateCigar", true},
                {new String[]{"ADD_PG_TAG_TO_READS=false"}, "MarkDuplicates", "MarkDuplicatesWithMateCigar", false},
                {new String[]{"PROGRAM_RECORD_ID=null"}, null, null, false},
        };
    }

    @Test(dataProvider = "programRecordSettings")
    public void testProgramRecordChaining(final String[] args, final String expectedId, final String expectedName,
                                          final boolean retagged) throws IOException {
        final Result result = run(programRecordInput(), args);

        final List<SAMProgramRecord> added = result.header.getProgramRecords().stream()
                .filter(pg -> !pg.getId().equals("aligner") && !pg.getId().equals("sorter"))
                .collect(Collectors.toList());
        final String newId;
        if (expectedId == null) {
            Assert.assertTrue(added.isEmpty());
            newId = null;
        } else {
            Assert.assertEquals(added.size(), 1);
            final SAMProgramRecord program = added.get(0);
            newId = program.getId();
            if (expectedId.endsWith(".")) {
                Assert.assertTrue(newId.startsWith(expectedId) && newId.length() > expectedId.length(), newId);
            } else {
                Assert.assertEquals(newId, expectedId);
            }
            Assert.assertEquals(program.getPreviousProgramGroupId(), "sorter");
            Assert.assertEquals(program.getProgramName(), expectedName);
        }

        for (final SAMRecord record : result.records) {
            final Object pg = record.getAttribute(SAMTag.PG.name());
            switch (record.getReadName()) {
                case "pgInter":
                    Assert.assertEquals(pg, "aligner");
                    break;
                case "pgLeaf":
                    Assert.assertEquals(pg, retagged ? newId : "sorter");
                    break;
                default:
                    Assert.assertNull(pg);
            }
        }
        assertDuplicates(result.records, List.of("pgLeaf/1", "pgLeaf/2", "pgNone/1", "pgNone/2"));
    }

    @Test
    public void testAssumeSortOrderOverridesHeader() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        addPair(builder, "p2", 0, 1000, 1200);
        final List<SAMRecord> records = new ArrayList<>(builder.getRecords());
        records.sort(new SAMRecordCoordinateComparator());
        final File input = write(builder.getHeader(), SAMFileHeader.SortOrder.unsorted, records);

        final Result result = run(input, "ASSUME_SORT_ORDER=coordinate");
        Assert.assertEquals(result.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        assertDuplicates(result.records, List.of("p2/1", "p2/2"));
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*coordinate.*")
    public void testRejectsQueryNameSortedInput() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        run(writeQueryNameSorted(builder.getHeader(), builder.getRecords()));
    }

    @Test(expectedExceptions = PicardException.class, expectedExceptionsMessageRegExp = "(?s).*SUM_OF_BASE_QUALITIES.*")
    public void testRejectsSumOfBaseQualities() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        run(writeCoordinateSorted(builder), "DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES");
    }
}
