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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests of {@link CheckDuplicateMarking}: for each MODE, which records sharing a query name are compared, the exit
 * status, and the query names written to OUTPUT, on queryname-sorted input read in file order and on coordinate-sorted
 * input that the program sorts, including a sort that spills to disk.
 */
public class CheckDuplicateMarkingBehaviorTest {

    /**
     * Paired records only, so that every mode can read them. Each group of records shares a query name:
     * consistent groups, a pair whose ends disagree, a secondary or supplementary record that disagrees with its
     * primaries, a pair whose unmapped end disagrees with its mapped end, and a pair that is not a proper pair
     * whose ends disagree.
     */
    private static List<SAMRecord> pairedRecords(final SAMRecordSetBuilder builder) {
        useReadGroup(builder, "rgA", "libA");
        final List<SAMRecord> records = new ArrayList<>();

        records.addAll(addPair(builder, "goodDup", 0, 1000, 1200));
        records.forEach(r -> r.setDuplicateReadFlag(true));
        records.addAll(addPair(builder, "goodNotDup", 0, 1000, 1200));

        final List<SAMRecord> badPair = addPair(builder, "badPair", 0, 2000, 2200);
        badPair.get(0).setDuplicateReadFlag(true);
        records.addAll(badPair);

        final List<SAMRecord> badSecondary = addPair(builder, "badSecondary", 0, 3000, 3200);
        final SAMRecord secondary = badSecondary.get(0).deepCopy();
        secondary.setSecondaryAlignment(true);
        secondary.setAlignmentStart(9000);
        secondary.setDuplicateReadFlag(true);
        records.addAll(badSecondary);
        records.add(secondary);

        final List<SAMRecord> badSupplementary = addPair(builder, "badSupplementary", 0, 4000, 4200);
        final SAMRecord supplementary = badSupplementary.get(1).deepCopy();
        supplementary.setSupplementaryAlignmentFlag(true);
        supplementary.setAlignmentStart(9500);
        supplementary.setDuplicateReadFlag(true);
        records.addAll(badSupplementary);
        records.add(supplementary);

        final List<SAMRecord> badUnmapped = addPairWithUnmappedMate(builder, "badUnmapped", 0, 5000);
        badUnmapped.get(0).setDuplicateReadFlag(true);
        records.addAll(badUnmapped);

        final List<SAMRecord> badImproper = addPair(builder, "badImproper", 0, 6000, 6200);
        badImproper.forEach(r -> r.setProperPairFlag(false));
        badImproper.get(1).setDuplicateReadFlag(true);
        records.addAll(badImproper);
        return records;
    }

    private static File input(final boolean querynameSorted) throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        final List<SAMRecord> records = pairedRecords(builder);
        return querynameSorted ? writeQueryNameSorted(builder.getHeader(), records) : writeCoordinateSorted(builder.getHeader(), records);
    }

    private static List<String> run(final File input, final int expectedStatus, final String... extraArgs) throws IOException {
        final File output = newTempFile("CheckDuplicateMarkingBehaviorTest", ".txt");
        final List<String> args = new ArrayList<>(List.of("INPUT=" + input.getAbsolutePath(), "OUTPUT=" + output.getAbsolutePath()));
        args.addAll(List.of(extraArgs));
        Assert.assertEquals(new CheckDuplicateMarking().instanceMain(args.toArray(new String[0])), expectedStatus);
        return Files.readAllLines(output.toPath());
    }

    @DataProvider(name = "modes")
    public Object[][] modes() {
        final List<Object[]> rows = new ArrayList<>();
        for (final boolean querynameSorted : new boolean[]{true, false}) {
            // In each disagreeing group only one record differs from the group's first record, so it is reported once.
            rows.add(new Object[]{querynameSorted, CheckDuplicateMarking.Mode.ALL,
                    List.of("badImproper", "badPair", "badSecondary", "badSupplementary", "badUnmapped")});
            rows.add(new Object[]{querynameSorted, CheckDuplicateMarking.Mode.PRIMARY_ONLY,
                    List.of("badImproper", "badPair", "badUnmapped")});
            rows.add(new Object[]{querynameSorted, CheckDuplicateMarking.Mode.PRIMARY_MAPPED_ONLY,
                    List.of("badImproper", "badPair")});
            rows.add(new Object[]{querynameSorted, CheckDuplicateMarking.Mode.PRIMARY_PROPER_PAIR_ONLY,
                    List.of("badPair")});
        }
        return rows.toArray(new Object[0][]);
    }

    @Test(dataProvider = "modes")
    public void testModes(final boolean querynameSorted, final CheckDuplicateMarking.Mode mode,
                          final List<String> expectedBadNames) throws IOException {
        final List<String> badNames = run(input(querynameSorted), 1, "MODE=" + mode.name());
        Assert.assertEquals(badNames, expectedBadNames);
    }

    @Test
    public void testDefaultModeIsAll() throws IOException {
        Assert.assertEquals(run(input(true), 1),
                List.of("badImproper", "badPair", "badSecondary", "badSupplementary", "badUnmapped"));
    }

    @Test
    public void testSortSpilledToDisk() throws IOException {
        Assert.assertEquals(run(input(false), 1, "MAX_RECORDS_IN_RAM=2"),
                List.of("badImproper", "badPair", "badSecondary", "badSupplementary", "badUnmapped"));
    }

    @Test
    public void testConsistentMarkingPasses() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        final List<SAMRecord> records = new ArrayList<>(pairedRecords(builder));
        records.removeIf(r -> r.getReadName().startsWith("bad"));
        for (final boolean querynameSorted : new boolean[]{true, false}) {
            final File input = querynameSorted ? writeQueryNameSorted(builder.getHeader(), records) : writeCoordinateSorted(builder.getHeader(), records);
            for (final CheckDuplicateMarking.Mode mode : CheckDuplicateMarking.Mode.values()) {
                Assert.assertEquals(run(input, 0, "MODE=" + mode.name()), List.of());
            }
        }
    }

    @Test
    public void testEmptyInputPasses() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        Assert.assertEquals(run(writeCoordinateSorted(builder), 0), List.of());
    }

    @Test
    public void testWithoutOutput() throws IOException {
        final File input = input(false);
        Assert.assertEquals(new CheckDuplicateMarking().instanceMain(new String[]{"INPUT=" + input.getAbsolutePath()}), 1);
        Assert.assertEquals(new CheckDuplicateMarking().instanceMain(new String[]{"INPUT=" + input.getAbsolutePath(),
                "MODE=PRIMARY_PROPER_PAIR_ONLY"}), 1);
    }

    @Test
    public void testReportsEveryDisagreementPastTheWarningLimit() throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final List<String> expected = new ArrayList<>();
        for (int i = 0; i < 120; i++) {
            final String name = "bad" + (1000 + i);
            addPair(builder, name, 0, 1000 + 10 * i, 1200 + 10 * i).get(i % 2).setDuplicateReadFlag(true);
            expected.add(name);
        }
        Assert.assertEquals(run(writeCoordinateSorted(builder), 1), expected);
    }

    /**
     * On queryname-sorted input each primary record that disagrees with the first record of its query name is
     * reported, so the name is listed once for each of the two primaries. The documentation does not specify this; the
     * expected value follows the current implementation.
     */
    @Test
    public void testQuerynameInputIsReadInFileOrder() throws IOException {
        // The first record of a query name sets the flag the others are compared with, so a duplicate secondary placed
        // first makes both primaries disagree.
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final List<SAMRecord> pair = addPair(builder, "grp", 0, 1000, 1200);
        final SAMRecord secondary = pair.get(0).deepCopy();
        secondary.setSecondaryAlignment(true);
        secondary.setAlignmentStart(9000);
        secondary.setDuplicateReadFlag(true);
        final File input = write(builder.getHeader(), SAMFileHeader.SortOrder.queryname, List.of(secondary, pair.get(0), pair.get(1)));

        Assert.assertEquals(run(input, 1), List.of("grp", "grp"));
        Assert.assertEquals(run(input, 0, "MODE=PRIMARY_ONLY"), List.of());
    }
}
