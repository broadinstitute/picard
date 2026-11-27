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
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.markduplicates.MarkDuplicates.DuplicateTaggingPolicy;
import picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.Result;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertMetrics;
import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.assertSameRecords;
import static picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.key;

/**
 * Tests the options of {@link MarkDuplicates} that change what is written rather than which reads are duplicates,
 * on coordinate-sorted input: TAGGING_POLICY and CLEAR_DT, REMOVE_DUPLICATES and REMOVE_SEQUENCING_DUPLICATES, the
 * program record options (PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, PROGRAM_GROUP_COMMAND_LINE, PROGRAM_GROUP_NAME,
 * ADD_PG_TAG_TO_READS), COMMENT, ASSUME_SORTED and ASSUME_SORT_ORDER, and that every other field of every record is
 * written unchanged and in input order.
 */
public class MarkDuplicatesOutputOptionsTest {
    private static final int CHR1 = 0, CHR2 = 1, CHR3 = 2;
    private static final boolean F = false, R = true;
    private static final String DT = MarkDuplicates.DUPLICATE_TYPE_TAG;
    private static final String SQ = MarkDuplicates.DUPLICATE_TYPE_SEQUENCING;
    private static final String LB = MarkDuplicates.DUPLICATE_TYPE_LIBRARY;
    private static final String PG = SAMTag.PG.name();

    // keeper K, optical duplicate O and library duplicate L share their 5' ends; fragment f2 duplicates f1; U is
    // unique. K/1, L/1 and U/1 carry DT tags in the input.
    private static final String K = "K:1:1101:1000:1000";
    private static final String O = "O:1:1101:1050:1050";
    private static final String L = "L:1:1101:5000:5000";
    private static final String F1 = "f1:1:1101:1000:1000";
    private static final String F2 = "f2:1:1101:1001:1001";
    private static final String U = "U:1:1101:9000:9000";

    private static MarkDuplicatesEndToEndRunner mixedInput() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        input.pair(K, "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30)[0].setAttribute(DT, LB);
        input.pair(O, "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.pair(L, "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 10)[0].setAttribute(DT, SQ);
        input.fragment(F1, "rg", CHR2, 100, F, "10M", 30);
        input.fragment(F2, "rg", CHR2, 100, F, "10M", 20);
        input.pair(U, "rg", CHR3, 100, F, "10M", CHR3, 291, R, "10M", 30)[0].setAttribute(DT, SQ);
        return input;
    }

    private static void assertMixedMetrics(final Result result) {
        assertMetrics(result.onlyMetrics(), "lib", 2, 4, 0, 0, 1, 2, 1);
        MarkDuplicatesCoordinateSortedTest.assertCoverageMultHistogram(result, 4, 2, 1);
        result.assertSetSizeHistograms(new double[]{1, 1, 3, 1}, new double[]{2, 1}, new double[]{1, 1, 2, 1});
    }

    private static Map<String, String> tags(final String... keysAndValues) {
        final Map<String, String> map = new TreeMap<>();
        for (int i = 0; i < keysAndValues.length; i += 2) {
            map.put(keysAndValues[i], keysAndValues[i + 1]);
        }
        return map;
    }

    private static Map<String, String> duplicateTypeTags(final Result result) {
        final Map<String, String> map = new TreeMap<>();
        for (final SAMRecord rec : result.records) {
            if (rec.getAttribute(DT) != null) {
                map.put(key(rec), rec.getStringAttribute(DT));
            }
        }
        return map;
    }

    @DataProvider(name = "taggingPolicies")
    public Object[][] taggingPolicies() {
        return new Object[][]{
                {null, null, tags()},
                {DuplicateTaggingPolicy.DontTag, true, tags()},
                {DuplicateTaggingPolicy.DontTag, false, tags(K + "/1", LB, L + "/1", SQ, U + "/1", SQ)},
                {DuplicateTaggingPolicy.OpticalOnly, true, tags(O + "/1", SQ, O + "/2", SQ)},
                {DuplicateTaggingPolicy.OpticalOnly, false, tags(O + "/1", SQ, O + "/2", SQ, K + "/1", LB, L + "/1", SQ, U + "/1", SQ)},
                {DuplicateTaggingPolicy.All, true, tags(O + "/1", SQ, O + "/2", SQ, L + "/1", LB, L + "/2", LB, F2, LB)},
                {DuplicateTaggingPolicy.All, null, tags(O + "/1", SQ, O + "/2", SQ, L + "/1", LB, L + "/2", LB, F2, LB)},
                {DuplicateTaggingPolicy.All, false, tags(O + "/1", SQ, O + "/2", SQ, L + "/1", LB, L + "/2", LB, F2, LB, K + "/1", LB, U + "/1", SQ)},
        };
    }

    /**
     * With CLEAR_DT=false and a policy that does not tag library duplicates (DontTag, OpticalOnly), the library
     * duplicate L/1 keeps the DT:Z:SQ it carries in the input. The documentation does not specify this; the expected
     * value follows the current implementation.
     */
    @Test(dataProvider = "taggingPolicies")
    public void testTaggingPolicyAndClearDt(final DuplicateTaggingPolicy policy, final Boolean clearDt, final Map<String, String> expectedTags) {
        final MarkDuplicatesEndToEndRunner input = mixedInput();
        final List<String> args = new ArrayList<>();
        if (policy != null) {
            args.add("TAGGING_POLICY=" + policy.name());
        }
        if (clearDt != null) {
            args.add("CLEAR_DT=" + clearDt);
        }
        final Result result = input.run(args.toArray(new String[0]));
        assertSameRecords(input.inputRecords(), result.records, PG, DT);
        result.assertDuplicates(O + "/1", O + "/2", L + "/1", L + "/2", F2);
        Assert.assertEquals(duplicateTypeTags(result), expectedTags);
        assertMixedMetrics(result);
    }

    @DataProvider(name = "removals")
    public Object[][] removals() {
        return new Object[][]{
                {false, false, new String[0]},
                {true, false, new String[]{O + "/1", O + "/2", L + "/1", L + "/2", F2}},
                {false, true, new String[]{O + "/1", O + "/2"}},
                {true, true, new String[]{O + "/1", O + "/2", L + "/1", L + "/2", F2}},
        };
    }

    @Test(dataProvider = "removals")
    public void testRemoveDuplicatesAndSequencingDuplicates(final boolean removeDuplicates, final boolean removeSequencingDuplicates,
                                                            final String[] removed) {
        for (final DuplicateTaggingPolicy policy : new DuplicateTaggingPolicy[]{DuplicateTaggingPolicy.DontTag, DuplicateTaggingPolicy.All}) {
            final MarkDuplicatesEndToEndRunner input = mixedInput();
            final Result result = input.run("REMOVE_DUPLICATES=" + removeDuplicates,
                    "REMOVE_SEQUENCING_DUPLICATES=" + removeSequencingDuplicates, "TAGGING_POLICY=" + policy.name());

            final Set<String> removedKeys = new HashSet<>(Arrays.asList(removed));
            final List<SAMRecord> kept = input.inputRecords().stream().filter(r -> !removedKeys.contains(key(r)))
                    .collect(Collectors.toList());
            assertSameRecords(kept, result.records, PG, DT);
            final List<String> stillFlagged = new ArrayList<>(Arrays.asList(O + "/1", O + "/2", L + "/1", L + "/2", F2));
            stillFlagged.removeAll(removedKeys);
            result.assertDuplicates(stillFlagged.toArray(new String[0]));
            // the metrics count the removed records
            assertMixedMetrics(result);
        }
    }

    @Test
    public void testBamInputAndOutput() {
        final MarkDuplicatesEndToEndRunner input = mixedInput();
        final File dir = IOUtil.createTempDir("MarkDuplicatesBam").toFile();
        final Result result;
        try {
            final File bam = input.write(new File(dir, "input.bam"));
            result = MarkDuplicatesEndToEndRunner.runOnNamedOutput(dir, Collections.singletonList(bam), "output.bam", "TAGGING_POLICY=All");
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
        Assert.assertEquals(result.exitCode, 0);
        assertSameRecords(input.inputRecords(), result.records, PG, DT);
        result.assertDuplicates(O + "/1", O + "/2", L + "/1", L + "/2", F2);
        Assert.assertEquals(duplicateTypeTags(result), tags(O + "/1", SQ, O + "/2", SQ, L + "/1", LB, L + "/2", LB, F2, LB));
        assertMixedMetrics(result);
    }

    private static MarkDuplicatesEndToEndRunner inputWithProgramGroups(final String... readProgramGroups) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        int i = 0;
        for (final String pg : readProgramGroups) {
            if (pg != null && input.header.getProgramRecord(pg) == null) {
                final SAMProgramRecord record = new SAMProgramRecord(pg);
                record.setProgramName(pg + "-name");
                input.header.addProgramRecord(record);
            }
            final SAMRecord[] pair = input.pair("p" + i, "rg", CHR1, 100 + 10 * i, F, "10M", CHR1, 291 + 10 * i, R, "10M", 30);
            for (final SAMRecord rec : pair) {
                rec.setAttribute(PG, pg);
            }
            i++;
        }
        return input;
    }

    /** The PG tag of each input record, by key. */
    private static Map<String, String> inputProgramGroups(final MarkDuplicatesEndToEndRunner input) {
        final Map<String, String> map = new HashMap<>();
        for (final SAMRecord rec : input.inputRecords()) {
            map.put(key(rec), rec.getStringAttribute(PG));
        }
        return map;
    }

    private static List<SAMProgramRecord> newProgramRecords(final MarkDuplicatesEndToEndRunner input, final Result result) {
        final Set<String> inputIds = input.header.getProgramRecords().stream().map(SAMProgramRecord::getId).collect(Collectors.toSet());
        return result.header.getProgramRecords().stream().filter(r -> !inputIds.contains(r.getId())).collect(Collectors.toList());
    }

    @Test
    public void testDefaultProgramRecordChainsToTheInputProgramGroup() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups("bwa", "bwa");
        final Result result = input.run();
        assertSameRecords(input.inputRecords(), result.records, PG);

        final List<SAMProgramRecord> added = newProgramRecords(input, result);
        Assert.assertEquals(added.size(), 1);
        final SAMProgramRecord pg = added.get(0);
        Assert.assertEquals(pg.getId(), "MarkDuplicates");
        Assert.assertEquals(pg.getProgramName(), "MarkDuplicates");
        Assert.assertEquals(pg.getPreviousProgramGroupId(), "bwa");
        Assert.assertNotNull(pg.getCommandLine());
        Assert.assertTrue(pg.getCommandLine().contains("output.sam"), pg.getCommandLine());
        Assert.assertEquals(result.header.getProgramRecord("bwa").getProgramName(), "bwa-name");
        for (final SAMRecord rec : result.records) {
            Assert.assertEquals(rec.getStringAttribute(PG), "MarkDuplicates", key(rec));
        }
    }

    @Test
    public void testProgramRecordForReadsWithoutProgramGroup() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups(null, null);
        input.unplacedFragment("unplaced1", "rg", 30);
        input.unplacedFragment("unplaced2", "rg", 30);
        final Result result = input.run();
        assertSameRecords(input.inputRecords(), result.records, PG);

        final List<SAMProgramRecord> added = newProgramRecords(input, result);
        Assert.assertEquals(added.size(), 1);
        Assert.assertEquals(added.get(0).getId(), "MarkDuplicates");
        Assert.assertNull(added.get(0).getPreviousProgramGroupId());
        for (final SAMRecord rec : result.records) {
            Assert.assertEquals(rec.getStringAttribute(PG), "MarkDuplicates", key(rec));
        }
    }

    @Test
    public void testProgramRecordIdThatIsAlreadyTaken() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups("MarkDuplicates");
        final Result result = input.run();

        final List<SAMProgramRecord> added = newProgramRecords(input, result);
        Assert.assertEquals(added.size(), 1);
        final String newId = added.get(0).getId();
        Assert.assertTrue(newId.startsWith("MarkDuplicates."), newId);
        Assert.assertEquals(added.get(0).getPreviousProgramGroupId(), "MarkDuplicates");
        for (final SAMRecord rec : result.records) {
            Assert.assertEquals(rec.getStringAttribute(PG), newId, key(rec));
        }
    }

    @Test
    public void testOneProgramRecordPerProgramGroupOfTheReads() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups("bwa", "other", null, "bwa", null);
        final Map<String, String> before = inputProgramGroups(input);
        final Result result = input.run();
        assertSameRecords(input.inputRecords(), result.records, PG);

        final List<SAMProgramRecord> added = newProgramRecords(input, result);
        Assert.assertEquals(added.size(), 3);
        Assert.assertEquals(added.stream().map(SAMProgramRecord::getPreviousProgramGroupId).collect(Collectors.toSet()),
                new HashSet<>(Arrays.asList("bwa", "other", null)));
        for (final SAMProgramRecord pg : added) {
            Assert.assertTrue(pg.getId().startsWith("MarkDuplicates"), pg.getId());
            Assert.assertEquals(pg.getProgramName(), "MarkDuplicates");
        }
        // each read points at the new record that follows its previous program group
        for (final SAMRecord rec : result.records) {
            final SAMProgramRecord pg = result.header.getProgramRecord(rec.getStringAttribute(PG));
            Assert.assertTrue(added.contains(pg), key(rec));
            Assert.assertEquals(pg.getPreviousProgramGroupId(), before.get(key(rec)), key(rec));
        }
    }

    @Test
    public void testProgramRecordArguments() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups("bwa", null);
        final Result result = input.run("PROGRAM_RECORD_ID=MD", "PROGRAM_GROUP_VERSION=1.2.3",
                "PROGRAM_GROUP_COMMAND_LINE=mark the duplicates", "PROGRAM_GROUP_NAME=Marker");

        final List<SAMProgramRecord> added = newProgramRecords(input, result);
        Assert.assertEquals(added.size(), 2);
        for (final SAMProgramRecord pg : added) {
            Assert.assertTrue(pg.getId().startsWith("MD"), pg.getId());
            Assert.assertEquals(pg.getProgramVersion(), "1.2.3");
            Assert.assertEquals(pg.getCommandLine(), "mark the duplicates");
            Assert.assertEquals(pg.getProgramName(), "Marker");
        }
        Assert.assertTrue(added.stream().anyMatch(pg -> pg.getId().equals("MD")));
        for (final SAMRecord rec : result.records) {
            Assert.assertTrue(rec.getStringAttribute(PG).startsWith("MD"), key(rec));
        }
    }

    @Test
    public void testNullProgramRecordIdAddsNoProgramRecord() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups("bwa", null);
        final Result result = input.run("PROGRAM_RECORD_ID=null");
        Assert.assertEquals(result.header.getProgramRecords(), input.header.getProgramRecords());
        // the PG tags of the reads are left as they were
        assertSameRecords(input.inputRecords(), result.records);
    }

    @Test
    public void testAddPgTagToReadsFalseLeavesReadsAlone() {
        final MarkDuplicatesEndToEndRunner input = inputWithProgramGroups("bwa", null);
        final Result result = input.run("ADD_PG_TAG_TO_READS=false");
        Assert.assertEquals(newProgramRecords(input, result).size(), 2);
        assertSameRecords(input.inputRecords(), result.records);
    }

    @Test
    public void testComments() {
        final MarkDuplicatesEndToEndRunner input = mixedInput();
        input.header.addComment("from the input");
        final Result result = input.run("COMMENT=first", "COMMENT=second comment");
        Assert.assertEquals(result.header.getComments(), Arrays.asList("@CO\tfrom the input", "@CO\tfirst", "@CO\tsecond comment"));
    }

    @Test
    public void testNoComments() {
        final MarkDuplicatesEndToEndRunner input = mixedInput();
        final Result result = input.run();
        Assert.assertEquals(result.header.getComments(), Collections.emptyList());
    }

    @DataProvider(name = "assumedCoordinateOrder")
    public Object[][] assumedCoordinateOrder() {
        return new Object[][]{
                {SAMFileHeader.SortOrder.unsorted, "ASSUME_SORTED=true"},
                {SAMFileHeader.SortOrder.unsorted, "ASSUME_SORT_ORDER=coordinate"},
                {SAMFileHeader.SortOrder.unknown, "ASSUME_SORTED=true"},
                {SAMFileHeader.SortOrder.unknown, "ASSUME_SORT_ORDER=coordinate"},
                {SAMFileHeader.SortOrder.coordinate, "ASSUME_SORT_ORDER=coordinate"},
                {SAMFileHeader.SortOrder.coordinate, "ASSUME_SORTED=false"},
        };
    }

    /** Input in coordinate order is marked as such when the header says so or when told to assume it. */
    @Test(dataProvider = "assumedCoordinateOrder")
    public void testAssumeCoordinateSortOrder(final SAMFileHeader.SortOrder headerOrder, final String arg) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner(headerOrder).readGroup("rg", "lib");
        input.pair("best", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 30);
        input.pair("dup", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20);
        input.pair("toOtherContig", "rg", CHR1, 200, F, "10M", CHR2, 291, R, "10M", 30);
        input.fragment("fragment", "rg", CHR2, 100, F, "10M", 30);
        input.fragment("fragmentDup", "rg", CHR2, 100, F, "10M", 20);
        input.pairWithUnmappedMate("mateUnmapped", "rg", CHR3, 100, F, "10M", 30);
        input.unplacedFragment("unplaced", "rg", 30);

        final Result result = input.run(arg);
        Assert.assertEquals(result.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
        assertSameRecords(input.inputRecords(), result.records, PG);
        result.assertDuplicates("dup/1", "dup/2", "fragmentDup");
        assertMetrics(result.onlyMetrics(), "lib", 3, 3, 0, 2, 1, 1, 0);
    }

    /** Every field and tag of every record is written unchanged apart from the duplicate flag and the PG tag. */
    @Test
    public void testOtherFieldsAndTagsAreUnchanged() {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", "lib");
        input.header.addProgramRecord(new SAMProgramRecord("bwa"));
        final List<SAMRecord[]> pairs = new ArrayList<>();
        pairs.add(input.pair("best", "rg", CHR1, 100, F, "4M1I5M", CHR1, 289, R, "3M2D7M", 30));
        pairs.add(input.pair("dup", "rg", CHR1, 100, F, "10M", CHR1, 291, R, "10M", 20));
        final List<SAMRecord> all = new ArrayList<>();
        for (final SAMRecord[] pair : pairs) {
            for (int end = 0; end < 2; end++) {
                pair[end].setProperPairFlag(true);
                pair[end].setInferredInsertSize(end == 0 ? 201 : -201);
                pair[end].setAttribute("MC", pair[1 - end].getCigarString());
                all.add(pair[end]);
            }
        }
        all.add(input.fragment("fragment", "rg", CHR2, 100, R, "2S6M2S", 25));
        all.addAll(Arrays.asList(input.pairWithUnmappedMate("mateUnmapped", "rg", CHR2, 500, R, "10M", 35)));
        all.add(input.nonPrimary(pairs.get(0)[0], true, CHR3, 700, R, "5H5M", 30));
        all.add(input.unplacedFragment("unplaced", "rg", 12));
        int i = 0;
        for (final SAMRecord rec : all) {
            i++;
            if (!rec.getReadUnmappedFlag()) {
                rec.setMappingQuality(10 + i);
                rec.setAttribute("NM", i);
                rec.setAttribute("MD", "5A" + i);
            }
            rec.setAttribute("XS", 1.5f * i);
            rec.setAttribute("ZB", new int[]{i, -i, 3});
            rec.setAttribute("ZC", 'c');
            rec.setAttribute("OQ", rec.getBaseQualityString());
            rec.setAttribute(PG, "bwa");
            final byte[] bases = rec.getReadBases();
            for (int b = 0; b < bases.length; b++) {
                bases[b] = (byte) "ACGTN".charAt((b + i) % 5);
            }
            rec.setReadBases(bases);
        }

        final Result result = input.run("TAGGING_POLICY=All");
        assertSameRecords(input.inputRecords(), result.records, PG, DT);
        result.assertDuplicates("dup/1", "dup/2");
        for (final SAMRecord rec : result.records) {
            Assert.assertEquals(rec.getStringAttribute(PG), "MarkDuplicates", key(rec));
            Assert.assertEquals(rec.getStringAttribute(DT), rec.getDuplicateReadFlag() ? LB : null, key(rec));
        }
    }
}
