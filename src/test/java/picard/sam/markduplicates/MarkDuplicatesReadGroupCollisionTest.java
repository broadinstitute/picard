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
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.stream.Collectors;

/**
 * Tests {@link MarkDuplicates} on two INPUT files whose read groups have the same ID. When the two read groups differ
 * in any attribute they are distinct read groups of the merged input: the output header holds both, the second under a
 * new ID, the records of the second file carry that new ID, and their library, their mates and their optical duplicates
 * follow their own read group. When the two read groups are identical they are one read group of the merged input.
 *
 * Every pair has its 5' ends at chr1:100 forward and chr1:300 reverse; fragments are forward at chr1:2000.
 */
public class MarkDuplicatesReadGroupCollisionTest {

    private static MarkDuplicatesEndToEndRunner file(final String library, final String sample) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner().readGroup("rg", library);
        input.header.getReadGroup("rg").setSample(sample);
        return input;
    }

    private static void addPair(final MarkDuplicatesEndToEndRunner input, final String name, final int quality) {
        input.pair(name, "rg", 0, 100, false, "10M", 0, 291, true, "10M", quality);
    }

    private static MarkDuplicatesEndToEndRunner.Result run(final MarkDuplicatesEndToEndRunner first,
                                                           final MarkDuplicatesEndToEndRunner second, final String... args) {
        final File dir = IOUtil.createTempDir("ReadGroupCollision").toFile();
        try {
            final List<File> inputs = Arrays.asList(first.write(new File(dir, "first.sam")), second.write(new File(dir, "second.sam")));
            final MarkDuplicatesEndToEndRunner.Result result = MarkDuplicatesEndToEndRunner.runOn(dir, inputs, args);
            Assert.assertEquals(result.exitCode, 0);
            return result;
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
    }

    /** ID, library and sample of each read group of the header, in header order. */
    private static List<String> readGroups(final SAMFileHeader header) {
        return header.getReadGroups().stream()
                .map(rg -> rg.getId() + " " + rg.getLibrary() + " " + rg.getSample())
                .collect(Collectors.toList());
    }

    /**
     * For each output record, keyed by its read group ID and its key: its duplicate flag and DT tag ("-" when it is
     * not a duplicate).
     */
    private static Map<String, String> flagsByReadGroup(final MarkDuplicatesEndToEndRunner.Result result) {
        final Map<String, String> flags = new TreeMap<>();
        for (final SAMRecord rec : result.records) {
            final String key = rec.getStringAttribute("RG") + " " + MarkDuplicatesEndToEndRunner.key(rec);
            final Object dt = rec.getAttribute(MarkDuplicates.DUPLICATE_TYPE_TAG);
            Assert.assertNull(flags.put(key, rec.getDuplicateReadFlag() ? String.valueOf(dt) : "-"), key);
        }
        return flags;
    }

    /** The ID the merged header gives the second file's read group, after checking it is not "rg". */
    private static String secondId(final SAMFileHeader header) {
        Assert.assertEquals(header.getReadGroups().size(), 2);
        Assert.assertEquals(header.getReadGroups().get(0).getId(), "rg");
        final String id = header.getReadGroups().get(1).getId();
        Assert.assertNotEquals(id, "rg");
        return id;
    }

    /**
     * The read groups are in libraries libA and libB: nothing is a duplicate across files, and pair "p1", present in
     * both files, is matched with its own mate in each.
     */
    @Test
    public void testCollidingReadGroupsOfDifferentLibraries() {
        final MarkDuplicatesEndToEndRunner first = file("libA", "sample");
        addPair(first, "p1", 30);
        addPair(first, "p2", 20);
        first.fragment("f1", "rg", 0, 2000, false, "10M", 30);
        final MarkDuplicatesEndToEndRunner second = file("libB", "sample");
        addPair(second, "p1", 20);
        addPair(second, "p3", 30);
        second.fragment("f2", "rg", 0, 2000, false, "10M", 30);

        final MarkDuplicatesEndToEndRunner.Result result = run(first, second, "TAGGING_POLICY=All");

        final String id = secondId(result.header);
        Assert.assertEquals(readGroups(result.header), List.of("rg libA sample", id + " libB sample"));
        final Map<String, String> expected = new TreeMap<>();
        expected.put("rg p1/1", "-");
        expected.put("rg p1/2", "-");
        expected.put("rg p2/1", "LB");
        expected.put("rg p2/2", "LB");
        expected.put("rg f1", "-");
        expected.put(id + " p1/1", "LB");
        expected.put(id + " p1/2", "LB");
        expected.put(id + " p3/1", "-");
        expected.put(id + " p3/2", "-");
        expected.put(id + " f2", "-");
        Assert.assertEquals(flagsByReadGroup(result), expected);
        MarkDuplicatesEndToEndRunner.assertMetrics(result.metrics("libA"), "libA", 1, 2, 0, 0, 0, 1, 0);
        MarkDuplicatesEndToEndRunner.assertMetrics(result.metrics("libB"), "libB", 1, 2, 0, 0, 0, 1, 0);
    }

    /**
     * The read groups differ only in their sample and share library libA, so pairs are duplicates across files. The
     * two "R:1:1101:1000:1000" pairs share a name and a location but not a read group, so they are not optical
     * duplicates of each other; "S:1:1101:1050:1000" is an optical duplicate of the first file's pair.
     */
    @Test
    public void testCollidingReadGroupsOfOneLibrary() {
        final MarkDuplicatesEndToEndRunner first = file("libA", "sample");
        addPair(first, "R:1:1101:1000:1000", 30);
        addPair(first, "S:1:1101:1050:1000", 15);
        final MarkDuplicatesEndToEndRunner second = file("libA", "other");
        addPair(second, "R:1:1101:1000:1000", 20);

        final MarkDuplicatesEndToEndRunner.Result result = run(first, second, "TAGGING_POLICY=All");

        final String id = secondId(result.header);
        Assert.assertEquals(readGroups(result.header), List.of("rg libA sample", id + " libA other"));
        final Map<String, String> expected = new TreeMap<>();
        expected.put("rg R:1:1101:1000:1000/1", "-");
        expected.put("rg R:1:1101:1000:1000/2", "-");
        expected.put("rg S:1:1101:1050:1000/1", "SQ");
        expected.put("rg S:1:1101:1050:1000/2", "SQ");
        expected.put(id + " R:1:1101:1000:1000/1", "LB");
        expected.put(id + " R:1:1101:1000:1000/2", "LB");
        Assert.assertEquals(flagsByReadGroup(result), expected);
        MarkDuplicatesEndToEndRunner.assertMetrics(result.metrics("libA"), "libA", 0, 3, 0, 0, 0, 2, 1);
    }

    /**
     * Identical read groups are merged into one, so a pair of the second file within 100 pixels of a pair of the first
     * is an optical duplicate of it.
     */
    @Test
    public void testIdenticalReadGroupsAreMerged() {
        final MarkDuplicatesEndToEndRunner first = file("libA", "sample");
        addPair(first, "R:1:1101:1000:1000", 30);
        final MarkDuplicatesEndToEndRunner second = file("libA", "sample");
        addPair(second, "S:1:1101:1050:1000", 20);

        final MarkDuplicatesEndToEndRunner.Result result = run(first, second, "TAGGING_POLICY=All");

        Assert.assertEquals(readGroups(result.header), List.of("rg libA sample"));
        final Map<String, String> expected = new TreeMap<>();
        expected.put("rg R:1:1101:1000:1000/1", "-");
        expected.put("rg R:1:1101:1000:1000/2", "-");
        expected.put("rg S:1:1101:1050:1000/1", "SQ");
        expected.put("rg S:1:1101:1050:1000/2", "SQ");
        Assert.assertEquals(flagsByReadGroup(result), expected);
        MarkDuplicatesEndToEndRunner.assertMetrics(result.metrics("libA"), "libA", 0, 2, 0, 0, 0, 1, 1);
    }
}
