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
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.ReadEndsForMarkDuplicatesWithBarcodes;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Hand-written SAM inputs for the query-grouped MarkDuplicates tests, and a runner that returns the output
 * records, header and metrics. Every read is 10 bases long on a 10,000 base chr1 or chr2, so the 5' end of
 * a forward read is its alignment start and that of a reverse read is its alignment start plus 9.
 */
final class MarkDuplicatesQueryNameTestSupport {

    static final int READ_LENGTH = 10;
    static final String BASES = "ACGTACGTAC";

    static final int PAIRED = 0x1;
    static final int PROPER = 0x2;
    static final int UNMAPPED = 0x4;
    static final int MATE_UNMAPPED = 0x8;
    static final int REVERSE = 0x10;
    static final int MATE_REVERSE = 0x20;
    static final int FIRST = 0x40;
    static final int SECOND = 0x80;
    static final int SECONDARY = 0x100;
    static final int DUPLICATE = 0x400;
    static final int SUPPLEMENTARY = 0x800;

    /** Tags MarkDuplicates may add, change or remove on a record; all other content must pass through unchanged. */
    static final List<String> TAGS_MARKDUPLICATES_MAY_CHANGE = Arrays.asList("DT", "DI", "DS", "MI", "PG");

    private MarkDuplicatesQueryNameTestSupport() {}

    /** Builds the text of a SAM file one line at a time, in exactly the order given. */
    static final class SamText {
        private final List<String> headerLines = new ArrayList<>();
        private final List<String> recordLines = new ArrayList<>();

        /** @param hdOrder the sort and group fields of the @HD line, e.g. "SO:queryname" or "SO:unsorted\tGO:query" */
        SamText(final String hdOrder) {
            headerLines.add("@HD\tVN:1.6\t" + hdOrder);
            headerLines.add("@SQ\tSN:chr1\tLN:10000");
            headerLines.add("@SQ\tSN:chr2\tLN:10000");
        }

        SamText readGroup(final String id, final String library) {
            headerLines.add("@RG\tID:" + id + "\tSM:sample\tPL:ILLUMINA" + (library == null ? "" : "\tLB:" + library));
            return this;
        }

        /**
         * Adds one record. A record with the unmapped flag and position 0 is unplaced; one with the unmapped flag
         * and a position is placed at that position. Mate fields are filled in from mateContig and matePos (0 for
         * no mate or an unplaced mate).
         */
        SamText record(final String name, final int flag, final String contig, final int pos, final String cigar,
                       final String mateContig, final int matePos, final char qual, final String... tags) {
            final boolean unmapped = (flag & UNMAPPED) != 0;
            final String rname = pos == 0 ? "*" : contig;
            final String rnext = matePos == 0 ? "*" : (mateContig.equals(rname) ? "=" : mateContig);
            final StringBuilder line = new StringBuilder();
            line.append(name).append('\t').append(flag).append('\t').append(rname).append('\t').append(pos)
                    .append('\t').append(unmapped ? 0 : 60).append('\t').append(unmapped ? "*" : cigar)
                    .append('\t').append(rnext).append('\t').append(matePos).append("\t0\t").append(BASES)
                    .append('\t').append(String.valueOf(qual).repeat(READ_LENGTH));
            for (final String tag : tags) {
                line.append('\t').append(tag);
            }
            recordLines.add(line.toString());
            return this;
        }

        /** Adds an unpaired mapped read. */
        SamText fragment(final String name, final String rg, final int pos, final boolean reverse, final char qual,
                         final String... tags) {
            return record(name, reverse ? REVERSE : 0, "chr1", pos, "10M", "*", 0, qual, withRg(rg, tags));
        }

        /** Adds read one then read two of a pair with both ends mapped to chr1, each with 10M and an MC tag. */
        SamText pair(final String name, final String rg, final int pos1, final boolean reverse1, final int pos2,
                     final boolean reverse2, final char qual, final String... tags) {
            return pair(name, rg, pos1, reverse1, pos2, reverse2, qual, qual, tags);
        }

        SamText pair(final String name, final String rg, final int pos1, final boolean reverse1, final int pos2,
                     final boolean reverse2, final char qual1, final char qual2, final String... tags) {
            return pair(name, rg, pos1, reverse1, pos2, reverse2, qual1, qual2, false, tags, tags);
        }

        /** Adds a pair with separate tags on each read, writing read two first when readTwoFirst is set. */
        SamText pair(final String name, final String rg, final int pos1, final boolean reverse1, final int pos2,
                     final boolean reverse2, final char qual1, final char qual2, final boolean readTwoFirst,
                     final String[] tags1, final String[] tags2) {
            if (readTwoFirst) {
                readTwo(name, rg, pos1, reverse1, pos2, reverse2, qual2, tags2);
                readOne(name, rg, pos1, reverse1, pos2, reverse2, qual1, tags1);
            } else {
                readOne(name, rg, pos1, reverse1, pos2, reverse2, qual1, tags1);
                readTwo(name, rg, pos1, reverse1, pos2, reverse2, qual2, tags2);
            }
            return this;
        }

        private void readOne(final String name, final String rg, final int pos1, final boolean reverse1, final int pos2,
                             final boolean reverse2, final char qual, final String[] tags) {
            record(name, PAIRED | PROPER | FIRST | (reverse1 ? REVERSE : 0) | (reverse2 ? MATE_REVERSE : 0),
                    "chr1", pos1, "10M", "chr1", pos2, qual, withRg(rg, append(tags, "MC:Z:10M")));
        }

        private void readTwo(final String name, final String rg, final int pos1, final boolean reverse1, final int pos2,
                             final boolean reverse2, final char qual, final String[] tags) {
            record(name, PAIRED | PROPER | SECOND | (reverse2 ? REVERSE : 0) | (reverse1 ? MATE_REVERSE : 0),
                    "chr1", pos2, "10M", "chr1", pos1, qual, withRg(rg, append(tags, "MC:Z:10M")));
        }

        /**
         * Adds a pair with one read mapped at pos (reverse if requested) and the other unmapped and placed at pos.
         * With mappedFirst, read one is the mapped read; otherwise read one is the unmapped read. Read one is
         * written first either way.
         */
        SamText pairWithUnmappedMate(final String name, final String rg, final int pos, final boolean reverse,
                                     final boolean mappedFirst, final char qual, final String... tags) {
            if (mappedFirst) {
                record(name, PAIRED | FIRST | MATE_UNMAPPED | (reverse ? REVERSE : 0), "chr1", pos, "10M",
                        "chr1", pos, qual, withRg(rg, tags));
                record(name, PAIRED | SECOND | UNMAPPED | (reverse ? MATE_REVERSE : 0), "chr1", pos, "*",
                        "chr1", pos, qual, withRg(rg, tags));
            } else {
                record(name, PAIRED | FIRST | UNMAPPED | (reverse ? MATE_REVERSE : 0), "chr1", pos, "*",
                        "chr1", pos, qual, withRg(rg, tags));
                record(name, PAIRED | SECOND | MATE_UNMAPPED | (reverse ? REVERSE : 0), "chr1", pos, "10M",
                        "chr1", pos, qual, withRg(rg, tags));
            }
            return this;
        }

        List<String> recordLines() {
            return Collections.unmodifiableList(recordLines);
        }

        File write() throws IOException {
            final File file = File.createTempFile("mdQueryName.", ".sam");
            file.deleteOnExit();
            final List<String> all = new ArrayList<>(headerLines);
            all.addAll(recordLines);
            Files.write(file.toPath(), (String.join("\n", all) + "\n").getBytes(StandardCharsets.UTF_8));
            return file;
        }

        private static String[] withRg(final String rg, final String[] tags) {
            return rg == null ? tags : append(tags, "RG:Z:" + rg);
        }

        private static String[] append(final String[] tags, final String tag) {
            final String[] out = Arrays.copyOf(tags, tags.length + 1);
            out[tags.length] = tag;
            return out;
        }
    }

    /** The output of one MarkDuplicates run. */
    static final class Result {
        final SAMFileHeader header;
        final List<SAMRecord> records;
        final List<String> recordLines;
        final List<DuplicationMetrics> metrics;

        Result(final SAMFileHeader header, final List<SAMRecord> records, final List<String> recordLines,
               final List<DuplicationMetrics> metrics) {
            this.header = header;
            this.records = records;
            this.recordLines = recordLines;
            this.metrics = metrics;
        }

        /** Values of the given tag on the output records, in output order (null where absent). */
        List<Object> tagValues(final String tag) {
            final List<Object> values = new ArrayList<>();
            records.forEach(r -> values.add(r.getAttribute(tag)));
            return values;
        }

        DuplicationMetrics metricsFor(final String library) {
            for (final DuplicationMetrics m : metrics) {
                if (library.equals(m.LIBRARY)) {
                    return m;
                }
            }
            throw new AssertionError("No metrics for library " + library);
        }
    }

    /** Runs MarkDuplicates on input with the given extra arguments, asserting that it returns 0. */
    static Result run(final File input, final String... args) throws IOException {
        final File dir = Files.createTempDirectory("mdQueryName").toFile();
        try {
            final File output = new File(dir, "out.sam");
            final File metricsFile = new File(dir, "metrics.txt");
            final List<String> argv = new ArrayList<>();
            argv.add("INPUT=" + input.getAbsolutePath());
            argv.add("OUTPUT=" + output.getAbsolutePath());
            argv.add("METRICS_FILE=" + metricsFile.getAbsolutePath());
            argv.add("TMP_DIR=" + dir.getAbsolutePath());
            argv.addAll(Arrays.asList(args));
            Assert.assertEquals(new MarkDuplicates().instanceMain(argv.toArray(new String[0])), 0);
            return read(output, metricsFile);
        } finally {
            IOUtil.recursiveDelete(dir.toPath());
        }
    }

    /** Runs MarkDuplicates on several inputs. */
    static Result run(final List<File> inputs, final String... args) throws IOException {
        final List<String> argv = new ArrayList<>();
        for (int i = 1; i < inputs.size(); i++) {
            argv.add("INPUT=" + inputs.get(i).getAbsolutePath());
        }
        argv.addAll(Arrays.asList(args));
        return run(inputs.get(0), argv.toArray(new String[0]));
    }

    /** Runs MarkDuplicates expecting it to throw, and returns what it threw. */
    static Throwable runExpectingFailure(final File input, final String... args) throws IOException {
        final File dir = Files.createTempDirectory("mdQueryName").toFile();
        try {
            final List<String> argv = new ArrayList<>();
            argv.add("INPUT=" + input.getAbsolutePath());
            argv.add("OUTPUT=" + new File(dir, "out.sam").getAbsolutePath());
            argv.add("METRICS_FILE=" + new File(dir, "metrics.txt").getAbsolutePath());
            argv.add("TMP_DIR=" + dir.getAbsolutePath());
            argv.addAll(Arrays.asList(args));
            try {
                new MarkDuplicates().instanceMain(argv.toArray(new String[0]));
            } catch (final Throwable t) {
                return t;
            }
            throw new AssertionError("MarkDuplicates did not throw");
        } finally {
            IOUtil.recursiveDelete(dir.toPath());
        }
    }

    private static Result read(final File output, final File metricsFile) throws IOException {
        final List<SAMRecord> records = new ArrayList<>();
        final SAMFileHeader header;
        try (SamReader reader = SamReaderFactory.makeDefault().open(output)) {
            header = reader.getFileHeader();
            reader.iterator().forEachRemaining(records::add);
        }
        final List<String> lines = new ArrayList<>();
        for (final String line : Files.readAllLines(output.toPath(), StandardCharsets.UTF_8)) {
            if (!line.startsWith("@")) {
                lines.add(line);
            }
        }
        final MetricsFile<DuplicationMetrics, Double> metrics = new MetricsFile<>();
        try (Reader reader = Files.newBufferedReader(metricsFile.toPath(), StandardCharsets.UTF_8)) {
            metrics.read(reader);
        }
        return new Result(header, records, lines, metrics.getMetrics());
    }

    /**
     * Asserts that each output line equals the input line at the same position once the duplicate flag bit is
     * cleared and the tags MarkDuplicates may change are dropped, which checks order, every field and every
     * other tag.
     */
    static void assertOnlyDuplicateFlagAndTagsChanged(final List<String> inputLines, final List<String> outputLines) {
        Assert.assertEquals(outputLines.size(), inputLines.size(), "record count");
        for (int i = 0; i < inputLines.size(); i++) {
            Assert.assertEquals(normalize(outputLines.get(i)), normalize(inputLines.get(i)), "record " + i);
        }
    }

    // htsjdk writes tags in its own order, so tags are compared as a sorted list.
    private static List<String> normalize(final String line) {
        final String[] fields = line.split("\t");
        final List<String> out = new ArrayList<>();
        final List<String> tags = new ArrayList<>();
        for (int i = 0; i < fields.length; i++) {
            if (i == 1) {
                out.add(String.valueOf(Integer.parseInt(fields[i]) & ~DUPLICATE));
            } else if (i < 11) {
                out.add(fields[i]);
            } else if (!TAGS_MARKDUPLICATES_MAY_CHANGE.contains(fields[i].substring(0, 2))) {
                tags.add(fields[i]);
            }
        }
        Collections.sort(tags);
        out.addAll(tags);
        return out;
    }

    /**
     * Asserts the duplicate flag of every output record, given as one character per record in output order:
     * 'D' for a duplicate and '.' otherwise.
     */
    static void assertDuplicateFlags(final Result result, final String expected) {
        final StringBuilder observed = new StringBuilder();
        result.records.forEach(r -> observed.append(r.getDuplicateReadFlag() ? 'D' : '.'));
        Assert.assertEquals(observed.toString(), expected);
    }

    /**
     * Asserts the counts of one library's metrics, and that PERCENT_DUPLICATION and ESTIMATED_LIBRARY_SIZE follow
     * from them as documented in DuplicationMetrics.
     */
    static void assertMetrics(final DuplicationMetrics m, final long unpaired, final long pairs,
                              final long secondaryOrSupplementary, final long unmapped, final long unpairedDuplicates,
                              final long pairDuplicates, final long opticalPairDuplicates) {
        Assert.assertEquals(m.UNPAIRED_READS_EXAMINED, unpaired, "UNPAIRED_READS_EXAMINED");
        Assert.assertEquals(m.READ_PAIRS_EXAMINED, pairs, "READ_PAIRS_EXAMINED");
        Assert.assertEquals(m.SECONDARY_OR_SUPPLEMENTARY_RDS, secondaryOrSupplementary, "SECONDARY_OR_SUPPLEMENTARY_RDS");
        Assert.assertEquals(m.UNMAPPED_READS, unmapped, "UNMAPPED_READS");
        Assert.assertEquals(m.UNPAIRED_READ_DUPLICATES, unpairedDuplicates, "UNPAIRED_READ_DUPLICATES");
        Assert.assertEquals(m.READ_PAIR_DUPLICATES, pairDuplicates, "READ_PAIR_DUPLICATES");
        Assert.assertEquals(m.READ_PAIR_OPTICAL_DUPLICATES, opticalPairDuplicates, "READ_PAIR_OPTICAL_DUPLICATES");
        final double examined = unpaired + 2.0 * pairs;
        final double percent = examined == 0 ? 0 : (unpairedDuplicates + 2.0 * pairDuplicates) / examined;
        Assert.assertEquals(m.PERCENT_DUPLICATION, percent, 1e-6, "PERCENT_DUPLICATION");
        Assert.assertEquals(m.ESTIMATED_LIBRARY_SIZE,
                expectedLibrarySize(pairs - opticalPairDuplicates, pairs - pairDuplicates), "ESTIMATED_LIBRARY_SIZE");
    }

    /**
     * The library size X solving the Lander-Waterman equation unique / X = 1 - exp(-pairs / X), truncated, or null
     * when there are no pairs or no duplicate pairs.
     */
    static Long expectedLibrarySize(final long pairs, final long unique) {
        if (pairs <= 0 || pairs - unique <= 0) {
            return null;
        }
        double lo = unique;
        double hi = unique;
        while (landerWaterman(hi, unique, pairs) > 0) {
            hi *= 2;
        }
        for (int i = 0; i < 200; i++) {
            final double mid = (lo + hi) / 2;
            if (landerWaterman(mid, unique, pairs) > 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return (long) Math.floor(lo);
    }

    private static double landerWaterman(final double x, final double unique, final double pairs) {
        return unique / x - 1 + Math.exp(-pairs / x);
    }

    /**
     * Arguments that make the read-end sorting collections hold at most two records in memory, so every input
     * with more than two read ends spills to disk. The ratio is derived from this JVM's heap, which MarkDuplicates
     * shares when run in-process.
     */
    static String[] spillArgs() {
        final double ratio = 2.5 * ReadEndsForMarkDuplicatesWithBarcodes.getSizeOf() / Runtime.getRuntime().maxMemory();
        return new String[]{"SORTING_COLLECTION_SIZE_RATIO=" + ratio, "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1"};
    }

    /** Concatenates argument arrays. */
    static String[] args(final String[] first, final String... rest) {
        final String[] out = Arrays.copyOf(first, first.length + rest.length);
        System.arraycopy(rest, 0, out, first.length, rest.length);
        return out;
    }
}
