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

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * Builds small SAM inputs record by record, runs {@link MarkDuplicates} on them through its command line and reads
 * back the output records, the output header and the metrics file.
 *
 * Records are written in the order (reference index, alignment start), unplaced unmapped reads last; records that
 * tie keep the order in which they were added. Each record is named in assertions by {@link #key(SAMRecord)}.
 */
final class MarkDuplicatesEndToEndRunner {
    static final int CONTIG_LENGTH = 100000;

    final SAMFileHeader header;
    private final List<SAMRecord> records = new ArrayList<>();

    /** A header with sort order coordinate, contigs chr1, chr2 and chr3 and no read groups. */
    MarkDuplicatesEndToEndRunner() {
        this(SAMFileHeader.SortOrder.coordinate);
    }

    MarkDuplicatesEndToEndRunner(final SAMFileHeader.SortOrder sortOrder) {
        header = new SAMFileHeader();
        header.setSortOrder(sortOrder);
        header.setSequenceDictionary(new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("chr1", CONTIG_LENGTH),
                new SAMSequenceRecord("chr2", CONTIG_LENGTH),
                new SAMSequenceRecord("chr3", CONTIG_LENGTH))));
    }

    /** Adds a read group; library may be null for a read group without an LB tag. */
    MarkDuplicatesEndToEndRunner readGroup(final String id, final String library) {
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(id);
        readGroup.setSample("sample");
        readGroup.setPlatform("ILLUMINA");
        if (library != null) {
            readGroup.setLibrary(library);
        }
        header.addReadGroup(readGroup);
        return this;
    }

    private SAMRecord newRecord(final String name, final String readGroup, final int contig, final int start,
                                final boolean negative, final String cigar, final int quality) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName(name);
        rec.setReferenceIndex(contig);
        rec.setAlignmentStart(start);
        rec.setReadNegativeStrandFlag(negative);
        rec.setMappingQuality(60);
        final Cigar parsed = TextCigarCodec.decode(cigar);
        rec.setCigar(parsed);
        fillBasesAndQualities(rec, parsed.getReadLength(), quality);
        if (readGroup != null) {
            rec.setAttribute(SAMTag.RG.name(), readGroup);
        }
        return rec;
    }

    private static void fillBasesAndQualities(final SAMRecord rec, final int length, final int quality) {
        final byte[] bases = new byte[length];
        final byte[] qualities = new byte[length];
        Arrays.fill(bases, (byte) 'A');
        Arrays.fill(qualities, (byte) quality);
        rec.setReadBases(bases);
        rec.setBaseQualities(qualities);
    }

    /** Adds an unpaired mapped read. */
    SAMRecord fragment(final String name, final String readGroup, final int contig, final int start,
                       final boolean negative, final String cigar, final int quality) {
        final SAMRecord rec = newRecord(name, readGroup, contig, start, negative, cigar, quality);
        records.add(rec);
        return rec;
    }

    /** Adds a pair with both ends mapped; returns {first of pair, second of pair}. */
    SAMRecord[] pair(final String name, final String readGroup,
                     final int contig1, final int start1, final boolean negative1, final String cigar1,
                     final int contig2, final int start2, final boolean negative2, final String cigar2,
                     final int quality) {
        final SAMRecord first = newRecord(name, readGroup, contig1, start1, negative1, cigar1, quality);
        final SAMRecord second = newRecord(name, readGroup, contig2, start2, negative2, cigar2, quality);
        first.setReadPairedFlag(true);
        first.setFirstOfPairFlag(true);
        second.setReadPairedFlag(true);
        second.setSecondOfPairFlag(true);
        setMate(first, second);
        setMate(second, first);
        records.add(first);
        records.add(second);
        return new SAMRecord[]{first, second};
    }

    private static void setMate(final SAMRecord rec, final SAMRecord mate) {
        rec.setMateReferenceIndex(mate.getReferenceIndex());
        rec.setMateAlignmentStart(mate.getAlignmentStart());
        rec.setMateNegativeStrandFlag(mate.getReadNegativeStrandFlag());
        rec.setMateUnmappedFlag(mate.getReadUnmappedFlag());
    }

    /**
     * Adds a pair whose second end is unmapped and placed at the position of the first end; returns
     * {mapped first of pair, unmapped second of pair}.
     */
    SAMRecord[] pairWithUnmappedMate(final String name, final String readGroup, final int contig, final int start,
                                     final boolean negative, final String cigar, final int quality) {
        final SAMRecord mapped = newRecord(name, readGroup, contig, start, negative, cigar, quality);
        final SAMRecord unmapped = new SAMRecord(header);
        unmapped.setReadName(name);
        unmapped.setReadUnmappedFlag(true);
        unmapped.setReferenceIndex(contig);
        unmapped.setAlignmentStart(start);
        unmapped.setMappingQuality(0);
        fillBasesAndQualities(unmapped, mapped.getReadLength(), quality);
        if (readGroup != null) {
            unmapped.setAttribute(SAMTag.RG.name(), readGroup);
        }
        mapped.setReadPairedFlag(true);
        mapped.setFirstOfPairFlag(true);
        unmapped.setReadPairedFlag(true);
        unmapped.setSecondOfPairFlag(true);
        setMate(mapped, unmapped);
        setMate(unmapped, mapped);
        records.add(mapped);
        records.add(unmapped);
        return new SAMRecord[]{mapped, unmapped};
    }

    /** Adds an unpaired, unmapped read with no position. */
    SAMRecord unplacedFragment(final String name, final String readGroup, final int quality) {
        final SAMRecord rec = new SAMRecord(header);
        rec.setReadName(name);
        rec.setReadUnmappedFlag(true);
        rec.setMappingQuality(0);
        fillBasesAndQualities(rec, 10, quality);
        if (readGroup != null) {
            rec.setAttribute(SAMTag.RG.name(), readGroup);
        }
        records.add(rec);
        return rec;
    }

    /** Adds a pair with both ends unmapped and no position. */
    SAMRecord[] unplacedPair(final String name, final String readGroup, final int quality) {
        final SAMRecord first = unplacedFragment(name, readGroup, quality);
        final SAMRecord second = unplacedFragment(name, readGroup, quality);
        first.setReadPairedFlag(true);
        first.setFirstOfPairFlag(true);
        first.setMateUnmappedFlag(true);
        second.setReadPairedFlag(true);
        second.setSecondOfPairFlag(true);
        second.setMateUnmappedFlag(true);
        return new SAMRecord[]{first, second};
    }

    /**
     * Adds a secondary (or, if supplementary is true, supplementary) alignment of a read already added, keeping its
     * name, read group, pairing flags and mate information.
     */
    SAMRecord nonPrimary(final SAMRecord primary, final boolean supplementary, final int contig, final int start,
                         final boolean negative, final String cigar, final int quality) {
        final SAMRecord rec = newRecord(primary.getReadName(), primary.getStringAttribute(SAMTag.RG.name()),
                contig, start, negative, cigar, quality);
        rec.setFlags(primary.getFlags());
        rec.setReadNegativeStrandFlag(negative);
        if (supplementary) {
            rec.setSupplementaryAlignmentFlag(true);
        } else {
            rec.setSecondaryAlignment(true);
        }
        if (primary.getReadPairedFlag()) {
            rec.setMateReferenceIndex(primary.getMateReferenceIndex());
            rec.setMateAlignmentStart(primary.getMateAlignmentStart());
        }
        records.add(rec);
        return rec;
    }

    /** Removes a record that was added, e.g. to leave a mate missing. */
    void remove(final SAMRecord rec) {
        Assert.assertTrue(records.remove(rec));
    }

    /** The records in the order they are written to the input file. */
    List<SAMRecord> inputRecords() {
        final List<SAMRecord> sorted = new ArrayList<>(records);
        sorted.sort(Comparator
                .comparingInt((SAMRecord r) -> r.getReferenceIndex() < 0 ? Integer.MAX_VALUE : r.getReferenceIndex())
                .thenComparingInt(SAMRecord::getAlignmentStart));
        return sorted;
    }

    /** Writes the input records to a SAM or BAM file, chosen by the file's extension. */
    File write(final File file) {
        try (final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, true, file, null)) {
            for (final SAMRecord rec : inputRecords()) {
                writer.addAlignment(rec);
            }
        }
        return file;
    }

    /** Writes the input, runs MarkDuplicates on it with the given extra arguments and requires success. */
    Result run(final String... args) {
        final File dir = IOUtil.createTempDir("MarkDuplicatesEndToEnd").toFile();
        try {
            final File input = write(new File(dir, "input.sam"));
            final Result result = runOn(dir, Collections.singletonList(input), args);
            Assert.assertEquals(result.exitCode, 0, "MarkDuplicates failed");
            return result;
        } finally {
            deleteTempDir(dir);
        }
    }

    /**
     * Runs MarkDuplicates on the given inputs, writing its outputs into dir. If it exits with 0 the output records,
     * header and metrics are read back into the result.
     */
    static Result runOn(final File dir, final List<File> inputs, final String... args) {
        return runOnNamedOutput(dir, inputs, "output.sam", args);
    }

    /** As {@link #runOn(File, List, String...)}, writing the output to the named file in dir. */
    static Result runOnNamedOutput(final File dir, final List<File> inputs, final String outputName, final String... args) {
        final File output = new File(dir, outputName);
        final File metricsFile = new File(dir, "metrics.txt");
        final List<String> argv = new ArrayList<>();
        for (final File input : inputs) {
            argv.add("INPUT=" + input.getAbsolutePath());
        }
        argv.add("OUTPUT=" + output.getAbsolutePath());
        argv.add("METRICS_FILE=" + metricsFile.getAbsolutePath());
        argv.addAll(Arrays.asList(args));
        final int exitCode = new MarkDuplicates().instanceMain(argv.toArray(new String[0]));
        if (exitCode != 0) {
            return new Result(exitCode, null, Collections.emptyList(), new MetricsFile<>());
        }

        final List<SAMRecord> outputRecords = new ArrayList<>();
        final SAMFileHeader outputHeader;
        try (final SamReader reader = SamReaderFactory.makeDefault().open(output)) {
            outputHeader = reader.getFileHeader();
            for (final SAMRecord rec : reader) {
                outputRecords.add(rec);
            }
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
        final MetricsFile<DuplicationMetrics, Double> metrics = new MetricsFile<>();
        try (final Reader reader = new FileReader(metricsFile)) {
            metrics.read(reader);
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
        return new Result(exitCode, outputHeader, outputRecords, metrics);
    }

    /**
     * Deletes a temporary directory. MarkDuplicates leaves BAM inputs open, and SAM inputs when it stops with an
     * exception; on Windows those files cannot be deleted, so the directory is then left for deletion on exit.
     */
    static void deleteTempDir(final File dir) {
        if (!IOUtil.deleteDirectoryTree(dir)) {
            IOUtil.deleteOnExit(dir.toPath());
        }
    }

    /** Names a record: the read name, "/1" or "/2" for the ends of a pair, and ":sec" or ":sup" for non-primary records. */
    static String key(final SAMRecord rec) {
        String key = rec.getReadName();
        if (rec.getReadPairedFlag()) {
            key += rec.getFirstOfPairFlag() ? "/1" : "/2";
        }
        if (rec.isSecondaryAlignment()) {
            key += ":sec";
        }
        if (rec.getSupplementaryAlignmentFlag()) {
            key += ":sup";
        }
        return key;
    }

    /** The SAM text of a record without the duplicate flag and without the given tags. */
    static String withoutDuplicateFlagAndTags(final SAMRecord rec, final String... tags) {
        final SAMRecord copy = rec.deepCopy();
        copy.setDuplicateReadFlag(false);
        for (final String tag : tags) {
            copy.setAttribute(tag, null);
        }
        return copy.getSAMString();
    }

    /**
     * Asserts that the output holds exactly the expected input records, in input order, identical apart from the
     * duplicate flag and the given tags.
     */
    static void assertSameRecords(final List<SAMRecord> expected, final List<SAMRecord> actual, final String... tagsThatMayDiffer) {
        final List<String> expectedText = expected.stream()
                .map(r -> withoutDuplicateFlagAndTags(r, tagsThatMayDiffer)).collect(Collectors.toList());
        final List<String> actualText = actual.stream()
                .map(r -> withoutDuplicateFlagAndTags(r, tagsThatMayDiffer)).collect(Collectors.toList());
        Assert.assertEquals(actualText, expectedText);
    }

    /**
     * Expected ESTIMATED_LIBRARY_SIZE from the Lander-Waterman equation C/X = 1 - exp(-N/X), solved here by bisection
     * on X to double precision; null when there are no duplicate pairs.
     */
    static Long expectedLibrarySize(final long readPairs, final long uniqueReadPairs) {
        if (readPairs <= 0 || readPairs - uniqueReadPairs <= 0) {
            return null;
        }
        final double n = readPairs;
        final double c = uniqueReadPairs;
        double lo = c;
        double hi = 2 * c;
        while (c / hi - 1 + Math.exp(-n / hi) > 0) {
            hi *= 2;
        }
        for (int i = 0; i < 200; i++) {
            final double mid = (lo + hi) / 2;
            if (c / mid - 1 + Math.exp(-n / mid) > 0) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        final double root = (lo + hi) / 2;
        final double fraction = root - Math.floor(root);
        Assert.assertTrue(fraction > 1e-6 && fraction < 1 - 1e-6, "library size too close to an integer to test: " + root);
        return (long) Math.floor(root);
    }

    /**
     * Asserts every column of one library's metrics row. PERCENT_DUPLICATION and ESTIMATED_LIBRARY_SIZE are derived
     * here from the counts.
     */
    static void assertMetrics(final DuplicationMetrics metrics, final String library,
                              final long unpairedReadsExamined, final long readPairsExamined,
                              final long secondaryOrSupplementary, final long unmappedReads,
                              final long unpairedReadDuplicates, final long readPairDuplicates,
                              final long readPairOpticalDuplicates) {
        Assert.assertEquals(metrics.LIBRARY, library);
        Assert.assertEquals(metrics.UNPAIRED_READS_EXAMINED, unpairedReadsExamined, "UNPAIRED_READS_EXAMINED");
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, readPairsExamined, "READ_PAIRS_EXAMINED");
        Assert.assertEquals(metrics.SECONDARY_OR_SUPPLEMENTARY_RDS, secondaryOrSupplementary, "SECONDARY_OR_SUPPLEMENTARY_RDS");
        Assert.assertEquals(metrics.UNMAPPED_READS, unmappedReads, "UNMAPPED_READS");
        Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, unpairedReadDuplicates, "UNPAIRED_READ_DUPLICATES");
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES, readPairDuplicates, "READ_PAIR_DUPLICATES");
        Assert.assertEquals(metrics.READ_PAIR_OPTICAL_DUPLICATES, readPairOpticalDuplicates, "READ_PAIR_OPTICAL_DUPLICATES");

        final long examinedReads = unpairedReadsExamined + 2 * readPairsExamined;
        final double percentDuplication = examinedReads == 0 ? 0.0 :
                (unpairedReadDuplicates + 2.0 * readPairDuplicates) / examinedReads;
        Assert.assertEquals(metrics.PERCENT_DUPLICATION, percentDuplication, 1e-6, "PERCENT_DUPLICATION");
        Assert.assertEquals(metrics.ESTIMATED_LIBRARY_SIZE,
                expectedLibrarySize(readPairsExamined - readPairOpticalDuplicates, readPairsExamined - readPairDuplicates),
                "ESTIMATED_LIBRARY_SIZE");
    }

    /** The output of one MarkDuplicates run. */
    static final class Result {
        final int exitCode;
        final SAMFileHeader header;
        final List<SAMRecord> records;
        final MetricsFile<DuplicationMetrics, Double> metricsFile;

        Result(final int exitCode, final SAMFileHeader header, final List<SAMRecord> records,
               final MetricsFile<DuplicationMetrics, Double> metricsFile) {
            this.exitCode = exitCode;
            this.header = header;
            this.records = records;
            this.metricsFile = metricsFile;
        }

        /** The keys of the output records flagged as duplicates, sorted. */
        Set<String> duplicates() {
            return records.stream().filter(SAMRecord::getDuplicateReadFlag).map(MarkDuplicatesEndToEndRunner::key)
                    .collect(Collectors.toCollection(TreeSet::new));
        }

        void assertDuplicates(final String... keys) {
            Assert.assertEquals(duplicates(), new TreeSet<>(Arrays.asList(keys)));
        }

        SAMRecord record(final String key) {
            final List<SAMRecord> matches = records.stream().filter(r -> key(r).equals(key)).collect(Collectors.toList());
            Assert.assertEquals(matches.size(), 1, "records named " + key);
            return matches.get(0);
        }

        List<String> keys() {
            return records.stream().map(MarkDuplicatesEndToEndRunner::key).collect(Collectors.toList());
        }

        DuplicationMetrics metrics(final String library) {
            final List<DuplicationMetrics> matches = metricsFile.getMetrics().stream()
                    .filter(m -> library.equals(m.LIBRARY)).collect(Collectors.toList());
            Assert.assertEquals(matches.size(), 1, "metrics rows for library " + library);
            return matches.get(0);
        }

        /** The only metrics row. */
        DuplicationMetrics onlyMetrics() {
            Assert.assertEquals(metricsFile.getMetrics().size(), 1);
            return metricsFile.getMetrics().get(0);
        }

        /** The histogram whose values are labelled valueLabel, or null if there is none. */
        Histogram<Double> histogram(final String valueLabel) {
            final List<Histogram<Double>> matches = metricsFile.getAllHistograms().stream()
                    .filter(h -> valueLabel.equals(h.getValueLabel())).collect(Collectors.toList());
            Assert.assertTrue(matches.size() <= 1, "histograms labelled " + valueLabel);
            return matches.isEmpty() ? null : matches.get(0);
        }

        Set<String> histogramLabels() {
            return metricsFile.getAllHistograms().stream().map(Histogram::getValueLabel).collect(Collectors.toCollection(HashSet::new));
        }

        /** The bins of the histogram labelled valueLabel whose value is not zero; empty if there is no such histogram. */
        Map<Double, Double> nonZeroBins(final String valueLabel) {
            final Map<Double, Double> bins = new TreeMap<>();
            final Histogram<Double> histogram = histogram(valueLabel);
            if (histogram != null) {
                for (final Histogram.Bin<Double> bin : histogram.values()) {
                    if (bin.getValue() != 0) {
                        bins.put(bin.getId(), bin.getValue());
                    }
                }
            }
            return bins;
        }

        /** Asserts the three duplicate set size histograms, each given as {set size, count, set size, count, ...}. */
        void assertSetSizeHistograms(final double[] allSets, final double[] opticalSets, final double[] nonOpticalSets) {
            Assert.assertEquals(nonZeroBins("all_sets"), bins(allSets), "all_sets");
            Assert.assertEquals(nonZeroBins("optical_sets"), bins(opticalSets), "optical_sets");
            Assert.assertEquals(nonZeroBins("non_optical_sets"), bins(nonOpticalSets), "non_optical_sets");
        }
    }

    static Map<Double, Double> bins(final double... sizeAndCount) {
        final Map<Double, Double> bins = new TreeMap<>();
        for (int i = 0; i < sizeAndCount.length; i += 2) {
            bins.put(sizeAndCount[i], sizeAndCount[i + 1]);
        }
        return bins;
    }
}
