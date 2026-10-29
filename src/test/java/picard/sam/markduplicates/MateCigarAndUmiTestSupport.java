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
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.testng.Assert;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Builds small SAM inputs, and reads back the records and metrics written by the mate-cigar and UMI-aware
 * duplicate-marking programs and by {@link CheckDuplicateMarking} and {@link CollectDuplicateMetrics}.
 */
final class MateCigarAndUmiTestSupport {
    static final int QUALITY = 30;

    private MateCigarAndUmiTestSupport() {
    }

    /** A builder with no default read group; records are kept in insertion order and sorted when written. */
    static SAMRecordSetBuilder newBuilder() {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.coordinate, false);
        builder.setReadLength(50);
        return builder;
    }

    /** Adds a read group with the given library to the builder's header and assigns it to the records added next. */
    static SAMReadGroupRecord useReadGroup(final SAMRecordSetBuilder builder, final String id, final String library) {
        final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(id);
        readGroup.setSample("sample");
        readGroup.setPlatform("ILLUMINA");
        readGroup.setLibrary(library);
        builder.setReadGroup(readGroup);
        return readGroup;
    }

    /** Adds a pair with both ends mapped to the same contig; the MC tags are set from the CIGARs. */
    static List<SAMRecord> addPair(final SAMRecordSetBuilder builder, final String name, final int contig,
                                   final int start1, final boolean negative1, final String cigar1,
                                   final int start2, final boolean negative2, final String cigar2) {
        return builder.addPair(name, contig, contig, start1, start2, false, false, cigar1, cigar2,
                negative1, negative2, false, false, QUALITY);
    }

    /** Adds an FR pair: the first end forward at start1, the second end reverse at start2, both 50M. */
    static List<SAMRecord> addPair(final SAMRecordSetBuilder builder, final String name, final int contig,
                                   final int start1, final int start2) {
        return addPair(builder, name, contig, start1, false, "50M", start2, true, "50M");
    }

    /** Adds a pair whose first end is mapped and whose second end is unmapped and placed at the first end's position. */
    static List<SAMRecord> addPairWithUnmappedMate(final SAMRecordSetBuilder builder, final String name, final int contig, final int start) {
        return builder.addPair(name, contig, contig, start, start, false, true, "50M", null,
                false, false, false, false, QUALITY);
    }

    /** Adds a pair with both ends unmapped and without a position. */
    static List<SAMRecord> addUnmappedPair(final SAMRecordSetBuilder builder, final String name) {
        return builder.addPair(name, -1, -1, -1, -1, true, true, null, null, false, false, false, false, QUALITY);
    }

    static SAMRecord addFragment(final SAMRecordSetBuilder builder, final String name, final int contig, final int start,
                                 final boolean negative, final String cigar) {
        return builder.addFrag(name, contig, start, negative, false, cigar, null, QUALITY);
    }

    static SAMRecord addSecondaryFragment(final SAMRecordSetBuilder builder, final String name, final int contig, final int start) {
        return builder.addFrag(name, contig, start, false, false, "50M", null, QUALITY, true, false);
    }

    static SAMRecord addSupplementaryFragment(final SAMRecordSetBuilder builder, final String name, final int contig, final int start) {
        return builder.addFrag(name, contig, start, false, false, "50M", null, QUALITY, false, true);
    }

    static final String MIXED_LIBRARY = "libA";

    /**
     * A single-library input with: FR pairs where the keeper wins on mapped reference length and the rest are ordered
     * by name, and a pair whose mate 5' end differs by one base; forward fragments that match only after soft clipping,
     * and a reverse fragment with the same 5' coordinate; a fragment and a pair with an unmapped mate that share the 5'
     * end of a mapped pair; secondary and supplementary records beside a primary fragment; pairs on a second contig at
     * the same coordinates as the first group; and an unmapped pair at the end.
     * The duplicates are listed in {@link #MIXED_DUPLICATES}.
     */
    static SAMRecordSetBuilder mixedInput() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", MIXED_LIBRARY);

        addPair(builder, "pA1", 0, 1000, 1200);
        addPair(builder, "pA2", 0, 1000, 1200);
        addPair(builder, "pA3", 0, 1000, false, "50M", 1195, true, "55M");
        addPair(builder, "pA4", 0, 1000, 1201);

        addFragment(builder, "fB1", 0, 3000, false, "50M");
        addFragment(builder, "fB2", 0, 3000, false, "50M");
        addFragment(builder, "fB3", 0, 3005, false, "5S45M");
        addFragment(builder, "fB5", 0, 2951, true, "50M");

        addPair(builder, "pC1", 0, 5000, 5300);
        addFragment(builder, "fC1", 0, 5000, false, "50M");
        addPairWithUnmappedMate(builder, "uC1", 0, 5000);

        addFragment(builder, "fD1", 0, 7000, false, "50M");
        addSecondaryFragment(builder, "sD1", 0, 7000);
        addSupplementaryFragment(builder, "sD2", 0, 7000);

        addPair(builder, "pE1", 1, 1000, 1200);
        addPair(builder, "pE2", 1, 1000, 1200);

        addUnmappedPair(builder, "zZ1");
        return builder;
    }

    static final List<String> MIXED_DUPLICATES = List.of(
            "pA1/1", "pA1/2", "pA2/1", "pA2/2", "fB2/0", "fB3/0", "fC1/0", "uC1/1", "pE2/1", "pE2/2");

    /** The metrics of {@link #mixedInput()} when no read name carries a location. */
    static void assertMixedMetrics(final DuplicationMetrics metrics) {
        assertDuplicationMetrics(metrics, MIXED_LIBRARY, 7, 7, 2, 3, 4, 3, 0);
    }

    /**
     * Four pairs at one position and two at another in library libA, named so that tile and x,y can be parsed, and
     * optionally a fifth pair at the first position in a second read group of the same library. Only o1 and o2 are
     * within 100 pixels of each other in the same tile and read group; o3 is 101 pixels from o2, o4 is in another tile,
     * o5 has o1's location but the other read group, and q2 is 101 pixels from q1. Keepers are the smallest names in
     * each set; the duplicates are listed by {@link #opticalDuplicates(boolean)}.
     */
    static File opticalInput(final boolean withSecondReadGroup) throws IOException {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "RUN:1:1101:1000:1000", 0, 1000, 1200); // o1
        addPair(builder, "RUN:1:1101:1100:1000", 0, 1000, 1200); // o2
        addPair(builder, "RUN:1:1101:1201:1000", 0, 1000, 1200); // o3
        addPair(builder, "RUN:1:1102:1000:1000", 0, 1000, 1200); // o4
        addPair(builder, "RUN:1:1101:5000:5000", 0, 3000, 3200); // q1
        addPair(builder, "RUN:1:1101:5000:5101", 0, 3000, 3200); // q2
        if (withSecondReadGroup) {
            useReadGroup(builder, "rgB", "libA");
            addPair(builder, "RUN:2:1101:1000:1000", 0, 1000, 1200); // o5
        }
        return writeCoordinateSorted(builder);
    }

    static List<String> opticalDuplicates(final boolean withSecondReadGroup) {
        final List<String> duplicates = new ArrayList<>(List.of(
                "RUN:1:1101:1100:1000/1", "RUN:1:1101:1100:1000/2",
                "RUN:1:1101:1201:1000/1", "RUN:1:1101:1201:1000/2",
                "RUN:1:1102:1000:1000/1", "RUN:1:1102:1000:1000/2",
                "RUN:1:1101:5000:5101/1", "RUN:1:1101:5000:5101/2"));
        if (withSecondReadGroup) {
            duplicates.addAll(List.of("RUN:2:1101:1000:1000/1", "RUN:2:1101:1000:1000/2"));
        }
        return duplicates;
    }

    /**
     * Identical pairs and fragments in libraries libA (read groups rgA and rgC) and libB, and a read group of library
     * libD without reads. Records are only duplicates of records in the same library.
     */
    static SAMRecordSetBuilder multiLibraryInput() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "a1", 0, 1000, 1200);
        addPair(builder, "a2", 0, 1000, 1200);
        addFragment(builder, "a3", 0, 3000, false, "50M");
        useReadGroup(builder, "rgB", "libB");
        addPair(builder, "b1", 0, 1000, 1200);
        addFragment(builder, "b2", 0, 3000, false, "50M");
        addFragment(builder, "b3", 0, 3000, false, "50M");
        useReadGroup(builder, "rgC", "libA");
        addPair(builder, "c1", 0, 1000, 1200);
        useReadGroup(builder, "rgD", "libD");
        return builder;
    }

    /** The duplicates and metrics of {@link #multiLibraryInput()}; a library without reads still gets a row of zeros. */
    static void assertMultiLibraryResult(final List<SAMRecord> records, final Map<String, DuplicationMetrics> metrics) {
        assertDuplicates(records, List.of("a2/1", "a2/2", "c1/1", "c1/2", "b3/0"));
        Assert.assertEquals(metrics.size(), 3);
        assertDuplicationMetrics(metrics.get("libA"), "libA", 1, 3, 0, 0, 0, 2, 0);
        assertDuplicationMetrics(metrics.get("libB"), "libB", 2, 1, 0, 0, 1, 0, 0);
        assertDuplicationMetrics(metrics.get("libD"), "libD", 0, 0, 0, 0, 0, 0, 0);
    }

    /**
     * Pairs with the first end forward at chr1:1000 and the second end reverse at 500 on another contig: x1 and x2 on
     * chr2, x3 on chr3. Only x2 is a duplicate.
     */
    static SAMRecordSetBuilder interContigInput() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        builder.addPair("x1", 0, 1, 1000, 500, false, false, "50M", "50M", false, true, false, false, QUALITY);
        builder.addPair("x2", 0, 1, 1000, 500, false, false, "50M", "50M", false, true, false, false, QUALITY);
        builder.addPair("x3", 0, 2, 1000, 500, false, false, "50M", "50M", false, true, false, false, QUALITY);
        return builder;
    }

    /** Asserts that the records with the given keys, and no others, have the duplicate flag. */
    static void assertDuplicates(final Collection<SAMRecord> records, final Collection<String> expectedDuplicates) {
        final List<String> duplicates = new ArrayList<>();
        records.stream().filter(SAMRecord::getDuplicateReadFlag).forEach(r -> duplicates.add(key(r)));
        Assert.assertEqualsNoOrder(duplicates.toArray(), expectedDuplicates.toArray(), "Duplicates were " + duplicates);
    }

    /** Writes the records coordinate sorted, with SO:coordinate in the header. */
    static File writeCoordinateSorted(final SAMFileHeader header, final Collection<SAMRecord> records) throws IOException {
        final List<SAMRecord> sorted = new ArrayList<>(records);
        sorted.sort(new SAMRecordCoordinateComparator());
        return write(header, SAMFileHeader.SortOrder.coordinate, sorted);
    }

    static File writeCoordinateSorted(final SAMRecordSetBuilder builder) throws IOException {
        return writeCoordinateSorted(builder.getHeader(), builder.getRecords());
    }

    /** Writes the records queryname sorted, with SO:queryname in the header. */
    static File writeQueryNameSorted(final SAMFileHeader header, final Collection<SAMRecord> records) throws IOException {
        final List<SAMRecord> sorted = new ArrayList<>(records);
        sorted.sort(new SAMRecordQueryNameComparator());
        return write(header, SAMFileHeader.SortOrder.queryname, sorted);
    }

    /** Writes the records in the given order, setting the given sort order in a copy of the header. */
    static File write(final SAMFileHeader header, final SAMFileHeader.SortOrder sortOrder, final List<SAMRecord> records) throws IOException {
        final SAMFileHeader outputHeader = header.clone();
        outputHeader.setSortOrder(sortOrder);
        final File file = newTempFile("input", ".sam");
        try (final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(outputHeader, true, file)) {
            records.forEach(writer::addAlignment);
        }
        return file;
    }

    static File newTempFile(final String prefix, final String suffix) throws IOException {
        final File file = File.createTempFile(prefix, suffix);
        file.deleteOnExit();
        return file;
    }

    static SAMFileHeader readHeader(final File file) throws IOException {
        try (final SamReader reader = SamReaderFactory.makeDefault().open(file)) {
            return reader.getFileHeader();
        }
    }

    static List<SAMRecord> readRecords(final File file) throws IOException {
        final List<SAMRecord> records = new ArrayList<>();
        try (final SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(file)) {
            for (final SAMRecord record : reader) {
                records.add(record);
            }
        }
        return records;
    }

    /** Identifies a record by name, end of pair, and secondary/supplementary status. */
    static String key(final SAMRecord record) {
        final StringBuilder key = new StringBuilder(record.getReadName()).append('/');
        if (!record.getReadPairedFlag()) {
            key.append('0');
        } else {
            key.append(record.getFirstOfPairFlag() ? '1' : '2');
        }
        if (record.isSecondaryAlignment()) {
            key.append('s');
        }
        if (record.getSupplementaryAlignmentFlag()) {
            key.append('S');
        }
        return key.toString();
    }

    static List<String> keys(final Collection<SAMRecord> records) {
        final List<String> keys = new ArrayList<>(records.size());
        records.forEach(r -> keys.add(key(r)));
        return keys;
    }

    static Map<String, SAMRecord> byKey(final Collection<SAMRecord> records) {
        final Map<String, SAMRecord> map = new HashMap<>();
        records.forEach(r -> map.put(key(r), r));
        return map;
    }

    static <M extends MetricBase> List<M> readMetrics(final File file) throws IOException {
        final MetricsFile<M, Comparable<?>> metricsFile = new MetricsFile<>();
        try (final Reader reader = new FileReader(file)) {
            metricsFile.read(reader);
        }
        return metricsFile.getMetrics();
    }

    /** Reads the histograms of a metrics file, keyed by value label and then by bin, keeping only non-empty bins. */
    static Map<String, Map<Double, Double>> readHistograms(final File file) throws IOException {
        final MetricsFile<MetricBase, Double> metricsFile = new MetricsFile<>();
        try (final Reader reader = new FileReader(file)) {
            metricsFile.read(reader);
        }
        final Map<String, Map<Double, Double>> histograms = new HashMap<>();
        for (final Histogram<Double> histogram : metricsFile.getAllHistograms()) {
            final Map<Double, Double> bins = new HashMap<>();
            for (final Histogram.Bin<Double> bin : histogram.values()) {
                if (bin.getValue() != 0) {
                    bins.put(bin.getId(), bin.getValue());
                }
            }
            histograms.put(histogram.getValueLabel(), bins);
        }
        return histograms;
    }

    static Map<String, DuplicationMetrics> readDuplicationMetrics(final File file) throws IOException {
        final Map<String, DuplicationMetrics> byLibrary = new HashMap<>();
        for (final DuplicationMetrics metrics : MateCigarAndUmiTestSupport.<DuplicationMetrics>readMetrics(file)) {
            Assert.assertNull(byLibrary.put(metrics.LIBRARY, metrics), "Two rows for library " + metrics.LIBRARY);
        }
        return byLibrary;
    }

    /**
     * Asserts every count in a row of duplication metrics and the two derived values: PERCENT_DUPLICATION is the fraction of
     * examined mapped reads that are duplicates, and ESTIMATED_LIBRARY_SIZE is the Lander-Waterman estimate from the
     * non-optical pairs (see {@link #librarySize(long, long)}).
     */
    static void assertDuplicationMetrics(final DuplicationMetrics actual, final String library,
                                         final long unpairedReadsExamined, final long readPairsExamined,
                                         final long secondaryOrSupplementary, final long unmappedReads,
                                         final long unpairedReadDuplicates, final long readPairDuplicates,
                                         final long readPairOpticalDuplicates) {
        Assert.assertNotNull(actual, "No metrics for library " + library);
        Assert.assertEquals(actual.LIBRARY, library);
        Assert.assertEquals(actual.UNPAIRED_READS_EXAMINED, unpairedReadsExamined, "UNPAIRED_READS_EXAMINED");
        Assert.assertEquals(actual.READ_PAIRS_EXAMINED, readPairsExamined, "READ_PAIRS_EXAMINED");
        Assert.assertEquals(actual.SECONDARY_OR_SUPPLEMENTARY_RDS, secondaryOrSupplementary, "SECONDARY_OR_SUPPLEMENTARY_RDS");
        Assert.assertEquals(actual.UNMAPPED_READS, unmappedReads, "UNMAPPED_READS");
        Assert.assertEquals(actual.UNPAIRED_READ_DUPLICATES, unpairedReadDuplicates, "UNPAIRED_READ_DUPLICATES");
        Assert.assertEquals(actual.READ_PAIR_DUPLICATES, readPairDuplicates, "READ_PAIR_DUPLICATES");
        Assert.assertEquals(actual.READ_PAIR_OPTICAL_DUPLICATES, readPairOpticalDuplicates, "READ_PAIR_OPTICAL_DUPLICATES");

        final long examinedReads = unpairedReadsExamined + 2 * readPairsExamined;
        final double expectedPercent = examinedReads == 0 ? 0.0 : (unpairedReadDuplicates + 2.0 * readPairDuplicates) / examinedReads;
        Assert.assertEquals(actual.PERCENT_DUPLICATION, expectedPercent, 1e-6, "PERCENT_DUPLICATION");
        Assert.assertEquals(actual.ESTIMATED_LIBRARY_SIZE,
                librarySize(readPairsExamined - readPairOpticalDuplicates, readPairsExamined - readPairDuplicates),
                "ESTIMATED_LIBRARY_SIZE");
    }

    /**
     * Solves the Lander-Waterman equation C/X = 1 - exp(-N/X) for the library size X, given N pairs of which C are unique,
     * and truncates it to a whole number of molecules. There is no finite estimate without duplicates, so this returns null then.
     */
    static Long librarySize(final long pairs, final long uniquePairs) {
        if (pairs <= 0 || pairs - uniquePairs <= 0) {
            return null;
        }
        final double n = pairs;
        final double c = uniquePairs;
        // g(x) = c / x - 1 + exp(-n / x) is positive just above c and negative for large x.
        double low = c;
        double high = 2 * c;
        while (c / high - 1 + Math.exp(-n / high) > 0) {
            high *= 2;
        }
        for (int i = 0; i < 200; i++) {
            final double middle = (low + high) / 2;
            if (c / middle - 1 + Math.exp(-n / middle) > 0) {
                low = middle;
            } else {
                high = middle;
            }
        }
        return (long) low;
    }
}
