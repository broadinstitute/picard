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

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;

import static picard.sam.markduplicates.MateCigarAndUmiTestSupport.*;

/**
 * Tests of {@link MarkDuplicatesWithMateCigarIterator} used directly: the records it returns and their order, its
 * counters and metrics after it is closed, the getters that require it to be closed, and the inputs it rejects.
 */
public class MarkDuplicatesWithMateCigarIteratorTest {

    private static CloseableIterator<SAMRecord> iteratorOver(final List<SAMRecord> records) {
        final Iterator<SAMRecord> iterator = records.iterator();
        return new CloseableIterator<SAMRecord>() {
            @Override
            public void close() {
            }

            @Override
            public boolean hasNext() {
                return iterator.hasNext();
            }

            @Override
            public SAMRecord next() {
                return iterator.next();
            }
        };
    }

    private static List<SAMRecord> coordinateSorted(final SAMRecordSetBuilder builder) {
        final List<SAMRecord> records = new ArrayList<>(builder.getRecords());
        records.sort(new SAMRecordCoordinateComparator());
        return records;
    }

    private static MarkDuplicatesWithMateCigarIterator newIterator(final SAMFileHeader header, final List<SAMRecord> records,
                                                                   final ScoringStrategy scoringStrategy,
                                                                   final boolean removeDuplicates,
                                                                   final boolean skipPairsWithNoMateCigar) {
        final File tmpDir = IOUtil.createTempDir("MarkDuplicatesWithMateCigarIteratorTest").toFile();
        tmpDir.deleteOnExit();
        return new MarkDuplicatesWithMateCigarIterator(header, iteratorOver(records), new OpticalDuplicateFinder(),
                scoringStrategy, 200, removeDuplicates, skipPairsWithNoMateCigar, 4, 2, List.of(tmpDir));
    }

    private static SAMFileHeader coordinateHeader(final SAMRecordSetBuilder builder) {
        builder.getHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return builder.getHeader();
    }

    private static List<SAMRecord> drain(final MarkDuplicatesWithMateCigarIterator iterator) {
        final List<SAMRecord> records = new ArrayList<>();
        while (iterator.hasNext()) {
            records.add(iterator.next());
        }
        return records;
    }

    @Test
    public void testReturnsAllRecordsInOrderAndCountsAfterClose() {
        final SAMRecordSetBuilder builder = mixedInput();
        final List<SAMRecord> input = coordinateSorted(builder);
        final List<String> inputOrder = keys(input);

        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(coordinateHeader(builder), input,
                ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true);
        final List<SAMRecord> output = drain(iterator);
        Assert.assertFalse(iterator.hasNext());
        Assert.assertThrows(NoSuchElementException.class, iterator::next);

        Assert.assertEquals(keys(output), inputOrder);
        assertDuplicates(output, MIXED_DUPLICATES);

        iterator.close();
        Assert.assertEquals(iterator.getNumDuplicates(), MIXED_DUPLICATES.size());
        Assert.assertEquals(iterator.getNumRecordsWithNoMateCigar(), 0);
        Assert.assertEquals(iterator.getOpticalDupesByLibraryId().getSumOfValues(), 0.0);

        // Pair counts are per read until the program halves them.
        final DuplicationMetrics metrics = iterator.getLibraryIdGenerator().getMetricsByLibrary(MIXED_LIBRARY);
        Assert.assertEquals(metrics.UNPAIRED_READS_EXAMINED, 7);
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, 14);
        Assert.assertEquals(metrics.SECONDARY_OR_SUPPLEMENTARY_RDS, 2);
        Assert.assertEquals(metrics.UNMAPPED_READS, 3);
        Assert.assertEquals(metrics.UNPAIRED_READ_DUPLICATES, 4);
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES, 6);
    }

    @Test
    public void testOpticalDuplicatesByLibrary() throws IOException {
        final File input = opticalInput(true);
        final SAMFileHeader header = readHeader(input);
        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(header, readRecords(input),
                ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true);
        assertDuplicates(drain(iterator), opticalDuplicates(true));
        iterator.close();

        final short libraryId = iterator.getLibraryIdGenerator().getLibraryIdsMap().get("libA");
        Assert.assertEquals(iterator.getOpticalDupesByLibraryId().get(libraryId).getValue(), 1.0);
        Assert.assertEquals(iterator.getOpticalDupesByLibraryId().getSumOfValues(), 1.0);
        Assert.assertEquals(iterator.getNumDuplicates(), 10);
    }

    @Test
    public void testRemoveDuplicates() {
        final SAMRecordSetBuilder builder = mixedInput();
        final List<SAMRecord> input = coordinateSorted(builder);
        final List<String> expected = keys(input).stream().filter(k -> !MIXED_DUPLICATES.contains(k)).collect(Collectors.toList());

        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(coordinateHeader(builder), input,
                ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, true, true);
        final List<SAMRecord> output = drain(iterator);
        iterator.close();

        Assert.assertEquals(keys(output), expected);
        Assert.assertEquals(iterator.getNumDuplicates(), MIXED_DUPLICATES.size());
    }

    @Test
    public void testEmptyInput() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(coordinateHeader(builder), List.of(),
                ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true);
        Assert.assertFalse(iterator.hasNext());
        Assert.assertThrows(NoSuchElementException.class, iterator::next);
        iterator.close();
        Assert.assertEquals(iterator.getNumDuplicates(), 0);
        Assert.assertNotNull(iterator.getLibraryIdGenerator().getMetricsByLibrary("libA"));
    }

    @Test
    public void testCountsRecordsWithoutMateCigar() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        addPair(builder, "p2", 0, 1000, 1200).forEach(r -> r.setAttribute(SAMTag.MC.name(), null));
        addPair(builder, "p3", 0, 1000, 1200);

        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(coordinateHeader(builder), coordinateSorted(builder),
                ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true);
        final List<SAMRecord> output = drain(iterator);
        iterator.close();

        assertDuplicates(output, List.of("p3/1", "p3/2"));
        Assert.assertEquals(iterator.getNumRecordsWithNoMateCigar(), 2);
        Assert.assertEquals(iterator.getNumDuplicates(), 2);
    }

    @Test
    public void testIteratorContract() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(coordinateHeader(builder), coordinateSorted(builder),
                ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true);

        Assert.assertSame(iterator.assertSorted(SAMFileHeader.SortOrder.coordinate), iterator);
        Assert.assertThrows(IllegalStateException.class, () -> iterator.assertSorted(SAMFileHeader.SortOrder.queryname));
        Assert.assertThrows(IllegalStateException.class, () -> iterator.assertSorted(SAMFileHeader.SortOrder.unsorted));
        Assert.assertThrows(UnsupportedOperationException.class, iterator::remove);

        // Counters are only available once the iterator is closed.
        Assert.assertThrows(PicardException.class, iterator::getNumDuplicates);
        Assert.assertThrows(PicardException.class, iterator::getNumRecordsWithNoMateCigar);
        Assert.assertThrows(PicardException.class, iterator::getLibraryIdGenerator);
        Assert.assertThrows(PicardException.class, iterator::getOpticalDupesByLibraryId);

        Assert.assertEquals(drain(iterator).size(), 2);
        iterator.close();
        Assert.assertEquals(iterator.getNumDuplicates(), 0);
    }

    @Test
    public void testRejectsUnsortedHeader() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        final SAMFileHeader header = builder.getHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        final PicardException exception = Assert.expectThrows(PicardException.class, () -> newIterator(header,
                coordinateSorted(builder), ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true));
        Assert.assertTrue(exception.getMessage().contains("coordinate sort order"), exception.getMessage());
    }

    @Test
    public void testRejectsSumOfBaseQualities() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200);
        final PicardException exception = Assert.expectThrows(PicardException.class, () -> newIterator(coordinateHeader(builder),
                coordinateSorted(builder), ScoringStrategy.SUM_OF_BASE_QUALITIES, false, true));
        Assert.assertTrue(exception.getMessage().contains("SUM_OF_BASE_QUALITIES"), exception.getMessage());
    }

    @Test
    public void testRejectsMissingMateCigarWhenNotSkipping() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        addPair(builder, "p1", 0, 1000, 1200).forEach(r -> r.setAttribute(SAMTag.MC.name(), null));
        final PicardException exception = Assert.expectThrows(PicardException.class, () -> newIterator(coordinateHeader(builder),
                coordinateSorted(builder), ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, false));
        Assert.assertTrue(exception.getMessage().contains("p1"), exception.getMessage());
    }

    @Test
    public void testRejectsDecreasingReferenceIndex() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final SAMRecord onSecondContig = addFragment(builder, "f2", 1, 1000, false, "50M");
        final SAMRecord onFirstContig = addFragment(builder, "f1", 0, 1000, false, "50M");
        final PicardException exception = Assert.expectThrows(PicardException.class, () -> newIterator(coordinateHeader(builder),
                List.of(onSecondContig, onFirstContig), ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true));
        Assert.assertTrue(exception.getMessage().contains("out of order"), exception.getMessage());
    }

    @Test
    public void testRejectsDecreasingPosition() {
        final SAMRecordSetBuilder builder = newBuilder();
        useReadGroup(builder, "rgA", "libA");
        final SAMRecord mapped = addFragment(builder, "f1", 0, 1000, false, "50M");
        // Unmapped reads placed on a contig skip the duplicate queue, so only the order check on output can catch this.
        final SAMRecord placedUnmapped = builder.addFrag("u1", 0, 500, false, true, null, null, QUALITY);

        final MarkDuplicatesWithMateCigarIterator iterator = newIterator(coordinateHeader(builder),
                List.of(mapped, placedUnmapped), ScoringStrategy.TOTAL_MAPPED_REFERENCE_LENGTH, false, true);
        Assert.assertTrue(iterator.hasNext());
        final PicardException exception = Assert.expectThrows(PicardException.class, iterator::next);
        Assert.assertTrue(exception.getMessage().contains("coordinate sort order"), exception.getMessage());
    }
}
