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

package picard.sam.markduplicates.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.sam.DuplicationMetrics;
import picard.sam.markduplicates.MarkDuplicates;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

/**
 * Tests the methods of {@link AbstractMarkDuplicatesCommandLineProgram} and
 * {@link AbstractOpticalDuplicateFinderCommandLineProgram} that MarkDuplicates and its relatives share:
 * trackOpticalDuplicates and addSingletonToCount (optical duplicate counts and the set-size histograms),
 * addReadToLibraryMetrics, getChainedPgIds, openInputs, and setupOpticalDuplicateFinder.
 */
public class AbstractMarkDuplicatesCommandLineProgramMethodsTest {

    private static SAMFileHeader header(final SAMFileHeader.SortOrder sortOrder, final String... readGroupsAndLibraries) {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(sortOrder);
        header.setSequenceDictionary(new SAMSequenceDictionary(Arrays.asList(
                new SAMSequenceRecord("chr1", 10000), new SAMSequenceRecord("chr2", 10000))));
        for (int i = 0; i < readGroupsAndLibraries.length; i += 2) {
            final SAMReadGroupRecord readGroup = new SAMReadGroupRecord(readGroupsAndLibraries[i]);
            readGroup.setSample("sample");
            readGroup.setLibrary(readGroupsAndLibraries[i + 1]);
            header.addReadGroup(readGroup);
        }
        return header;
    }

    private static ReadEndsForMarkDuplicates end(final short libraryId, final byte opticalOrientation, final int x, final int y) {
        final ReadEndsForMarkDuplicates end = new ReadEndsForMarkDuplicates();
        end.libraryId = libraryId;
        end.orientation = ReadEnds.FR;
        end.orientationForOpticalDuplicates = opticalOrientation;
        end.readGroup = 0;
        end.setTile((short) 1101);
        end.setX(x);
        end.setY(y);
        return end;
    }

    private static Map<Double, Double> bins(final Histogram<Double> histogram) {
        final Map<Double, Double> bins = new TreeMap<>();
        for (final Histogram.Bin<Double> bin : histogram.values()) {
            bins.put(bin.getId(), bin.getValue());
        }
        return bins;
    }

    private static Map<Double, Double> expectedBins(final double... sizeAndCount) {
        final Map<Double, Double> bins = new TreeMap<>();
        for (int i = 0; i < sizeAndCount.length; i += 2) {
            bins.put(sizeAndCount[i], sizeAndCount[i + 1]);
        }
        return bins;
    }

    private static double opticalCount(final LibraryIdGenerator generator, final short libraryId) {
        final Histogram.Bin<Short> bin = generator.getOpticalDuplicatesByLibraryIdMap().get(libraryId);
        return bin == null ? 0 : bin.getValue();
    }

    @Test
    public void testTrackOpticalDuplicates() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(SAMFileHeader.SortOrder.coordinate));
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, null);

        // keeper, one end within 100 pixels of it and one far from both
        final List<ReadEndsForMarkDuplicates> set = Arrays.asList(
                end((short) 1, ReadEnds.FR, 1000, 1000), end((short) 1, ReadEnds.FR, 1100, 900), end((short) 1, ReadEnds.FR, 5000, 5000));
        AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(set, set.get(0), finder, generator);
        Assert.assertFalse(set.get(0).isOpticalDuplicate);
        Assert.assertTrue(set.get(1).isOpticalDuplicate);
        Assert.assertFalse(set.get(2).isOpticalDuplicate);
        Assert.assertEquals(opticalCount(generator, (short) 1), 1.0);
        Assert.assertEquals(generator.getNumberOfOpticalDuplicateClusters(), 1);
        Assert.assertEquals(bins(generator.getDuplicateCountHist()), expectedBins(3, 1));
        Assert.assertEquals(bins(generator.getOpticalDuplicateCountHist()), expectedBins(2, 1));
        Assert.assertEquals(bins(generator.getNonOpticalDuplicateCountHist()), expectedBins(2, 1));

        // a second set, in another library, with no optical duplicates
        final List<ReadEndsForMarkDuplicates> second = Arrays.asList(
                end((short) 2, ReadEnds.FR, 1000, 1000), end((short) 2, ReadEnds.FR, 1101, 1000));
        AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(second, second.get(0), finder, generator);
        Assert.assertFalse(second.get(1).isOpticalDuplicate);
        Assert.assertEquals(opticalCount(generator, (short) 2), 0.0);
        Assert.assertEquals(bins(generator.getDuplicateCountHist()), expectedBins(2, 1, 3, 1));
        Assert.assertEquals(bins(generator.getOpticalDuplicateCountHist()), expectedBins(2, 1));
        Assert.assertEquals(bins(generator.getNonOpticalDuplicateCountHist()), expectedBins(2, 2));

        // a third set in the first library: two of four ends are optical duplicates
        final List<ReadEndsForMarkDuplicates> third = Arrays.asList(
                end((short) 1, ReadEnds.FR, 1000, 1000), end((short) 1, ReadEnds.FR, 1050, 1000),
                end((short) 1, ReadEnds.FR, 1150, 1000), end((short) 1, ReadEnds.FR, 3000, 1000));
        AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(third, third.get(0), finder, generator);
        Assert.assertEquals(opticalCount(generator, (short) 1), 3.0);
        Assert.assertEquals(bins(generator.getDuplicateCountHist()), expectedBins(2, 1, 3, 1, 4, 1));
        Assert.assertEquals(bins(generator.getOpticalDuplicateCountHist()), expectedBins(2, 1, 3, 1));
        Assert.assertEquals(bins(generator.getNonOpticalDuplicateCountHist()), expectedBins(2, 3));
    }

    /** Ends whose first of pair is on different strands are only compared with ends like them. */
    @Test
    public void testTrackOpticalDuplicatesSeparatesFirstOfPairStrands() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(SAMFileHeader.SortOrder.coordinate));
        final OpticalDuplicateFinder finder = new OpticalDuplicateFinder(OpticalDuplicateFinder.DEFAULT_READ_NAME_REGEX, 100, null);
        final List<ReadEndsForMarkDuplicates> set = Arrays.asList(
                end((short) 1, ReadEnds.FR, 1000, 1000), end((short) 1, ReadEnds.RF, 1010, 1000),
                end((short) 1, ReadEnds.FR, 1020, 1000), end((short) 1, ReadEnds.RF, 1030, 1000));
        AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates(set, set.get(0), finder, generator);
        Assert.assertFalse(set.get(0).isOpticalDuplicate);
        Assert.assertTrue(set.get(2).isOpticalDuplicate);
        // the two RF ends are optical duplicates of each other, so exactly one of them is flagged
        Assert.assertTrue(set.get(1).isOpticalDuplicate ^ set.get(3).isOpticalDuplicate);
        Assert.assertEquals(opticalCount(generator, (short) 1), 2.0);
        Assert.assertEquals(bins(generator.getDuplicateCountHist()), expectedBins(4, 1));
        Assert.assertEquals(bins(generator.getOpticalDuplicateCountHist()), expectedBins(3, 1));
        Assert.assertEquals(bins(generator.getNonOpticalDuplicateCountHist()), expectedBins(2, 1));
    }

    @Test
    public void testAddSingletonToCount() {
        final LibraryIdGenerator generator = new LibraryIdGenerator(header(SAMFileHeader.SortOrder.coordinate));
        AbstractMarkDuplicatesCommandLineProgram.addSingletonToCount(generator);
        AbstractMarkDuplicatesCommandLineProgram.addSingletonToCount(generator);
        Assert.assertEquals(bins(generator.getDuplicateCountHist()), expectedBins(1, 2));
        Assert.assertEquals(bins(generator.getNonOpticalDuplicateCountHist()), expectedBins(1, 2));
        Assert.assertEquals(bins(generator.getOpticalDuplicateCountHist()), expectedBins());
    }

    @Test
    public void testAddReadToLibraryMetrics() {
        final SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, "rg1", "libA", "rg2", "libB");
        final LibraryIdGenerator generator = new LibraryIdGenerator(header);

        final SAMRecord fragment = new SAMRecord(header);
        fragment.setAttribute(SAMTag.RG.name(), "rg1");
        fragment.setReferenceIndex(0);
        fragment.setAlignmentStart(10);
        final SAMRecord pairEnd = fragment.deepCopy();
        pairEnd.setReadPairedFlag(true);
        final SAMRecord mateUnmapped = pairEnd.deepCopy();
        mateUnmapped.setMateUnmappedFlag(true);
        final SAMRecord secondary = pairEnd.deepCopy();
        secondary.setSecondaryAlignment(true);
        final SAMRecord supplementary = fragment.deepCopy();
        supplementary.setSupplementaryAlignmentFlag(true);
        final SAMRecord unmapped = pairEnd.deepCopy();
        unmapped.setReadUnmappedFlag(true);
        final SAMRecord unmappedSecondary = unmapped.deepCopy();
        unmappedSecondary.setSecondaryAlignment(true);

        for (final SAMRecord rec : Arrays.asList(fragment, pairEnd, pairEnd, mateUnmapped, secondary, supplementary, unmapped, unmappedSecondary)) {
            AbstractMarkDuplicatesCommandLineProgram.addReadToLibraryMetrics(rec, header, generator, false);
        }
        final DuplicationMetrics libA = generator.getMetricsByLibrary("libA");
        Assert.assertEquals(libA.UNPAIRED_READS_EXAMINED, 2);
        Assert.assertEquals(libA.SECONDARY_OR_SUPPLEMENTARY_RDS, 2);
        Assert.assertEquals(libA.UNMAPPED_READS, 2);
        Assert.assertEquals(generator.getMetricsByLibrary("libB").UNPAIRED_READS_EXAMINED, 0);

        // a read without a read group gets metrics for an unknown library
        final SAMRecord noReadGroup = fragment.deepCopy();
        noReadGroup.setAttribute(SAMTag.RG.name(), null);
        AbstractMarkDuplicatesCommandLineProgram.addReadToLibraryMetrics(noReadGroup, header, generator, false);
        final DuplicationMetrics unknown = generator.getMetricsByLibrary("Unknown Library");
        Assert.assertEquals(unknown.LIBRARY, "Unknown Library");
        Assert.assertEquals(unknown.UNPAIRED_READS_EXAMINED, 1);
        Assert.assertEquals(generator.getMetricsByLibraryMap().size(), 3);
    }

    @Test
    public void testGetChainedPgIds() {
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.PROGRAM_GROUP_VERSION = "1.0";
        markDuplicates.PROGRAM_GROUP_COMMAND_LINE = "command";
        markDuplicates.pgIdsSeen.addAll(Arrays.asList("bwa", null));
        final SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate);
        header.addProgramRecord(new SAMProgramRecord("bwa"));

        final Map<String, String> chained = markDuplicates.getChainedPgIds(header);
        Assert.assertEquals(chained.keySet(), new HashSet<>(Arrays.asList("bwa", null)));
        Assert.assertEquals(new HashSet<>(chained.values()).size(), 2);
        Assert.assertTrue(chained.containsValue("MarkDuplicates"));
        Assert.assertEquals(header.getProgramRecords().size(), 3);
        for (final Map.Entry<String, String> entry : chained.entrySet()) {
            final SAMProgramRecord pg = header.getProgramRecord(entry.getValue());
            Assert.assertTrue(entry.getValue().startsWith("MarkDuplicates"));
            Assert.assertEquals(pg.getPreviousProgramGroupId(), entry.getKey());
            Assert.assertEquals(pg.getProgramName(), "MarkDuplicates");
            Assert.assertEquals(pg.getProgramVersion(), "1.0");
            Assert.assertEquals(pg.getCommandLine(), "command");
        }
    }

    @Test
    public void testGetChainedPgIdsWithoutProgramRecordId() {
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.PROGRAM_RECORD_ID = null;
        markDuplicates.pgIdsSeen.add("bwa");
        final SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate);
        Assert.assertNull(markDuplicates.getChainedPgIds(header));
        Assert.assertTrue(header.getProgramRecords().isEmpty());
    }

    private static File write(final File file, final SAMFileHeader header, final String readGroup, final int... starts) {
        try (final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(header, true, file)) {
            for (final int start : starts) {
                final SAMRecord rec = new SAMRecord(header);
                rec.setReadName(readGroup + "-" + start);
                rec.setReferenceIndex(0);
                rec.setAlignmentStart(start);
                rec.setCigarString("4M");
                rec.setReadString("ACGT");
                rec.setBaseQualityString("IIII");
                rec.setAttribute(SAMTag.RG.name(), readGroup);
                writer.addAlignment(rec);
            }
        }
        return file;
    }

    private static List<String> drain(final AbstractMarkDuplicatesCommandLineProgram.SamHeaderAndIterator headerAndIterator) {
        final List<String> names = new ArrayList<>();
        headerAndIterator.iterator.forEachRemaining(r -> names.add(r.getReadName()));
        headerAndIterator.iterator.close();
        return names;
    }

    @Test
    public void testOpenInputsAssumeSorted() {
        final File dir = IOUtil.createTempDir("OpenInputs").toFile();
        try {
            final File in = write(new File(dir, "in.sam"), header(SAMFileHeader.SortOrder.unsorted, "rg1", "libA"), "rg1", 5, 9);
            final MarkDuplicates asIsReader = new MarkDuplicates();
            asIsReader.INPUT = Collections.singletonList(in.getAbsolutePath());
            final AbstractMarkDuplicatesCommandLineProgram.SamHeaderAndIterator asIs = asIsReader.openInputs(false);
            Assert.assertEquals(asIs.header.getSortOrder(), SAMFileHeader.SortOrder.unsorted);
            Assert.assertEquals(drain(asIs), Arrays.asList("rg1-5", "rg1-9"));

            final MarkDuplicates assumeSorted = new MarkDuplicates();
            assumeSorted.INPUT = Collections.singletonList(in.getAbsolutePath());
            assumeSorted.ASSUME_SORTED = true;
            final AbstractMarkDuplicatesCommandLineProgram.SamHeaderAndIterator opened = assumeSorted.openInputs(true);
            Assert.assertEquals(opened.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
            Assert.assertEquals(drain(opened), Arrays.asList("rg1-5", "rg1-9"));

            final MarkDuplicates assumeQueryname = new MarkDuplicates();
            assumeQueryname.INPUT = Collections.singletonList(in.getAbsolutePath());
            assumeQueryname.ASSUME_SORT_ORDER = SAMFileHeader.SortOrder.queryname;
            final AbstractMarkDuplicatesCommandLineProgram.SamHeaderAndIterator queryname = assumeQueryname.openInputs(false);
            Assert.assertEquals(queryname.header.getSortOrder(), SAMFileHeader.SortOrder.queryname);
            drain(queryname);
        } finally {
            IOUtil.recursiveDelete(dir.toPath());
        }
    }

    @Test
    public void testOpenInputsMergesSeveralInputs() {
        final File dir = IOUtil.createTempDir("OpenInputs").toFile();
        try {
            final File first = write(new File(dir, "first.sam"), header(SAMFileHeader.SortOrder.coordinate, "rg1", "libA"), "rg1", 5, 20);
            final File second = write(new File(dir, "second.sam"), header(SAMFileHeader.SortOrder.coordinate, "rg2", "libB"), "rg2", 1, 10, 30);
            final MarkDuplicates markDuplicates = new MarkDuplicates();
            markDuplicates.INPUT = Arrays.asList(first.getAbsolutePath(), second.getAbsolutePath());
            final AbstractMarkDuplicatesCommandLineProgram.SamHeaderAndIterator opened = markDuplicates.openInputs(false);
            Assert.assertEquals(opened.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
            Assert.assertNotNull(opened.header.getReadGroup("rg1"));
            Assert.assertNotNull(opened.header.getReadGroup("rg2"));
            Assert.assertEquals(drain(opened), Arrays.asList("rg2-1", "rg1-5", "rg2-10", "rg1-20", "rg2-30"));
        } finally {
            IOUtil.recursiveDelete(dir.toPath());
        }
    }

    @Test
    public void testSetupOpticalDuplicateFinderUsesTheArguments() {
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.READ_NAME_REGEX = "([0-9]+)_([0-9]+)_([0-9]+)";
        markDuplicates.OPTICAL_DUPLICATE_PIXEL_DISTANCE = 5;
        markDuplicates.MAX_OPTICAL_DUPLICATE_SET_SIZE = 2;
        markDuplicates.setupOpticalDuplicateFinder();

        final OpticalDuplicateFinder finder = markDuplicates.opticalDuplicateFinder;
        Assert.assertNotNull(finder);
        Assert.assertEquals(finder.opticalDuplicatePixelDistance, 5);
        final ReadEndsForMarkDuplicates location = new ReadEndsForMarkDuplicates();
        Assert.assertTrue(finder.addLocationInformation("12_345_678", location));
        Assert.assertEquals(location.getTile(), 12);
        Assert.assertEquals(location.getX(), 345);
        Assert.assertEquals(location.getY(), 678);
        Assert.assertFalse(finder.addLocationInformation("RUN:1:12:345:678", location));

        // two ends 5 pixels apart are optical duplicates, three are more than MAX_OPTICAL_DUPLICATE_SET_SIZE
        final List<ReadEndsForMarkDuplicates> two = Arrays.asList(end((short) 1, ReadEnds.FR, 100, 100), end((short) 1, ReadEnds.FR, 105, 95));
        Assert.assertEquals(finder.findOpticalDuplicates(two, two.get(0)), new boolean[]{false, true});
        final List<ReadEndsForMarkDuplicates> three = Arrays.asList(end((short) 1, ReadEnds.FR, 100, 100),
                end((short) 1, ReadEnds.FR, 105, 95), end((short) 1, ReadEnds.FR, 101, 101));
        Assert.assertEquals(finder.findOpticalDuplicates(three, three.get(0)), new boolean[]{false, false, false});
        final List<ReadEndsForMarkDuplicates> tooFar = Arrays.asList(end((short) 1, ReadEnds.FR, 100, 100), end((short) 1, ReadEnds.FR, 106, 100));
        Assert.assertEquals(finder.findOpticalDuplicates(tooFar, tooFar.get(0)), new boolean[]{false, false});
    }
}
