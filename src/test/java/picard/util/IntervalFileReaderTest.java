package picard.util;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.List;

/**
 * Tests for {@link IntervalFileReader}: format detection and generic interval loading.
 */
public class IntervalFileReaderTest {

    private static final String TEST_DATA_DIR = "testdata/picard/util/BedToIntervalListTest";

    /** Build a minimal dictionary with chr1..chr8 each of length 1 000 000. */
    private static SAMSequenceDictionary buildDictionary() {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (int i = 1; i <= 8; i++) {
            dict.addSequence(new SAMSequenceRecord("chr" + i, 1_000_000));
        }
        return dict;
    }

    private static SAMFileHeader buildHeader() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(buildDictionary());
        return header;
    }

    private static BufferedReader readerOf(final String content) {
        return new BufferedReader(new StringReader(content));
    }

    // -------------------------------------------------------------------------
    // Format detection tests
    // -------------------------------------------------------------------------

    @Test
    public void testDetectBedFormat() throws IOException {
        try (final BufferedReader reader = readerOf("chr1\t100\t200\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.BED);
        }
    }

    @Test
    public void testDetectIntervalListFormat() throws IOException {
        try (final BufferedReader reader = readerOf("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000000\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.INTERVAL_LIST);
        }
    }

    @Test
    public void testDetectUnknownFormat() throws IOException {
        try (final BufferedReader reader = readerOf("not_valid_data\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.UNKNOWN);
            Assert.assertEquals(result.firstLine(), "not_valid_data");
        }
    }

    @Test
    public void testDetectSkipsComments() throws IOException {
        try (final BufferedReader reader = readerOf("# comment\n\nchr1\t100\t200\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.BED);
        }
    }

    @Test
    public void testDetectEmptyFile() throws IOException {
        try (final BufferedReader reader = readerOf("")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.UNKNOWN);
            Assert.assertNull(result.firstLine());
        }
    }

    @Test
    public void testDetectOnlyComments() throws IOException {
        try (final BufferedReader reader = readerOf("# just a comment\n# another\n")) {
            final IntervalFileReader.FormatDetectionResult result = IntervalFileReader.detectIntervalFormat(reader);
            Assert.assertEquals(result.format(), IntervalFileReader.IntervalFileFormat.UNKNOWN);
            Assert.assertNull(result.firstLine());
        }
    }

    // -------------------------------------------------------------------------
    // fromBed(File, ...) tests
    // -------------------------------------------------------------------------

    @Test
    public void testFromBedFileSimple() {
        final File bedFile = new File(TEST_DATA_DIR, "simple.bed");
        final IntervalList result = IntervalFileReader.fromBed(bedFile, buildHeader(), false, false, null);
        final List<Interval> intervals = result.getIntervals();
        Assert.assertEquals(intervals.size(), 2);
        // simple.bed: chr1 100 2000 and chr1 3000 4000
        // BEDCodec: 0-based half-open -> 1-based closed: (100,2000) -> (101,2000), (3000,4000) -> (3001,4000)
        Assert.assertEquals(intervals.get(0).getContig(), "chr1");
        Assert.assertEquals(intervals.get(0).getStart(), 101);
        Assert.assertEquals(intervals.get(0).getEnd(), 2000);
        Assert.assertEquals(intervals.get(1).getContig(), "chr1");
        Assert.assertEquals(intervals.get(1).getStart(), 3001);
        Assert.assertEquals(intervals.get(1).getEnd(), 4000);
    }

    @Test
    public void testFromBedFileExtendedFields() {
        final File bedFile = new File(TEST_DATA_DIR, "extended.bed");
        final IntervalList result = IntervalFileReader.fromBed(bedFile, buildHeader(), false, false, null);
        Assert.assertFalse(result.getIntervals().isEmpty());
    }

    // -------------------------------------------------------------------------
    // loadIntervals(File, ...) tests
    // -------------------------------------------------------------------------

    @Test
    public void testLoadIntervalsFromBedFile() {
        final File bedFile = new File(TEST_DATA_DIR, "simple.bed");
        final IntervalList result = IntervalFileReader.loadIntervals(bedFile, buildDictionary());
        Assert.assertEquals(result.getIntervals().size(), 2);
        // simple.bed: chr1 100 2000  and  chr1 3000 4000
        Assert.assertEquals(result.getIntervals().get(0).getStart(), 101);
        Assert.assertEquals(result.getIntervals().get(0).getEnd(), 2000);
        Assert.assertEquals(result.getIntervals().get(1).getStart(), 3001);
        Assert.assertEquals(result.getIntervals().get(1).getEnd(), 4000);
    }

    @Test
    public void testLoadIntervalsFromIntervalListFile() {
        final File ilFile = new File(TEST_DATA_DIR, "seq_dict_test.dictionary.interval_list");
        // interval_list format doesn't actually need the dictionary arg (it has its own header),
        // but we pass one anyway since loadIntervals requires it for the BED code path.
        final IntervalList result = IntervalFileReader.loadIntervals(ilFile, buildDictionary());
        Assert.assertFalse(result.getIntervals().isEmpty());
    }
}
