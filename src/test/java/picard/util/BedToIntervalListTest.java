package picard.util;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.IOUtil;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

/**
 * @author nhomer
 */
public class BedToIntervalListTest {

    private static final String TEST_DATA_DIR = "testdata/picard/util/BedToIntervalListTest";

    private void doTest(final String inputBed, final String header) throws IOException, SAMException {
        doTest(inputBed, header, true);
    }

    private void doTest(final String inputBed, final String header, final boolean keepLengthZero) throws IOException, SAMException {
        final File outputFile  = File.createTempFile("bed_to_interval_list_test.", ".interval_list");
        outputFile.deleteOnExit();
        final BedToIntervalList program = new BedToIntervalList();
        final File inputBedFile = new File(TEST_DATA_DIR, inputBed);
        program.INPUT = inputBedFile;
        program.SEQUENCE_DICTIONARY = new File(TEST_DATA_DIR, header);
        program.OUTPUT = outputFile;
        program.UNIQUE = true;
        program.KEEP_LENGTH_ZERO_INTERVALS = keepLengthZero;
        program.doWork();

        // Assert they are equal
        IOUtil.assertFilesEqual(new File(inputBedFile.getAbsolutePath() + ".interval_list"), outputFile);
    }


    @Test(dataProvider = "testBedToIntervalListDataProvider")
    public void testBedToIntervalList(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam");
    }

    @Test(dataProvider = "testLengthZeroIntervalsSkippedProvider")
    public void testLengthZeroIntervalsSkipped(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam", false);
    }

    // test a fixed bed file using different dictionaries
    @Test(dataProvider = "testBedToIntervalListSequenceDictionaryDataProvider")
    public void testBedToIntervalListSequenceDictionary(final String dictionary) throws IOException {
        doTest("seq_dict_test.bed", dictionary);
    }

    // test for bad dictionaries - we expect these to throw exceptions
    @Test(dataProvider = "testBedToIntervalListSequenceDictionaryBadDataProvider",
          expectedExceptions = {SAMException.class, PicardException.class})
    public void testBedToIntervalListBadSequenceDictionary(final String dictionary) throws IOException {
        doTest("seq_dict_test.bed", dictionary);
    }

    @Test(dataProvider = "testBedToIntervalListOutOfBoundsDataProvider", expectedExceptions = PicardException.class)
    public void testBedToIntervalListOutOfBounds(final String inputBed) throws IOException {
        doTest(inputBed, "header.sam");
    }

    @Test(expectedExceptions = PicardException.class)
    public void testRejectIntervalListInput() throws IOException {
        // Feeding an interval_list file (which starts with @ SAM headers) to BedToIntervalList
        // should throw a clear PicardException rather than a confusing NumberFormatException.
        doTest("seq_dict_test.dictionary.interval_list", "header.sam");
    }

    @Test
    public void testStdinSupport() throws IOException {
        // Feed the contents of simple.bed through System.in and verify the output matches
        // the expected interval_list — proving that /dev/stdin is now fully supported.
        final File simpleBed = new File(TEST_DATA_DIR, "simple.bed");
        final byte[] bedBytes = Files.readAllBytes(simpleBed.toPath());

        final InputStream originalIn = System.in;
        try {
            System.setIn(new ByteArrayInputStream(bedBytes));

            final File outputFile = File.createTempFile("bed_to_interval_list_stdin_test.", ".interval_list");
            outputFile.deleteOnExit();

            final BedToIntervalList program = new BedToIntervalList();
            program.INPUT = new File("/dev/stdin");
            program.SEQUENCE_DICTIONARY = new File(TEST_DATA_DIR, "header.sam");
            program.OUTPUT = outputFile;
            program.UNIQUE = true;

            program.doWork();

            IOUtil.assertFilesEqual(new File(simpleBed.getAbsolutePath() + ".interval_list"), outputFile);
        } finally {
            System.setIn(originalIn);
        }
    }

    @DataProvider
    public Object[][] testBedToIntervalListDataProvider() {
        return new Object[][]{
                {"simple.bed"},
                {"overlapping.bed"},
                {"extended.bed"},
                {"one_base_interval.bed"},
                {"zero_base_interval.bed"},
                {"first_base_in_contig.bed"},
                {"zero_length_interval_at_first_position_in_contig.bed"},
                {"last_base_in_contig.bed"},
                {"multi_contig.bed"}
        };
    }

    // test for each of the file categories supported by SAMSequenceDictionaryExtractor
    @DataProvider
    public Object[][] testBedToIntervalListSequenceDictionaryDataProvider() {
        return new Object[][]{
                {"seq_dict_test.dictionary.interval_list"},
                {"seq_dict_test.dictionary.fasta"},
                {"seq_dict_test.dictionary.dict"},
                {"seq_dict_test.dictionary.sam"},
                {"seq_dict_test.dictionary.vcf"}
        };
    }

    @DataProvider
    public Object[][] testBedToIntervalListSequenceDictionaryBadDataProvider() {
        return new Object[][]{
                // a file that does not represent a dictionary
                {"seq_dict_test.dictionary.bad"},
                // a file that is a valid dictionary, but missing contigs in the dictionary
                {"seq_dict_test.dictionary.bad.vcf"}
        };
    }

    @DataProvider
    public Object[][] testBedToIntervalListOutOfBoundsDataProvider() {
        return new Object[][]{
                {"end_after_chr.bed"},
                {"end_before_chr.bed"},
                {"missing_chr.bed"},
                {"start_after_chr.bed"},
                {"start_before_chr.bed"},
                {"off_by_one_interval.bed"}
        };
    }

    @DataProvider
    public Object[][] testLengthZeroIntervalsSkippedProvider() {
        return new Object[][]{
                {"zero_length_test.bed"}
        };
    }
}
