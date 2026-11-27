package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.utils.ValidationUtils;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;
import picard.nio.PicardHtsPath;
import picard.sam.util.SamTestUtil;
import picard.util.TestNGUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class PositionBasedDownsampleSamTest extends CommandLineProgramTest {
    final static String sample = "TestSample";
    final static String readGroupId = "TestReadGroup";
    final static String platform = "ILLUMINA";
    final static String library = "TestLibrary";

    final static File TEST_DIR = new File("testdata/picard/sam/PositionalDownsampleSam/");
    final File dict = new File(TEST_DIR, "header.dict");

    File tempDir;
    File tempSamFile;

    SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(readGroupId);

    @Override
    public String getCommandLineProgramName() {
        return PositionBasedDownsampleSam.class.getSimpleName();
    }

    //create a samfile with one tile and randomly placed reads within it.
    @BeforeTest
    void setupBuilder() throws IOException {
        final int numReads = 10000;
        final String flowCellBarcode = "TESTBARCODE";
        final int maxX = 10000;
        final int maxY = 20000;
        final int minX = 1000;
        final int minY = 2000;

        final String separator = ":";
        final int lane = 1;
        final int tile = 2203;

        //use fixed random seed to eliminate possibility of intermittent failures
        final Random rg = new Random(31);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);

        for (int i = 0; i < numReads; i++) {
            final int x = rg.nextInt(maxX) + minX;
            final int y = rg.nextInt(maxY) + minY;

            final String readName = flowCellBarcode + separator + lane + separator + tile + separator + x + separator + y;
            setBuilder.addPair(readName, 1, 1, 100);
        }

        tempDir = IOUtil.createTempDir("pds_test_PositionalDownsampling").toFile();
        tempSamFile = File.createTempFile("PositionalDownsampleSam", ".bam", tempDir);

        BufferedLineReader bufferedLineReader = null;
        try {
            bufferedLineReader = new BufferedLineReader(new FileInputStream(dict));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        final SAMFileHeader header = codec.decode(bufferedLineReader, dict.toString());

        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);

        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        header.addReadGroup(readGroupRecord);

        setBuilder.setHeader(header);

        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(header, false, tempSamFile);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }
        writer.close();
    }

    @AfterTest
    private void tearDown() {
        IOUtil.deleteDirectoryTree(tempDir);
    }

    @Test //this is actually a test of SAMRecordSetBuilder
    public void TestBuilder() {
        final ValidateSamFile validateSamFile = new ValidateSamFile();

        validateSamFile.INPUT = new PicardHtsPath(tempSamFile);
        Assert.assertEquals(validateSamFile.doWork(), 0);
    }

    @DataProvider(name = "ValidArgumentsTestProvider")
    public Object[][] ValidArgumentsTestProvider() {
        final List<Object[]> objects = new ArrayList<Object[]>();
        for (double i = 0.3; i <= 1; i += .1) {
            final Object[] array = {i};
            objects.add(array);
        }
        return objects.toArray(new Object[1][]);
    }

    // test removing some reads from a sparse, single tile bam
    @Test(dataProvider = "ValidArgumentsTestProvider")
    public void testDownsampleSingleTile(final double fraction) throws IOException {
        testDownsampleWorker(tempSamFile, fraction);
    }

    public void testDownsampleWorker(final File samFile, final double fraction) throws IOException {

        final File downsampled = File.createTempFile("PositionalDownsampleSam", ".bam", tempDir);
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + downsampled.getAbsolutePath(),
                "FRACTION=" + fraction,
                "CREATE_INDEX=true"
        };

        // make sure results is successful
        Assert.assertEquals(runPicardCommandLine(args), 0);

        // make sure that the resulting BAM is valid.
        final ValidateSamFile validateSamFile = new ValidateSamFile();

        validateSamFile.INPUT = new PicardHtsPath(downsampled);
        Assert.assertEquals(validateSamFile.doWork(), 0);

        //make sure that the total number of record in the resulting file in in the ballpark:
        TestNGUtil.assertGreaterThan(SamTestUtil.countSamTotalRecord(downsampled.toPath()), fraction * .8 * SamTestUtil.countSamTotalRecord(samFile.toPath()));
        TestNGUtil.assertLessThan(SamTestUtil.countSamTotalRecord(downsampled.toPath()), fraction * 1.2 * SamTestUtil.countSamTotalRecord(samFile.toPath()));
    }

    
    @DataProvider(name="allowTwiceData")
    public Object[][] allowTwiceData(){
        return new Object[][]{{true},{false}};
    }

    // test that downsampling twice yields an error.
    @Test(dataProvider = "allowTwiceData")
    public void TestInvalidTwice(final boolean allowMultiple) throws IOException {
        final File samFile = tempSamFile;
        final File downsampled = File.createTempFile("PositionalDownsampleSam", ".bam", tempDir);
        downsampled.deleteOnExit();
        final double fraction = .1;
        downsampled.deleteOnExit();
        final String[] args1 = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + downsampled.getAbsolutePath(),
                "FRACTION=" + fraction
        };
        runPicardCommandLine(args1);

        final File downsampledTwice = File.createTempFile("PositionalDownsampleSam", ".bam", tempDir);
        downsampledTwice.deleteOnExit();

        final List<String> args2 = new ArrayList<String>();

        args2.add("INPUT=" + downsampled.getAbsolutePath());
        args2.add("OUTPUT=" + downsampledTwice.getAbsolutePath());
        args2.add("FRACTION=" + fraction);

        if(allowMultiple) args2.add("ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS=true");

        //should blow up due to bad inputs
        if(allowMultiple)
            runPicardCommandLine(args2);
        else
            try {
                runPicardCommandLine(args2);
                throw new PicardException("Should have thrown an error!!");
            } catch (final PicardException ignored){

            }
    }

    // test that program fails on p<0  or p>1
    @DataProvider(name = "InvalidArgumentsTestProvider")
    public Object[][] InvalidArgumentsTestProvider() {
        return new Object[][]{{-1.0}, {-.00001}, {-5.0}, {1.00001}, {5.0}, {50.0}, {Double.MAX_VALUE}, {Double.POSITIVE_INFINITY}, {Double.NEGATIVE_INFINITY}};
    }

    @Test(dataProvider = "InvalidArgumentsTestProvider")
    public void TestInvalidArguments(final double fraction) throws IOException {
        final File samFile = tempSamFile;

        final File downsampled = File.createTempFile("PositionalDownsampleSam", ".bam", tempDir);
        downsampled.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "OUTPUT=" + downsampled.getAbsolutePath(),
                "FRACTION=" + fraction
        };
        //should blow up due to bad inputs
        Assert.assertNotEquals(runPicardCommandLine(args), 0);
    }

    /** Writes a coordinate-agnostic SAM/BAM whose read names are FC:lane:tile:x:y, optionally with a UMI suffix. */
    private File writeReadNameFixture(final File dir, final String name, final int numReads, final long seed,
                                      final boolean withUmi, final int lanes, final int tiles) throws IOException {
        final SAMFileHeader header = new SAMFileHeader();
        header.addSequence(new htsjdk.samtools.SAMSequenceRecord("chr1", 1000000));
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        final File file = new File(dir, name);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, file);
        final Random rg = new Random(seed);
        final String[] umis = {"ACGT", "TTGCA", "N", "AC+GT", "ACG-T"};
        for (int i = 0; i < numReads; i++) {
            final int lane = 1 + rg.nextInt(lanes);
            final int tile = 1101 + rg.nextInt(tiles);
            final int x = 1000 + rg.nextInt(30000);
            final int y = 1000 + rg.nextInt(30000);
            String readName = "FC123:" + lane + ":" + tile + ":" + x + ":" + y;
            if (withUmi) {
                readName += ":" + umis[rg.nextInt(umis.length)];
            }
            final SAMRecord rec = new SAMRecord(header);
            rec.setReadName(readName);
            rec.setReferenceIndex(0);
            rec.setAlignmentStart(1 + i % 900000);
            rec.setCigarString("10M");
            rec.setReadString("ACGTACGTAC");
            rec.setBaseQualityString("IIIIIIIIII");
            rec.setReadUnmappedFlag(false);
            writer.addAlignment(rec);
        }
        writer.close();
        return file;
    }

    /** Sorted list of "name/flag" for every record in the file: a strict, order-independent view of the kept reads. */
    private static List<String> readKeptRecords(final File file) throws IOException {
        final List<String> kept = new ArrayList<>();
        try (final htsjdk.samtools.SamReader reader = htsjdk.samtools.SamReaderFactory.makeDefault().open(file)) {
            for (final SAMRecord rec : reader) {
                kept.add(rec.getReadName() + "/" + rec.getFlags() + "/" + rec.getAlignmentStart());
            }
        }
        java.util.Collections.sort(kept);
        return kept;
    }

    private static final String ILLUMINA_UMI_EQUIVALENT_REGEX = "(?:.*:)?([0-9]+):([0-9]+):([0-9]+)(:[ACGTN+-]+)?$";

    @Test
    public void testUmiReadNameFormatMatchesEquivalentRegex() throws IOException {
        final File dir = IOUtil.createTempDir("pds_umi").toFile();
        try {
            final File input = writeReadNameFixture(dir, "umi_in.bam", 20000, 11L, true, 2, 3);
            final File viaRegex = new File(dir, "umi_regex.bam");
            final File viaFormat = new File(dir, "umi_format.bam");

            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "INPUT=" + input, "OUTPUT=" + viaRegex, "FRACTION=0.3",
                    "READ_NAME_REGEX=" + ILLUMINA_UMI_EQUIVALENT_REGEX}), 0);
            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "INPUT=" + input, "OUTPUT=" + viaFormat, "FRACTION=0.3",
                    "READ_NAME_FORMAT=ILLUMINA_UMI"}), 0);

            final List<String> expected = readKeptRecords(viaRegex);
            Assert.assertFalse(expected.isEmpty());
            Assert.assertEquals(readKeptRecords(viaFormat), expected);
        } finally {
            IOUtil.deleteDirectoryTree(dir);
        }
    }

    @Test
    public void testReadNameFormatAndRegexAreExclusive() {
        final String[] args = new String[]{
                "INPUT=" + tempSamFile.getAbsolutePath(),
                "OUTPUT=" + tempSamFile.getAbsolutePath() + ".unused.bam",
                "FRACTION=0.5",
                "READ_NAME_FORMAT=ILLUMINA_UMI",
                "READ_NAME_REGEX=" + ILLUMINA_UMI_EQUIVALENT_REGEX};
        Assert.assertNotEquals(runPicardCommandLine(args), 0);
    }

    /** Single-pass is approximate in general. On this input (few reads per tile, stable extremes) it keeps the same reads. */
    @Test
    public void testSinglePassMatchesTwoPassOnSyntheticInput() throws IOException {
        final File dir = IOUtil.createTempDir("pds_single").toFile();
        try {
            final File input = writeReadNameFixture(dir, "in.bam", 20000, 21L, true, 2, 3);
            final File twoPass = new File(dir, "two.bam");
            final File singlePass = new File(dir, "single.bam");

            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "INPUT=" + input, "OUTPUT=" + twoPass, "FRACTION=0.4",
                    "READ_NAME_FORMAT=ILLUMINA_UMI", "TWO_PASS=true"}), 0);
            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "INPUT=" + input, "OUTPUT=" + singlePass, "FRACTION=0.4",
                    "READ_NAME_FORMAT=ILLUMINA_UMI", "READ_COUNT=20000"}), 0);

            Assert.assertFalse(readKeptRecords(twoPass).isEmpty());
            Assert.assertEquals(readKeptRecords(singlePass), readKeptRecords(twoPass));
        } finally {
            IOUtil.deleteDirectoryTree(dir);
        }
    }

    /** Without READ_COUNT, the default (single-pass) run must fall back to the exact two-pass result. */
    @Test
    public void testDefaultWithoutReadCountIsTwoPass() throws IOException {
        final File dir = IOUtil.createTempDir("pds_default").toFile();
        try {
            final File input = writeReadNameFixture(dir, "in.bam", 5000, 23L, false, 2, 2);
            final File defaults = new File(dir, "default.bam");
            final File twoPass = new File(dir, "two.bam");

            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "INPUT=" + input, "OUTPUT=" + defaults, "FRACTION=0.4"}), 0);
            Assert.assertEquals(runPicardCommandLine(new String[]{
                    "INPUT=" + input, "OUTPUT=" + twoPass, "FRACTION=0.4", "TWO_PASS=true"}), 0);

            Assert.assertEquals(readKeptRecords(defaults), readKeptRecords(twoPass));
        } finally {
            IOUtil.deleteDirectoryTree(dir);
        }
    }

    /** A READ_COUNT far from the real read count scales the widening wrongly, so single-pass must refuse to finish. */
    @Test
    public void testSinglePassRejectsWrongReadCount() throws IOException {
        final File dir = IOUtil.createTempDir("pds_single_count").toFile();
        try {
            final File input = writeReadNameFixture(dir, "in.bam", 20000, 22L, false, 1, 2);
            final File out = new File(dir, "out.bam");
            try {
                final int rc = runPicardCommandLine(new String[]{
                        "INPUT=" + input, "OUTPUT=" + out, "FRACTION=0.4",
                        "READ_COUNT=10000"});
                Assert.assertNotEquals(rc, 0, "a READ_COUNT 50% off must fail");
            } catch (final PicardException expected) {
                // acceptable: the failure surfaces as an exception
            }
        } finally {
            IOUtil.deleteDirectoryTree(dir);
        }
    }

}