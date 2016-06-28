package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.testng.Assert.assertEquals;

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

        tempDir = IOUtil.createTempDir("pds_test", "PositionalDownsampling");
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

        validateSamFile.INPUT = tempSamFile;
        assertEquals(validateSamFile.doWork(), 0);
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
        assertEquals(runPicardCommandLine(args), 0);

        // make sure that the resulting BAM is valid.
        final ValidateSamFile validateSamFile = new ValidateSamFile();

        validateSamFile.INPUT = downsampled;
        assertEquals(validateSamFile.doWork(), 0);

        //make sure that the total number of record in the resulting file in in the ballpark:
        assertGreaterThan(countSamTotalRecord(downsampled), fraction * .8 * countSamTotalRecord(samFile));
        assertLessThan(countSamTotalRecord(downsampled), fraction * 1.2 * countSamTotalRecord(samFile));
    }

    private long countSamTotalRecord(final File samFile) {
        final SamReader reader = SamReaderFactory.make().open(samFile);
        assert reader.hasIndex();
        long total = 0;

        for (int i = 0; i < reader.getFileHeader().getSequenceDictionary().size(); i++) {
            total += reader.indexing().getIndex().getMetaData(i).getAlignedRecordCount();
            total += reader.indexing().getIndex().getMetaData(i).getUnalignedRecordCount();
        }
        return total;
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
        assert runPicardCommandLine(args) != 0;
    }

    // these fit in a TestNG.Utils class, but there isn't one in picard-public...
    static public void assertGreaterThan(final double lhs, final double rhs) {
        Assert.assertTrue(lhs > rhs, String.format("Expected inequality is not true: %g > %g", lhs, rhs));
    }

    static public void assertLessThan(final double lhs, final double rhs) {
        Assert.assertTrue(lhs < rhs, String.format("Expected inequality is not true: %g < %g", lhs, rhs));
    }

}