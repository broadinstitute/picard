package picard.analysis;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.SortSam;

import java.io.*;
import java.util.Random;

/**
 * Tests for methods in CollectWgsMetrics
 *
 *
 * @author Kylee Bergin
 */

public class CollectWgsMetricsTest extends CommandLineProgramTest {

    private final static File REF_DICT_DIR = new File("testdata/picard/sam/CollectGcBiasMetrics/");
    private final static File TEST_DIR = new File("testdata/picard/sam/");
    private final File referenceDict = new File(REF_DICT_DIR, "MSmallHeader.dict");
    private File tempSamFile;
    private File outfile;

    private final static int LENGTH = 99;
    private final static String SAMPLE = "TestSample1";
    private final static String READ_GROUP_ID = "TestReadGroup1";
    private final static String PLATFORM = "ILLUMINA";
    private final static String LIBRARY = "TestLibrary1";
    private final static int NUM_READS = 40000;

    public String getCommandLineProgramName() {
        return CollectWgsMetrics.class.getSimpleName();
    }

    @DataProvider(name = "wgsDataProvider")
    public Object[][] wgsDataProvider() {
        final String referenceFile = "testdata/picard/quality/chrM.reference.fasta";

        return new Object[][] {
                {tempSamFile, outfile, referenceFile}
        };
    }

    @Test(dataProvider = "wgsDataProvider")
    public void testMetricsFromWGS(final File input, final File outfile, final String referenceFile) throws IOException {
        outfile.deleteOnExit();
        final int sampleSize = 1000;

        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "SAMPLE_SIZE=" + sampleSize
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetrics.WgsMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_COVERAGE, 13.985155, .02);
            Assert.assertEquals(metrics.PCT_EXC_OVERLAP, 0.0);  // 52 of 606 bases
            Assert.assertEquals(metrics.PCT_EXC_BASEQ, 0.399906, .02);    // 114 of 606 bases
            Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.0);    // 202 of 606 bases
            Assert.assertEquals(metrics.SD_COVERAGE, 57.364434, .02);
            Assert.assertEquals(metrics.MEDIAN_COVERAGE, 0.0);
            Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.0);
            Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.0);
            Assert.assertEquals(metrics.PCT_EXC_CAPPED, 0.519542, .001);
            Assert.assertEquals(metrics.PCT_EXC_TOTAL, 0.919537, .001);
            Assert.assertEquals(metrics.PCT_1X, 0.056364, .0001);
            Assert.assertEquals(metrics.PCT_5X, 0.056364, .0001);
            Assert.assertEquals(metrics.PCT_10X, 0.056364, .0001);
            Assert.assertEquals(metrics.PCT_15X, 0.056364, .0001);
            Assert.assertEquals(metrics.PCT_20X, 0.056364, .0001);
            Assert.assertEquals(metrics.PCT_25X, 0.056303, .0001);
            Assert.assertEquals(metrics.PCT_30X, 0.056303, .0001);
            Assert.assertEquals(metrics.PCT_40X, 0.056243, .0001);
            Assert.assertEquals(metrics.PCT_50X, 0.056243, .0001);
            Assert.assertEquals(metrics.PCT_60X, 0.056182, .0001);
            Assert.assertEquals(metrics.PCT_70X, 0.056182, .0001);
            Assert.assertEquals(metrics.PCT_80X, 0.056122, .0001);
            Assert.assertEquals(metrics.PCT_90X, 0.056062, .0001);
            Assert.assertEquals(metrics.PCT_100X, 0.056062, .0001);
            Assert.assertEquals(metrics.HET_SNP_SENSITIVITY, 0.056362, .02);
            Assert.assertEquals(metrics.HET_SNP_Q, 0.0);

        }
    }

    //create a samfile for testing.
    @BeforeTest
    void setupBuilder() throws IOException {
        final String readName = "TESTBARCODE";

        //Create Sam Files
        tempSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        final File tempSamFileUnsorted = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        tempSamFileUnsorted.deleteOnExit();
        tempSamFile.deleteOnExit();
        final SAMFileHeader header = new SAMFileHeader();

        //Check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(referenceDict));
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        } catch (final SAMException e) {
            e.printStackTrace();
        }

        //Set readGroupRecord
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(READ_GROUP_ID);
        readGroupRecord.setSample(SAMPLE);
        readGroupRecord.setPlatform(PLATFORM);
        readGroupRecord.setLibrary(LIBRARY);
        readGroupRecord.setPlatformUnit(READ_GROUP_ID);
        header.addReadGroup(readGroupRecord);

        //Add to setBuilder
        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);
        setBuilder.setHeader(header);

        //Read settings
        final String separator = ":";
        final int ID = 1;
        final int maxReadStart = 800;
        final int minReadStart = 1;
        final Random rg = new Random(5);

        for (int i = 0; i < NUM_READS; i++) {
            final int start = rg.nextInt(maxReadStart) + minReadStart;
            final String newReadName = readName + separator + ID + separator + i;
            setBuilder.addPair(newReadName, 0, start + ID, start + ID + LENGTH);
        }

        //Write SAM file
        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(header, false, tempSamFileUnsorted);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }
        writer.close();

        //sort the temp file
        final SortSam sorter = new SortSam();
        final String[] args = new String[]{
                "INPUT=" + tempSamFileUnsorted.getAbsolutePath(),
                "OUTPUT=" + tempSamFile.getAbsolutePath(),
                "SORT_ORDER=coordinate"
        };

        sorter.instanceMain(args);

        //create output files for tests
        outfile = File.createTempFile("testWgsMetrics", ".txt");
        outfile.deleteOnExit();
    }
}