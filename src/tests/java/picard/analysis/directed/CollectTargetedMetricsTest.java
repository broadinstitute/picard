package picard.analysis.directed;

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

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;

public class CollectTargetedMetricsTest extends CommandLineProgramTest {
    private final static File TEST_DIR = new File("testdata/picard/sam/CollectGcBiasMetrics/");
    private final File dict = new File(TEST_DIR, "Mheader.dict");
    private File tempSamFile;
    private File outfile;
    private File perTargetOutfile;
    private final static int LENGTH = 99;

    private final static String sample = "TestSample1";
    private final static String readGroupId = "TestReadGroup1";
    private final static String platform = "ILLUMINA";
    private final static String library = "TestLibrary1";
    private final static int numReads = 40000;

    @Override
    public String getCommandLineProgramName() {
        return CollectTargetedPcrMetrics.class.getSimpleName();
    }

    //create a samfile with 40000 reads for testing whether a cap is found.
    @BeforeTest
    void setupBuilder() throws IOException {
        final String readName = "TESTBARCODE";

        //Create Sam Files
        tempSamFile = File.createTempFile("CollectTargetedMetrics", ".bam", TEST_DIR);
        final File tempSamFileUnsorted = File.createTempFile("CollectTargetedMetrics", ".bam", TEST_DIR);
        tempSamFileUnsorted.deleteOnExit();
        tempSamFile.deleteOnExit();
        final SAMFileHeader header = new SAMFileHeader();

        //Check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(dict));
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        } catch (final SAMException e) {
            e.printStackTrace();
        }

        //Set readGroupRecord
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(readGroupId);
        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);
        readGroupRecord.setPlatformUnit(readGroupId);
        header.addReadGroup(readGroupRecord);

        //Add to setBuilder
        final SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);
        setBuilder.setHeader(header);

        //Read settings
        final String separator = ":";
        final int ID = 1;
        final int max = 15000;
        final int min = 1;
        final Random rg = new Random(5);

        for (int i = 0; i < numReads; i++) {
            final int start = rg.nextInt(max) + min;
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
        outfile = File.createTempFile("test", ".TargetedMetrics_Coverage");
        perTargetOutfile = File.createTempFile("perTarget", ".perTargetCoverage");
        outfile.deleteOnExit();
        perTargetOutfile.deleteOnExit();
    }

    @DataProvider(name = "targetedIntervalDataProvider")
    public Object[][] targetedIntervalDataProvider() {
        final String referenceFile = "testdata/picard/quality/chrM.reference.fasta";
        final String emptyIntervals = "testdata/picard/quality/chrM.empty.interval_list";
        final String singleIntervals = "testdata/picard/quality/chrM.single.interval_list";

        return new Object[][] {
                {tempSamFile, outfile, perTargetOutfile, referenceFile, singleIntervals, 1000},
                {tempSamFile, outfile, perTargetOutfile, referenceFile, emptyIntervals, 1000}
        };
    }

    @Test(dataProvider = "targetedIntervalDataProvider")
    public void runCollectTargetedMetricsTest(final File input, final File outfile, final File perTargetOutfile, final String referenceFile,
                                final String targetIntervals, final int sampleSize) throws IOException {

        final String[] args = new String[] {
                "TARGET_INTERVALS=" + targetIntervals,
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "PER_TARGET_COVERAGE=" + perTargetOutfile.getAbsolutePath(),
                "LEVEL=ALL_READS",
                "AMPLICON_INTERVALS=" + targetIntervals,
                "SAMPLE_SIZE=" + sampleSize
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<TargetedPcrMetrics, Comparable<?>> output = new MetricsFile<TargetedPcrMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final TargetedPcrMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_READS, numReads * 2);
            Assert.assertEquals(metrics.HET_SNP_SENSITIVITY, .997972, .02);
        }
    }
}