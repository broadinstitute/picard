package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.SortSam;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Random;
import static org.testng.Assert.assertEquals;

/**
 * Created by kbergin on 3/26/15 to test GcBias MultiLevel Collector.
 */
public class CollectGcBiasMetricsTest extends CommandLineProgramTest {
    private final static String sample1 = "TestSample1";
    private final static String sample2 = "TestSample2";
    private final static String readGroupId1 = "TestReadGroup1";
    private final static String readGroupId2 = "TestReadGroup2";
    private final static String readGroupId3 = "TestReadGroup3";
    private final static String platform = "ILLUMINA";
    private final static String library1 = "TestLibrary1";
    private final static String library2 = "TestLibrary2";
    private final static String library3 = "TestLibrary3";

    private final static File TEST_DIR = new File("testdata/picard/sam/CollectGcBiasMetrics/");
    private final File dict = new File(TEST_DIR, "Mheader.dict");

    File tempSamFile;

    SAMRecordSetBuilder setBuilder1 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMRecordSetBuilder setBuilder2 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMRecordSetBuilder setBuilder3 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMReadGroupRecord readGroupRecord1 = new SAMReadGroupRecord(readGroupId1);
    SAMReadGroupRecord readGroupRecord2 = new SAMReadGroupRecord(readGroupId2);
    SAMReadGroupRecord readGroupRecord3 = new SAMReadGroupRecord(readGroupId3);

    //create a samfile with different samples, read groups and libraries that overlap for testing.
    @BeforeTest
    void setupBuilder() throws IOException {
        final int numReads = 100;
        final String flowCellBarcode = "TESTBARCODE";
        final String readName = flowCellBarcode;

        tempSamFile = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        File tempSamFileUnsorted = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        tempSamFileUnsorted.deleteOnExit();
        tempSamFile.deleteOnExit();

        BufferedLineReader bufferedLineReader = null;
        try {
            bufferedLineReader = new BufferedLineReader(new FileInputStream(dict));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
        final SAMFileHeader header = codec.decode(bufferedLineReader, dict.toString());
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        //build different levels to put into the same bam file for testing multi level collection
        setup(numReads, readName, 1, readGroupId1, readGroupRecord1, sample1, library1, header, setBuilder1); //Sample 1, Library 1, RG 1
        setup(numReads, readName, 2, readGroupId2, readGroupRecord2, sample1, library2, header, setBuilder2); //Sample 1, Library 2, RG 2
        setup(numReads, readName, 3, readGroupId3, readGroupRecord3, sample2, library3, header, setBuilder3); //Sample 2, Library 3, RG 3

        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(header, false, tempSamFileUnsorted);

        for (final SAMRecord record : setBuilder1) {
            writer.addAlignment(record);
        }
        for (final SAMRecord record : setBuilder2) {
            writer.addAlignment(record);
        }
        for (final SAMRecord record : setBuilder3) {
            writer.addAlignment(record);
        }
        writer.close();

        //sort the temp file
        final SortSam sorter = new SortSam();
        final String[] args = new String[]{"INPUT=" + tempSamFileUnsorted.getAbsolutePath(), "OUTPUT=" + tempSamFile.getAbsolutePath(), "SORT_ORDER=coordinate"};

        sorter.instanceMain(args);
    }

    void setup(final int numReads,
               final String readName,
               final int ID,
               final String readGroupId,
               final SAMReadGroupRecord readGroupRecord,
               final String sample,
               final String library,
               final SAMFileHeader header,
               final SAMRecordSetBuilder setBuilder) throws IOException {

        final String separator = ":";
        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);
        readGroupRecord.setPlatformUnit(readGroupId);
        header.addReadGroup(readGroupRecord);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);

        setBuilder.setHeader(header);

        final int max = 15000;
        final int min = 1;
        final Random rg = new Random(5);

        for (int i = 0; i < numReads; i++) {
            final int start = rg.nextInt(max) + min;
            final String newReadName = readName + separator + ID + separator + i;
            setBuilder.addPair(newReadName, 0, start+ID, start+ID+99);
        }
    }

    public String getCommandLineProgramName() {
        return CollectGcBiasMetrics.class.getSimpleName();
    }

    @Test //test all collection levels
    public void test() throws IOException{
        runTest(tempSamFile);
    }

    public void runTest(final File input) throws IOException {
        final File outfile = File.createTempFile("test", ".gc_bias_summary_metrics");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias_detail_metrics");
        final File pdf = File.createTempFile("test", ".pdf");
        final String referenceFile = "testdata/picard/quality/chrM.reference.fasta";
        final int windowSize = 100;
        final double minGenFraction = 1.0E-5;
        final boolean biSulfiteSeq = false;
        final boolean assumeSorted = false;
        outfile.deleteOnExit();
        detailsOutfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + detailsOutfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "SUMMARY_OUTPUT=" + outfile.getAbsolutePath(),
                "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                "WINDOW_SIZE=" + windowSize,
                "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                "LEVEL=ALL_READS",
                "LEVEL=SAMPLE",
                "LEVEL=READ_GROUP",
                "ASSUME_SORTED=" + assumeSorted
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<GcBiasSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals("All Reads")) { //ALL_READS level
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 300);
                Assert.assertEquals(metrics.ALIGNED_READS, 600);
                Assert.assertEquals(metrics.AT_DROPOUT, 7.234062);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.086217);
            } else if (metrics.READ_GROUP != null && metrics.READ_GROUP.equals("TestReadGroup1")) { //Library 1
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 9.20674);
                Assert.assertEquals(metrics.GC_DROPOUT, 3.834244);
            } else if (metrics.READ_GROUP != null && metrics.READ_GROUP.equals("TestReadGroup2")) {//Library 2
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 10.144505);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.08986);
            } else if (metrics.READ_GROUP != null && metrics.READ_GROUP.equals("TestReadGroup3")) {//Library 3
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 9.229205);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.977838);
            } else if (metrics.SAMPLE != null && metrics.SAMPLE.equals("TestSample1")) {//Library 1 and 2
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 200);
                Assert.assertEquals(metrics.ALIGNED_READS, 400);
                Assert.assertEquals(metrics.AT_DROPOUT, 7.410747);
                Assert.assertEquals(metrics.GC_DROPOUT, 3.83986);
            } else if (metrics.SAMPLE != null && metrics.SAMPLE.equals("TestSample2")) {//Library 3
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 9.229205);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.977838);
            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }
}
