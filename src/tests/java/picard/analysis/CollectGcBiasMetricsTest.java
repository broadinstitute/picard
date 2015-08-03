package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.SAMException;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.SortSam;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

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
    private final static int LENGTH = 99;
    private final static int NUM_READS = 100;
    private final static String READ_NAME = "TESTBARCODE";

    private final static File TEST_DIR = new File("testdata/picard/sam/CollectGcBiasMetrics/");
    private final File dict = new File(TEST_DIR, "MNOheader.dict");

    File tempSamFileChrM_O;
    File tempSamFileAllChr;

    SAMRecordSetBuilder setBuilder1 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMRecordSetBuilder setBuilder2 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMRecordSetBuilder setBuilder3 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    SAMRecordSetBuilder setBuilder4 = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);

    SAMReadGroupRecord readGroupRecord1 = new SAMReadGroupRecord(readGroupId1);
    SAMReadGroupRecord readGroupRecord2 = new SAMReadGroupRecord(readGroupId2);
    SAMReadGroupRecord readGroupRecord3 = new SAMReadGroupRecord(readGroupId3);

    /////////////////////////////////////////////////////////////////////////////
    //create two Sam Files.
    //One with different samples, read groups and libraries that overlap for runGcBiasMultiLevelTest. Reads will align to chrM and O, not N.
    //Second Sam file is one sample/read group/library but has reads that align to all three chr (M,N,O). For runWindowsComparisonTest.
    /////////////////////////////////////////////////////////////////////////////

    @BeforeTest
    void setupBuilder() throws IOException {
        tempSamFileChrM_O = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        tempSamFileAllChr = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        tempSamFileChrM_O.deleteOnExit();
        tempSamFileAllChr.deleteOnExit();

        final File tempSamFileUnsorted = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        tempSamFileUnsorted.deleteOnExit();


        final SAMFileHeader header = new SAMFileHeader();

        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(dict));
            header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        } catch (final SAMException e) {
            e.printStackTrace();
        }

        //build different levels to put into the same bam file for testing multi level collection
        setupTest1(1, readGroupId1, readGroupRecord1, sample1, library1, header, setBuilder1); //Sample 1, Library 1, RG 1
        setupTest1(2, readGroupId2, readGroupRecord2, sample1, library2, header, setBuilder2); //Sample 1, Library 2, RG 2
        setupTest1(3, readGroupId3, readGroupRecord3, sample2, library3, header, setBuilder3); //Sample 2, Library 3, RG 3

        //build one last readgroup for comparing that window count stays the same whether you use all contigs or not
        setupTest2(1, readGroupId1, readGroupRecord1, sample1, library1, header, setBuilder4);

        final List<SAMRecordSetBuilder> test1Builders = new ArrayList<SAMRecordSetBuilder>();
        test1Builders.add(setBuilder1);
        test1Builders.add(setBuilder2);
        test1Builders.add(setBuilder3);

        final List<SAMRecordSetBuilder> test2Builders = new ArrayList<SAMRecordSetBuilder>();
        test2Builders.add(setBuilder4);

        tempSamFileChrM_O = build(test1Builders, tempSamFileUnsorted, header);
        tempSamFileAllChr = build(test2Builders, tempSamFileUnsorted, header);
    }

    public String getCommandLineProgramName() {
        return CollectGcBiasMetrics.class.getSimpleName();
    }

    /////////////////////////////////////////////////////////////////////////////
    //This test checks the functionality of the gc bias code. Compares values from running a generated temporary Sam file through
    // CollectGcBiasMetrics to manually-calculated values.
    /////////////////////////////////////////////////////////////////////////////
    @Test
    public void runGcBiasMultiLevelTest() throws IOException {
        final File outfile = File.createTempFile("test", ".gc_bias.summary_metrics");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias.detail_metrics");
        outfile.deleteOnExit();
        detailsOutfile.deleteOnExit();

        runGcBias(tempSamFileChrM_O, outfile, detailsOutfile);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<GcBiasSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals("All Reads")) { //ALL_READS level
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 300);
                Assert.assertEquals(metrics.ALIGNED_READS, 600);
                Assert.assertEquals(metrics.AT_DROPOUT, 21.624498);
                Assert.assertEquals(metrics.GC_DROPOUT, 3.525922);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.831374);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.049672);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else if (metrics.READ_GROUP != null && metrics.READ_GROUP.equals("TestReadGroup1")) { //Library 1
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 23.627784);
                Assert.assertEquals(metrics.GC_DROPOUT, 2.582877);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.793584);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.060382);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else if (metrics.READ_GROUP != null && metrics.READ_GROUP.equals("TestReadGroup2")) {//Library 2
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 23.784958);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.025922);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.816258);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.053956);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else if (metrics.READ_GROUP != null && metrics.READ_GROUP.equals("TestReadGroup3")) {//Library 3
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 21.962578);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.559328);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.88428);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.034676);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else if (metrics.SAMPLE != null && metrics.SAMPLE.equals("TestSample1")) {//Library 1 and 2
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 200);
                Assert.assertEquals(metrics.ALIGNED_READS, 400);
                Assert.assertEquals(metrics.AT_DROPOUT, 23.194597);
                Assert.assertEquals(metrics.GC_DROPOUT, 3.275922);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.804921);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.057169);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else if (metrics.SAMPLE != null && metrics.SAMPLE.equals("TestSample2")) {//Library 3
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 100);
                Assert.assertEquals(metrics.ALIGNED_READS, 200);
                Assert.assertEquals(metrics.AT_DROPOUT, 21.962578);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.559328);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.88428);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.034676);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    //Compare GcBiasDetailMetrics output file from test1 which only has reads that align to chrM and chrO, but not chrN (in the middle) to
    // a GcBiasDetailMetrics output file that has reads aligned to all three chromosomes in this reference file. The number of 100bp windows
    // found across the whole reference should be the same regardless of where records align.
    //This test ensures that there is not a bug in calculating the gc windows.
    /////////////////////////////////////////////////////////////////////////////
    @Test
    public void runWindowsComparisonTest() throws IOException {
        final File outfile = File.createTempFile("test", ".gc_bias_summary_metrics");
        final File allChrOutFile = File.createTempFile("testAllChr", ".gc_bias_summary_metrics");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias_detail_metrics");
        final File allChrDetailsOutfile = File.createTempFile("testAllChrDetails", ".gc_bias_detail_metrics");
        outfile.deleteOnExit();
        allChrOutFile.deleteOnExit();
        detailsOutfile.deleteOnExit();
        allChrDetailsOutfile.deleteOnExit();

        runGcBias(tempSamFileChrM_O, outfile, detailsOutfile);
        runGcBias(tempSamFileAllChr, allChrOutFile, allChrDetailsOutfile);

        final MetricsFile<GcBiasDetailMetrics, Comparable<?>> outputDetails = new MetricsFile<GcBiasDetailMetrics, Comparable<?>>();
        outputDetails.read(new FileReader(detailsOutfile));
        final List<GcBiasDetailMetrics> details = outputDetails.getMetrics();

        final MetricsFile<GcBiasDetailMetrics, Comparable<?>> outputAllChrDetails = new MetricsFile<GcBiasDetailMetrics, Comparable<?>>();
        outputAllChrDetails.read(new FileReader(allChrDetailsOutfile));

        int i = 0;

        //Output for the two sam files are only the same for the "All Reads" level
        for (final GcBiasDetailMetrics metrics : outputAllChrDetails.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals("All Reads")) {
                Assert.assertEquals(metrics.WINDOWS, details.get(i).WINDOWS);
                i++;
            }
            else {break;}
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // Writes the setBuilders to a SAMFileWriter and sorts the sam.
    // Takes in a list of SAMRecordSetBuilders because of the multi-level collection: setBuilders cannot take in more than one read group
    // or library or sample, so there are separate ones for each type when testing multi-level collection.
    /////////////////////////////////////////////////////////////////////////////
    public File build (final List<SAMRecordSetBuilder> setBuilder, final File unsortedSam, final SAMFileHeader header) throws IOException {
        final File sortedSam = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        sortedSam.deleteOnExit();

        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(header, false, unsortedSam);

        for( final SAMRecordSetBuilder subSetBuilder : setBuilder){
            for (final SAMRecord record : subSetBuilder) {
                writer.addAlignment(record);
            }
        }
        writer.close();

        final SortSam sorter = new SortSam();
        final String[] args = new String[] {
                "INPUT=" + unsortedSam.getAbsolutePath(),
                "OUTPUT=" + sortedSam.getAbsolutePath(),
                "SORT_ORDER=coordinate"
        };

        sorter.instanceMain(args);

        return sortedSam;
    }

    /////////////////////////////////////////////////////////////////////////////
    // Runs CollectGcBias with input Sam file and outputs details and summary files for truth assertion.
    /////////////////////////////////////////////////////////////////////////////
    public void runGcBias (final File input, final File outfile, final File detailsOutfile) throws IOException {
        final String referenceFile = "testdata/picard/metrics/chrMNO.reference.fasta";
        final File pdf = File.createTempFile("test", ".pdf");
        pdf.deleteOnExit();

        final int windowSize = 100;
        final double minGenFraction = 1.0E-5;
        final boolean biSulfiteSeq = false;
        final boolean assumeSorted = false;

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + detailsOutfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "SUMMARY_OUTPUT=" + outfile.getAbsolutePath(),
                "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                "SCAN_WINDOW_SIZE=" + windowSize,
                "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                "LEVEL=ALL_READS",
                "LEVEL=SAMPLE",
                "LEVEL=READ_GROUP",
                "ASSUME_SORTED=" + assumeSorted
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
    }

    /////////////////////////////////////////////////////////////////////////////
    //Used to generate the Sam Record Sets with SamRecordSetBuilder.addPair().
    //testNumber 1: runGcBiasMultiLevelTest, generates records aligning to chrM and chrO
    //testNumber 2: runWindowsComparisonTest, generates records aligning to chrM,N,O.
    /////////////////////////////////////////////////////////////////////////////
    public void setupTest1(final int ID, final String readGroupId, final SAMReadGroupRecord readGroupRecord, final String sample,
                      final String library, final SAMFileHeader header, final SAMRecordSetBuilder setBuilder)
            throws IOException {

        final String separator = ":";
        final int contig1 = 0;
        final int contig2 = 1;
        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);
        readGroupRecord.setPlatformUnit(readGroupId);
        header.addReadGroup(readGroupRecord);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);

        setBuilder.setHeader(header);

        final int max = 800;
        final int min = 1;
        final Random rg = new Random(5);

        //add records that align to chrM and O but not N
        for (int i = 0; i < NUM_READS; i++) {
            final int start = rg.nextInt(max) + min;
            final String newReadName = READ_NAME + separator + ID + separator + i;

            if (i != NUM_READS - 1) {
                setBuilder.addPair(newReadName, contig1, start + ID, start + ID + LENGTH);
            } else {
                setBuilder.addPair(newReadName, contig2, start + ID, start + ID + LENGTH);
            }
        }
    }

    public void setupTest2(final int ID, final String readGroupId, final SAMReadGroupRecord readGroupRecord, final String sample,
                           final String library, final SAMFileHeader header, final SAMRecordSetBuilder setBuilder)
            throws IOException {

        final String separator = ":";
        final int contig1 = 0;
        final int contig2 = 1;
        final int contig3 = 2;
        readGroupRecord.setSample(sample);
        readGroupRecord.setPlatform(platform);
        readGroupRecord.setLibrary(library);
        readGroupRecord.setPlatformUnit(readGroupId);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);

        setBuilder.setHeader(header);

        final int max = 800;
        final int min = 1;
        final Random rg = new Random(5);

        //add records that align to all 3 chr in reference file
        for (int i = 0; i < NUM_READS; i++) {
            final int start = rg.nextInt(max) + min;
            final String newReadName = READ_NAME + separator + ID + separator + i;

            if (i<=NUM_READS/3) {
                setBuilder.addPair(newReadName, contig1, start + ID, start + ID + LENGTH);
            } else if (i< (NUM_READS - (NUM_READS/3))) {
                setBuilder.addPair(newReadName, contig2, start + ID, start + ID + LENGTH);
            } else {
                setBuilder.addPair(newReadName, contig3, start + ID, start + ID + LENGTH);
            }
        }
    }
}
