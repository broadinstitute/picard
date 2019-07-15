package picard.analysis;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.BufferedLineReader;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.SortSam;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import static picard.analysis.GcBiasMetricsCollector.PerUnitGcBiasMetricsCollector.ACCUMULATION_LEVEL_ALL_READS;

/**
 * Tests the two default "programs" that have tests in CollectMultipleMetrics
 *
 * @author Yossi farjoun
 */
public class CollectMultipleMetricsTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    public String getCommandLineProgramName() {
        return CollectMultipleMetrics.class.getSimpleName();
    }

    @Test
    public void testAlignmentSummaryViaMultipleMetrics() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile = File.createTempFile("alignmentMetrics", "");
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "PROGRAM=null",
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectAlignmentSummaryMetrics.name(),
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectInsertSizeMetrics.name(),
                "EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH= 58",
                "EXTRA_ARGUMENT=CollectInsertSizeMetrics::METRIC_ACCUMULATION_LEVEL=LIBRARY",
                "EXTRA_ARGUMENT=CollectInsertSizeMetrics::METRIC_ACCUMULATION_LEVEL=READ_GROUP",

        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<InsertSizeMetrics, Comparable<?>> outputISM = new MetricsFile<>();
        outputISM.read(new FileReader(outfile + ".insert_size_metrics"));
        Assert.assertEquals(outputISM.getMetrics().size(), 3);

        for (final InsertSizeMetrics metrics : outputISM.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_INSERT_SIZE, 40D);
        }

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile + ".alignment_summary_metrics"));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
                case FIRST_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 7);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 3);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 59);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 19.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 303);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*58D/303D*/0.191419);
                    Assert.assertEquals(metrics.BAD_CYCLES, 19);
                    break;
                case SECOND_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 9);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 7);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 239);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 707);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*19D/707D*/0.026874);
                    Assert.assertEquals(metrics.BAD_CYCLES, 3);
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 18);
                    Assert.assertEquals(metrics.PF_READS, 16);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 10);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 298);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 1010);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*77D/1010D*/0.076238);
                    Assert.assertEquals(metrics.BAD_CYCLES, 22);
                    break;
                case UNPAIRED:
                default:
                    Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @DataProvider
    Object[][] extraArgumentValue() {
                List<Object[]> tests = new ArrayList<>();

                // this is actually legal after conversion to the non-legacy parser
                if (CommandLineProgram.useLegacyParser(getClass())) {
                    tests.add(new Object[]{"CollectInsertSizeMetrics::OUTPUT =hi.out"});
                }

                tests.add(new Object[]{"BLAH::HISTOGRAM_WIDTH=58"});
                tests.add(new Object[]{"QualityScoreDistribution::HISTOGRAM_WIDTH=58"});
                tests.add(new Object[]{"CollectInsertSizeMetrics::BLAH=58"});
                tests.add(new Object[]{"CollectInsertSizeMetrics:HISTOGRAM_WIDTH=58"});
                tests.add(new Object[]{"CollectInsertSizeMetrics::HISTOGRAM_WIDTH=5a8"});
                tests.add(new Object[]{"CollectInsertSizeMetrics::HISTOGRAM_WIDTH="});
                tests.add(new Object[]{"CollectInsertSizeMetrics::HISTOGRAM_WIDTH=hello"});
                tests.add(new Object[]{"CollectInsertSizeMetrics::REFERENCE=whyNot.fasta"});
                tests.add(new Object[]{"CollectInsertSizeMetrics::OUTPUT="});

        return tests.toArray(new Object[0][]);
    }

    @Test(expectedExceptions = {SAMException.class, CommandLineException.class}, dataProvider = "extraArgumentValue")
    public void testMultipleMetricsFailing(final String extra) throws IOException {

        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile = File.createTempFile("alignmentMetrics", "");
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "PROGRAM=null",
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectAlignmentSummaryMetrics.name(),
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectInsertSizeMetrics.name(),
                "EXTRA_ARGUMENT=" + extra
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
    }

    @Test
    public void testInsertSize() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File outfile = File.createTempFile("test", "");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File pdf = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "PROGRAM=null",
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectAlignmentSummaryMetrics.name(),
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectInsertSizeMetrics.name()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);
        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile + ".insert_size_metrics"));
        for (final InsertSizeMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.PAIR_ORIENTATION.name(), "FR");
            if (metrics.LIBRARY == null) {  // SAMPLE or ALL_READS level
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 13);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            } else if (metrics.LIBRARY.equals("Solexa-41753")) { // one LIBRARY and one READ_GROUP
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.READ_PAIRS, 2);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);

            } else if (metrics.LIBRARY.equals("Solexa-41748") && metrics.READ_GROUP == null) {
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 40);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 9);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            } else if (metrics.LIBRARY.equals("Solexa-41734") && metrics.READ_GROUP == null) {
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 26);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 9);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);
            } else if (metrics.READ_GROUP.equals("62A79AAXX100907.7")) {
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 4);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);
            } else if (metrics.READ_GROUP.equals("62A79AAXX100907.6")) {
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 38);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 5);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 9);
            } else if (metrics.READ_GROUP.equals("62A79AAXX100907.5")) {
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 1);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);
            } else if (metrics.READ_GROUP.equals("62A79AAXX100907.3")) {
                Assert.assertEquals((int) metrics.MEDIAN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.READ_PAIRS, 1);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);
            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }

    @Test
    public void testRnaSeqMetricsViaMultipleMetrics() throws Exception {
        final File input = new File(TEST_DATA_DIR, "rna_seq_metrics.sam");
        final File outfile = File.createTempFile("tmp.", ".rna_metrics");
        final File intervals = new File(TEST_DATA_DIR, "rna_seq_metrics.interval_list");
        final File refflat = new File(TEST_DATA_DIR, "rna_seq_metrics.refflat");
        final Set<String> ignoreSequence = new HashSet<>();
        ignoreSequence.add("chrM");
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "INTERVALS=" + intervals.getAbsolutePath(),
                "REF_FLAT=" + refflat.getAbsolutePath(),
                "IGNORE_SEQUENCE=" + StringUtils.join(ignoreSequence, ','),
                "PROGRAM=" + CollectMultipleMetrics.Program.RnaSeqMetrics.name(),
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<RnaSeqMetrics, Comparable<?>> output = new MetricsFile<RnaSeqMetrics, Comparable<?>>();
        output.read(new FileReader(outfile + ".rna_metrics"));
        final RnaSeqMetrics metrics = output.getMetrics().get(0);

        Assert.assertEquals(metrics.PF_ALIGNED_BASES, 396);
        Assert.assertEquals(metrics.PF_BASES, 432);
        Assert.assertEquals(metrics.RIBOSOMAL_BASES.longValue(), 108L);
        Assert.assertEquals(metrics.CODING_BASES, 136);
        Assert.assertEquals(metrics.UTR_BASES, 51);
        Assert.assertEquals(metrics.INTRONIC_BASES, 50);
        Assert.assertEquals(metrics.INTERGENIC_BASES, 51);
        Assert.assertEquals(metrics.CORRECT_STRAND_READS, 3);
        Assert.assertEquals(metrics.INCORRECT_STRAND_READS, 4);
        Assert.assertEquals(metrics.IGNORED_READS, 1);
        Assert.assertEquals(metrics.NUM_R1_TRANSCRIPT_STRAND_READS, 1);
        Assert.assertEquals(metrics.NUM_R2_TRANSCRIPT_STRAND_READS, 2);
        Assert.assertEquals(metrics.NUM_UNEXPLAINED_READS, 2);
        Assert.assertEquals(metrics.PCT_R1_TRANSCRIPT_STRAND_READS, 0.333333);
        Assert.assertEquals(metrics.PCT_R2_TRANSCRIPT_STRAND_READS, 0.666667);
    }

    @Test //test all gcBias collection levels
    public void testGcBiasMetrics() throws IOException {
        runGcTest(tempSamFile);
    }

    public void runGcTest(final File input) throws IOException {
        final File outfile = File.createTempFile("test", "");
        outfile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + CHR_M_REFERENCE.getAbsolutePath(),
                "METRIC_ACCUMULATION_LEVEL=" + MetricAccumulationLevel.ALL_READS.name(),
                "PROGRAM=null",
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectAlignmentSummaryMetrics.name(),
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectInsertSizeMetrics.name(),
                "PROGRAM=" + CollectMultipleMetrics.Program.CollectGcBiasMetrics.name()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<GcBiasSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(outfile + ".gc_bias.summary_metrics"));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS)) { //ALL_READS level
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 300);
                Assert.assertEquals(metrics.ALIGNED_READS, 600);
                Assert.assertEquals(metrics.AT_DROPOUT, 7.234062);
                Assert.assertEquals(metrics.GC_DROPOUT, 4.086217);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 1.06826);
                Assert.assertEquals(metrics.GC_NC_40_59, 0.987036);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }

    //gcBias multi level collector test creates a sam file from chrM for testing purposes
    //more variables needed for gcbias test to create temp sam file
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
    private final File dict = CHR_M_DICT;
    private File tempSamFile;
    private final SAMRecordSetBuilder setBuilder1 =
            new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    private final SAMRecordSetBuilder setBuilder2 =
            new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    private final SAMRecordSetBuilder setBuilder3 =
            new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);

    private final SAMReadGroupRecord readGroupRecord1 = new SAMReadGroupRecord(readGroupId1);
    private final SAMReadGroupRecord readGroupRecord2 = new SAMReadGroupRecord(readGroupId2);
    private final SAMReadGroupRecord readGroupRecord3 = new SAMReadGroupRecord(readGroupId3);

    //create a samfile with different samples, read groups and libraries that overlap for testing.
    @BeforeTest
    void setupBuilder() throws IOException {
        final int numReads = 100;
        final String flowCellBarcode = "TESTBARCODE";
        tempSamFile = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        final File tempSamIndex = new File(tempSamFile.getAbsolutePath().replace("bam", "bai"));
        final File tempSamFileUnsorted = File.createTempFile("CollectGcBias", ".bam", TEST_DIR);
        tempSamFileUnsorted.deleteOnExit();
        tempSamIndex.deleteOnExit();
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
        setup(numReads, flowCellBarcode, 1, readGroupId1, readGroupRecord1, sample1, library1, header, setBuilder1); //Sample 1, Library 1, RG 1
        setup(numReads, flowCellBarcode, 2, readGroupId2, readGroupRecord2, sample1, library2, header, setBuilder2); //Sample 1, Library 2, RG 2
        setup(numReads, flowCellBarcode, 3, readGroupId3, readGroupRecord3, sample2, library3, header, setBuilder3); //Sample 2, Library 3, RG 3

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
        final String[] args = new String[]{"INPUT=" + tempSamFileUnsorted.getAbsolutePath(),
                "OUTPUT=" + tempSamFile.getAbsolutePath(),
                "SORT_ORDER=coordinate"};
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
            setBuilder.addPair(newReadName, 0, start + ID, start + ID + 99);
        }
    }
}
