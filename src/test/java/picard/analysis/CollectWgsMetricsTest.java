/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

package picard.analysis;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
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

    private final static int READ_PAIR_DISTANCE = 99;
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
                {tempSamFile, outfile, referenceFile, "false"},
                {tempSamFile, outfile, referenceFile, "true"},
        };
    }

    @Test(dataProvider = "wgsDataProvider")
    public void testMetricsFromWGS(final File input, final File outfile, final String referenceFile,
            final String useFastAlgorithm) throws IOException {
        outfile.deleteOnExit();
        final int sampleSize = 1000;

        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "SAMPLE_SIZE=" + sampleSize,
                "USE_FAST_ALGORITHM=" + useFastAlgorithm
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
        final File tempSamIndex = new File(tempSamFile.getAbsolutePath().replace("bam", "bai"));
        final File tempSamFileUnsorted = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        tempSamFileUnsorted.deleteOnExit();
        tempSamIndex.deleteOnExit();
        tempSamFile.deleteOnExit();
        final File sortedSamIdx = new File(TEST_DIR, tempSamFile.getName() + ".idx");
        sortedSamIdx.deleteOnExit();

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
            setBuilder.addPair(newReadName, 0, start + ID, start + ID + READ_PAIR_DISTANCE);
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
        outfile = File.createTempFile("testWgsMetrics", ".txt", TEST_DIR);
        outfile.deleteOnExit();
    }

    @DataProvider(name = "wgsAlgorithm")
    public Object[][] wgsAlgorithm() {
        return new Object[][] {
                {"false"},
                {"true"},
        };
    }
    @Test(dataProvider = "wgsAlgorithm")
    public void testLargeIntervals(final String useFastAlgorithm) throws IOException {
        final File input = new File(TEST_DIR, "forMetrics.sam");
        final File outfile = File.createTempFile("test", ".wgs_metrics");
        final File ref = new File(TEST_DIR, "merger.fasta");
        final File intervals = new File(TEST_DIR, "largeIntervals.interval_list");
        final int sampleSize = 1000;
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + ref.getAbsolutePath(),
                "INTERVALS=" + intervals.getAbsolutePath(),
                "SAMPLE_SIZE=" + sampleSize,
                "USE_FAST_ALGORITHM=" + useFastAlgorithm
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final CollectWgsMetrics.WgsMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.GENOME_TERRITORY, 404);
            Assert.assertEquals(metrics.PCT_EXC_MAPQ, 0.271403);
            Assert.assertEquals(metrics.PCT_EXC_DUPE, 0.182149);
            Assert.assertEquals(metrics.PCT_EXC_UNPAIRED, 0.091075);
        }
    }

    @Test(dataProvider = "wgsAlgorithm")
    public void testExclusions(final String useFastAlgorithm) throws IOException {
        final File reference = new File("testdata/picard/sam/merger.fasta");
        final File tempSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        tempSamFile.deleteOnExit();


        final SAMRecordSetBuilder setBuilder = CollectWgsMetricsTestUtils.createTestSAMBuilder(reference, READ_GROUP_ID, SAMPLE, PLATFORM, LIBRARY);

        setBuilder.setReadLength(10);

        int expectedSingletonCoverage = 0;

        expectedSingletonCoverage += 13;
        setBuilder.addPair("overlappingReads", 0, 2, 5, false, false, "10M", "10M", true, false, 30);

        expectedSingletonCoverage += 2 * 5; // 5 bases for each mate are good (see AAA!!!AA!! below).
        setBuilder.addPair("poorQualityReads", 1, 2, 20, false, false, "10M", "10M", true, false, -1);

        for(int i = 1; i < 5; i++) {
            setBuilder.addPair("deepStack-" + i, 2, 2, 20, false, false, "10M", "10M", true, false, 30);
        }

        // modify quality of reads
        setBuilder.getRecords().stream()
                .filter(samRecord -> samRecord.getReadName().equals("poorQualityReads"))
                .forEach(record -> record.setBaseQualityString("AAA!!!AA!!"));

        setBuilder.getSamReader();

        // Write SAM file
        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, tempSamFile);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }
        writer.close();

        // create output files for tests
        final File outfile = File.createTempFile("testWgsMetrics", ".txt");
        outfile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT="  + tempSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=true",
                "COVERAGE_CAP=3",
                "USE_FAST_ALGORITHM=" + useFastAlgorithm
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));
        final CollectWgsMetrics.WgsMetrics metrics = output.getMetrics().get(0);

        final Histogram<Integer> highQualityDepthHistogram = output.getAllHistograms().get(0);
        final Histogram<Integer> baseQHistogram = output.getAllHistograms().get(1);

        Assert.assertEquals((long) highQualityDepthHistogram.getSumOfValues(), metrics.GENOME_TERRITORY);
        Assert.assertEquals((long) highQualityDepthHistogram.get(1).getValue(), expectedSingletonCoverage);
        Assert.assertEquals((long) highQualityDepthHistogram.get(3).getValue(), 2*10);
    }

    @Test(dataProvider = "wgsAlgorithm")
    public void testPoorQualityBases(final String useFastAlgorithm) throws IOException {
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final File testSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        testSamFile.deleteOnExit();

        /**
         *  Our test SAM looks as follows:
         *
         *   ----------   <- reads with great base qualities (60) ->  ----------
         *   ----------                                               ----------
         *   ----------                                               ----------
         *   **********   <- reads with poor base qualities (10) ->   **********
         *   **********                                               **********
         *   **********                                               **********
         *
         *  We exclude half of the bases because they are low quality.
         *  We do not exceed the coverage cap (3), thus none of the bases is excluded as such.
         *
         */

        final SAMRecordSetBuilder setBuilder = CollectWgsMetricsTestUtils.createTestSAMBuilder(reference, READ_GROUP_ID, SAMPLE, PLATFORM, LIBRARY);
        setBuilder.setReadLength(10);
        for (int i = 0; i < 3; i++){
            setBuilder.addPair("GreatBQRead:" + i, 0, 1, 30, false, false, "10M", "10M", false, true, 40);
        }

        for (int i = 0; i < 3; i++){
            setBuilder.addPair("PoorBQRead:" + i, 0, 1, 30, false, false, "10M", "10M", false, true, 10);
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, testSamFile);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }

        writer.close();

        final File outfile = File.createTempFile("testExcludedBases-metrics", ".txt");
        outfile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT="  + testSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=true",
                "COVERAGE_CAP=3",
                "USE_FAST_ALGORITHM=" + useFastAlgorithm
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final CollectWgsMetrics.WgsMetrics metrics = output.getMetrics().get(0);

        Assert.assertEquals(metrics.PCT_EXC_BASEQ, 0.5);
        Assert.assertEquals(metrics.PCT_EXC_CAPPED, 0.0);
    }

    @Test(dataProvider = "wgsAlgorithm")
    public void testGiantDeletion(final String useFastAlgorithm) throws IOException {
        final File reference = new File("testdata/picard/quality/chrM.reference.fasta");
        final File testSamFile = File.createTempFile("CollectWgsMetrics", ".bam", TEST_DIR);
        testSamFile.deleteOnExit();

        final SAMRecordSetBuilder setBuilder = CollectWgsMetricsTestUtils.createTestSAMBuilder(reference, READ_GROUP_ID, SAMPLE, PLATFORM, LIBRARY);
        setBuilder.setReadLength(10);
        for (int i = 0; i < 3; i++){
            setBuilder.addPair("Read:" + i, 0, 1, 30, false, false, "5M10000D5M", "1000D10M", false, true, 40);
        }

        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(setBuilder.getHeader(), false, testSamFile);

        for (final SAMRecord record : setBuilder) {
            writer.addAlignment(record);
        }

        writer.close();

        final File outfile = File.createTempFile("testGiantDeletion", ".txt");
        outfile.deleteOnExit();

        final String[] args = new String[] {
                "INPUT="  + testSamFile.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "INCLUDE_BQ_HISTOGRAM=true",
                "COVERAGE_CAP=3",
                "USE_FAST_ALGORITHM=" + useFastAlgorithm
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetrics.WgsMetrics, Integer> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        final CollectWgsMetrics.WgsMetrics metrics = output.getMetrics().get(0);

        Assert.assertEquals(metrics.PCT_EXC_BASEQ, 0.0);
        Assert.assertEquals(metrics.PCT_EXC_CAPPED, 0.0);
    }
}
