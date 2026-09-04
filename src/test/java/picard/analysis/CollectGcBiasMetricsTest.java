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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.SAMException;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.SortSam;
import picard.util.RExecutor;
import picard.vcf.VcfTestUtils;

import static picard.analysis.GcBiasMetricsCollector.PerUnitGcBiasMetricsCollector.*;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
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
    private final String REFERENCE_FILE_1 = "testdata/picard/metrics/chrMNO.reference.fasta";

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
        tempSamFileChrM_O = VcfTestUtils.createTemporaryIndexedFile("CollectGcBias", ".bam");
        tempSamFileAllChr = VcfTestUtils.createTemporaryIndexedFile("CollectGcBias", ".bam");

        final File tempSamFileUnsorted = VcfTestUtils.createTemporaryIndexedFile("CollectGcBias", ".bam");

        final SAMFileHeader header = new SAMFileHeader();

        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(dict.toPath()));
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

        final List<SAMRecordSetBuilder> test1Builders = new ArrayList<>();
        test1Builders.add(setBuilder1);
        test1Builders.add(setBuilder2);
        test1Builders.add(setBuilder3);

        final List<SAMRecordSetBuilder> test2Builders = new ArrayList<>();
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

        runGcBias(tempSamFileChrM_O, REFERENCE_FILE_1, outfile, detailsOutfile, false);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS)) { //ALL_READS level
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

        runGcBias(tempSamFileChrM_O, REFERENCE_FILE_1, outfile, detailsOutfile, false);
        runGcBias(tempSamFileAllChr, REFERENCE_FILE_1, allChrOutFile, allChrDetailsOutfile, false);

        final MetricsFile<GcBiasDetailMetrics, Comparable<?>> outputDetails = new MetricsFile<>();
        outputDetails.read(new FileReader(detailsOutfile));
        final List<GcBiasDetailMetrics> details = outputDetails.getMetrics();

        final MetricsFile<GcBiasDetailMetrics, Comparable<?>> outputAllChrDetails = new MetricsFile<>();
        outputAllChrDetails.read(new FileReader(allChrDetailsOutfile));

        int i = 0;

        //Output for the two sam files are only the same for the "All Reads" level
        for (final GcBiasDetailMetrics metrics : outputAllChrDetails.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS)) {
                Assert.assertEquals(metrics.WINDOWS, details.get(i).WINDOWS);
                i++;
            } else {
                break;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////
    // tempSamFileAllChr has reads spread across all three contigs in the chrMNO reference (chrM, chrN,
    // chrO). Splitting it into a "MO" group (chrM + chrO) and an "N" group (chrN) should produce two
    // sections in the output whose ALIGNED_READS sum to the ungrouped total, since every read's contig
    // belongs to exactly one of those groups.
    /////////////////////////////////////////////////////////////////////////////
    @Test
    public void runContigGroupsTest() throws IOException {
        final File summaryOutfile = File.createTempFile("test", ".gc_bias.summary_metrics");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias.detail_metrics");
        summaryOutfile.deleteOnExit();
        detailsOutfile.deleteOnExit();
        final File pdf = File.createTempFile("test", ".pdf");
        pdf.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + tempSamFileAllChr.getAbsolutePath(),
                "OUTPUT=" + detailsOutfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + REFERENCE_FILE_1,
                "SUMMARY_OUTPUT=" + summaryOutfile.getAbsolutePath(),
                "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                "SCAN_WINDOW_SIZE=100",
                "MINIMUM_GENOME_FRACTION=1.0E-5",
                "IS_BISULFITE_SEQUENCED=false",
                "LEVEL=ALL_READS",
                "ASSUME_SORTED=false",
                "CONTIG_GROUP=MO:chrM,chrO",
                "CONTIG_GROUP=N:chrN"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<GcBiasSummaryMetricsWithContigGroup, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(summaryOutfile));

        long moAligned = 0;
        long nAligned = 0;
        long moClusters = 0;
        long nClusters = 0;
        boolean sawMO = false;
        boolean sawN = false;
        for (final GcBiasSummaryMetricsWithContigGroup metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.ACCUMULATION_LEVEL, ACCUMULATION_LEVEL_ALL_READS);
            if ("MO".equals(metrics.CONTIG_GROUP)) {
                sawMO = true;
                moAligned = metrics.ALIGNED_READS;
                moClusters = metrics.TOTAL_CLUSTERS;
            } else if ("N".equals(metrics.CONTIG_GROUP)) {
                sawN = true;
                nAligned = metrics.ALIGNED_READS;
                nClusters = metrics.TOTAL_CLUSTERS;
            } else {
                Assert.fail("Unexpected CONTIG_GROUP: " + metrics.CONTIG_GROUP);
            }
        }
        Assert.assertTrue(sawMO);
        Assert.assertTrue(sawN);
        Assert.assertTrue(moAligned > 0);
        Assert.assertTrue(nAligned > 0);

        final File ungroupedSummary = File.createTempFile("test", ".gc_bias.summary_metrics");
        final File ungroupedDetails = File.createTempFile("test", ".gc_bias.detail_metrics");
        ungroupedSummary.deleteOnExit();
        ungroupedDetails.deleteOnExit();
        runGcBias(tempSamFileAllChr, REFERENCE_FILE_1, ungroupedSummary, ungroupedDetails, false);
        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> ungroupedOutput = new MetricsFile<>();
        ungroupedOutput.read(new FileReader(ungroupedSummary));
        long totalAligned = 0;
        long totalClusters = 0;
        for (final GcBiasSummaryMetrics metrics : ungroupedOutput.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS)) {
                totalAligned = metrics.ALIGNED_READS;
                totalClusters = metrics.TOTAL_CLUSTERS;
                // The ungrouped run never used CONTIG_GROUP, so its output is written as plain
                // GcBiasSummaryMetrics (no CONTIG_GROUP column at all), unlike the grouped run above.
                Assert.assertFalse(metrics instanceof GcBiasSummaryMetricsWithContigGroup);
            }
        }
        Assert.assertEquals(moAligned + nAligned, totalAligned);
        // TOTAL_CLUSTERS is a property of the accumulation level, not of a contig group, so it stays
        // genome-wide and identical across the groups rather than being split between them.
        Assert.assertEquals(moClusters, totalClusters);
        Assert.assertEquals(nClusters, totalClusters);

        // The detail file should carry the CONTIG_GROUP column too, with 101 GC bins per group.
        final MetricsFile<GcBiasDetailMetricsWithContigGroup, Comparable<?>> outputDetails = new MetricsFile<>();
        outputDetails.read(new FileReader(detailsOutfile));
        long moWindowRows = 0;
        long nWindowRows = 0;
        for (final GcBiasDetailMetricsWithContigGroup detail : outputDetails.getMetrics()) {
            if ("MO".equals(detail.CONTIG_GROUP)) moWindowRows++;
            else if ("N".equals(detail.CONTIG_GROUP)) nWindowRows++;
        }
        Assert.assertEquals(moWindowRows, 101);
        Assert.assertEquals(nWindowRows, 101);
    }

    /////////////////////////////////////////////////////////////////////////////
    // CONTIG_GROUP specs parse into a contig -> group names map, and a contig may belong to
    // more than one group.
    /////////////////////////////////////////////////////////////////////////////
    @Test
    public void testParseContigGroups() {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        for (final String name : new String[]{"chr1", "chr2", "chrX", "chrY"}) {
            dict.addSequence(new SAMSequenceRecord(name, 1000));
        }

        final Map<String, List<String>> contigToGroupNames = GcBiasUtils.parseContigGroups(
                Arrays.asList("autosomes:chr1, chr2", "sex:chrX,chrY", "everything:chr1,chr2,chrX,chrY"), dict);

        Assert.assertEquals(contigToGroupNames.get("chr1"), Arrays.asList("autosomes", "everything"));
        Assert.assertEquals(contigToGroupNames.get("chr2"), Arrays.asList("autosomes", "everything"));
        Assert.assertEquals(contigToGroupNames.get("chrX"), Arrays.asList("sex", "everything"));
        Assert.assertEquals(contigToGroupNames.get("chrY"), Arrays.asList("sex", "everything"));
    }

    @DataProvider(name = "badContigGroups")
    public Object[][] badContigGroups() {
        return new Object[][]{
                {"no colon separator", "autosomes chr1"},
                {"empty group name", ":chr1"},
                {"no contigs listed", "empty:"},
                {"contig absent from the sequence dictionary", "bad:chrDoesNotExist"},
        };
    }

    /////////////////////////////////////////////////////////////////////////////
    // Malformed CONTIG_GROUP specs, and contigs absent from the sequence dictionary, should be
    // hard errors rather than being silently ignored.
    /////////////////////////////////////////////////////////////////////////////
    @Test(dataProvider = "badContigGroups", expectedExceptions = PicardException.class)
    public void testBadContigGroupsThrow(final String description, final String spec) {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 1000));

        GcBiasUtils.parseContigGroups(Collections.singletonList(spec), dict);
    }

    /////////////////////////////////////////////////////////////////////////////
    // The same group name given twice is a user error, not a silent merge.
    /////////////////////////////////////////////////////////////////////////////
    @Test(expectedExceptions = PicardException.class)
    public void testDuplicateContigGroupNameThrows() {
        final SAMSequenceDictionary dict = new SAMSequenceDictionary();
        dict.addSequence(new SAMSequenceRecord("chr1", 1000));
        dict.addSequence(new SAMSequenceRecord("chr2", 1000));

        GcBiasUtils.parseContigGroups(Arrays.asList("dup:chr1", "dup:chr2"), dict);
    }

    /////////////////////////////////////////////////////////////////////////////
    // Writes the setBuilders to a SAMFileWriter and sorts the sam.
    // Takes in a list of SAMRecordSetBuilders because of the multi-level collection: setBuilders cannot take in more than one read group
    // or library or sample, so there are separate ones for each type when testing multi-level collection.
    /////////////////////////////////////////////////////////////////////////////
    public File build (final List<SAMRecordSetBuilder> setBuilder, final File unsortedSam, final SAMFileHeader header) throws IOException {
        final File sortedSam = VcfTestUtils.createTemporaryIndexedFile("CollectGcBias", ".bam");

        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeBAMWriter(header, false, unsortedSam);

        for (final SAMRecordSetBuilder subSetBuilder : setBuilder) {
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
    public void runGcBias (final File input, final String referenceFile, final File summaryOutfile, final File detailsOutfile,
                           final boolean nonDups) throws IOException {
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
                "SUMMARY_OUTPUT=" + summaryOutfile.getAbsolutePath(),
                "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                "SCAN_WINDOW_SIZE=" + windowSize,
                "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                "LEVEL=ALL_READS",
                "LEVEL=SAMPLE",
                "LEVEL=READ_GROUP",
                "ASSUME_SORTED=" + assumeSorted,
                "ALSO_IGNORE_DUPLICATES=" + nonDups
        };
        runPicardCommandLine(args);
    }

    /**
     * Compares metric's results by summary files without duplicates.
     * @throws IOException
     */
    @Test
    public void runNonDupsComparisonTest() throws IOException {
        final File inputFileWithDuplicates = new File("testdata/picard/metrics/chrMReads.sam");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias_detail_metrics");
        final File summaryOutfile = File.createTempFile("test", ".gc_bias_summary_metrics");
        detailsOutfile.deleteOnExit();
        summaryOutfile.deleteOnExit();

        runGcBias(inputFileWithDuplicates, CHR_M_REFERENCE.getAbsolutePath(), summaryOutfile, detailsOutfile, true);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> outputSummary = new MetricsFile<>();
        outputSummary.read(new FileReader(summaryOutfile));

        for (final GcBiasSummaryMetrics summary : outputSummary.getMetrics()) {
            if (summary.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS) && summary.READS_USED.equals(READS_USED_UNIQUE)) { //ALL_READS level for case without duplicates
                Assert.assertEquals(summary.TOTAL_CLUSTERS, 3);
                Assert.assertEquals(summary.ALIGNED_READS, 3);
                Assert.assertEquals(summary.AT_DROPOUT, 79.180328);
                Assert.assertEquals(summary.GC_DROPOUT, 12.28901);
                Assert.assertEquals(summary.GC_NC_0_19, 0.0);
                Assert.assertEquals(summary.GC_NC_20_39, 0.0);
                Assert.assertEquals(summary.GC_NC_40_59, 1.246783);
                Assert.assertEquals(summary.GC_NC_60_79, 0.0);
                Assert.assertEquals(summary.GC_NC_80_100, 0.0);
            }
            if (summary.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS) && summary.READS_USED.equals(READS_USED_ALL)) { //ALL_READS level
                Assert.assertEquals(summary.TOTAL_CLUSTERS, 5);
                Assert.assertEquals(summary.ALIGNED_READS, 5);
                Assert.assertEquals(summary.AT_DROPOUT, 79.180328);
                Assert.assertEquals(summary.GC_DROPOUT, 10.37037);
                Assert.assertEquals(summary.GC_NC_0_19, 0.0);
                Assert.assertEquals(summary.GC_NC_20_39, 0.0);
                Assert.assertEquals(summary.GC_NC_40_59, 1.246783);
                Assert.assertEquals(summary.GC_NC_60_79, 0.0);
                Assert.assertEquals(summary.GC_NC_80_100, 0.0);
            }
        }
    }

    /**
     * If SAM/BAM file with '*' in SEQ field omit this read.
     */
    @Test
    public void runCheckingNoSEQTest() throws IOException {
        final File input = new File("testdata/picard/metrics/chrM_NO_SEQ.sam");
        final File summaryOutfile = File.createTempFile("test", ".gc_bias.summary_metrics");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias.detail_metrics");
        summaryOutfile.deleteOnExit();
        detailsOutfile.deleteOnExit();

        runGcBias(input, CHR_M_REFERENCE.getAbsolutePath(), summaryOutfile, detailsOutfile, false);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(summaryOutfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS)) { //ALL_READS level
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 0);
                Assert.assertEquals(metrics.ALIGNED_READS, 1);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.682453);
                Assert.assertEquals(metrics.GC_DROPOUT, 14.693382);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.0);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.246783);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            }
        }
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

    @Test
    public void testChartFailureGATKLite () throws IOException {
        final PrintStream stderr = System.err;
        final String gatkLiteDockerProperty = System.getProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR);

        try {
            final ByteArrayOutputStream stdoutCapture = new ByteArrayOutputStream();
            System.setErr(new PrintStream(stdoutCapture));

            System.setProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR, "true");

            final File input = new File("testdata/picard/metrics/chrM_NO_SEQ.sam");
            final File summaryOutfile = File.createTempFile("test", ".gc_bias.summary_metrics");
            final File detailsOutfile = File.createTempFile("test", ".gc_bias.detail_metrics");
            summaryOutfile.deleteOnExit();
            detailsOutfile.deleteOnExit();

            final File pdf = File.createTempFile("test", ".pdf");
            pdf.deleteOnExit();

            final int windowSize = 100;
            final double minGenFraction = 1.0E-5;
            final boolean biSulfiteSeq = false;
            final boolean assumeSorted = false;

            final String[] args = new String[]{
                    "INPUT=" + input.getAbsolutePath(),
                    "OUTPUT=" + detailsOutfile.getAbsolutePath(),
                    "REFERENCE_SEQUENCE=" + CHR_M_REFERENCE.getAbsolutePath(),
                    "SUMMARY_OUTPUT=" + summaryOutfile.getAbsolutePath(),
                    "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                    "SCAN_WINDOW_SIZE=" + windowSize,
                    "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                    "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                    "LEVEL=ALL_READS",
                    "LEVEL=SAMPLE",
                    "LEVEL=READ_GROUP",
                    "ASSUME_SORTED=" + assumeSorted,
                    "ALSO_IGNORE_DUPLICATES=" + false
            };
            Assert.assertEquals(runPicardCommandLine(args), 1);

            Assert.assertTrue(stdoutCapture.toString().contains("The histogram file cannot be written because it requires R, which is not available in the GATK Lite Docker image."));  
        }
        finally {
            System.setErr(stderr);
            if(gatkLiteDockerProperty != null) {
                System.setProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR, gatkLiteDockerProperty);
            }
            else{
                System.clearProperty(RExecutor.GATK_LITE_DOCKER_ENV_VAR);
            } 
        }
    }

    @Test
    public void testWithoutChartOutput() throws IOException {
        final File input = new File("testdata/picard/metrics/chrM_NO_SEQ.sam");
        final File summaryOutfile = File.createTempFile("test", ".gc_bias.summary_metrics");
        final File detailsOutfile = File.createTempFile("test", ".gc_bias.detail_metrics");
        summaryOutfile.deleteOnExit();
        detailsOutfile.deleteOnExit();

        final int windowSize = 100;
        final double minGenFraction = 1.0E-5;
        final boolean biSulfiteSeq = false;
        final boolean assumeSorted = false;

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + detailsOutfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + CHR_M_REFERENCE.getAbsolutePath(),
                "SUMMARY_OUTPUT=" + summaryOutfile.getAbsolutePath(),
                "SCAN_WINDOW_SIZE=" + windowSize,
                "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                "LEVEL=ALL_READS",
                "LEVEL=SAMPLE",
                "LEVEL=READ_GROUP",
                "ASSUME_SORTED=" + assumeSorted,
                "ALSO_IGNORE_DUPLICATES=" + false
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(summaryOutfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals(ACCUMULATION_LEVEL_ALL_READS)) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 0);
                Assert.assertEquals(metrics.ALIGNED_READS, 1);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.682453);
                Assert.assertEquals(metrics.GC_DROPOUT, 14.693382);
                Assert.assertEquals(metrics.GC_NC_0_19, 0.0);
                Assert.assertEquals(metrics.GC_NC_20_39, 0.0);
                Assert.assertEquals(metrics.GC_NC_40_59, 1.246783);
                Assert.assertEquals(metrics.GC_NC_60_79, 0.0);
                Assert.assertEquals(metrics.GC_NC_80_100, 0.0);
            }
        }
    }
}
