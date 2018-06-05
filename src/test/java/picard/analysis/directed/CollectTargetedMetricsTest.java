package picard.analysis.directed;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.analysis.TheoreticalSensitivity;
import picard.analysis.TheoreticalSensitivityMetrics;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.SortSam;
import picard.util.TestNGUtil;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class CollectTargetedMetricsTest extends CommandLineProgramTest {

    private final String TEST_DATA_DIR = "testdata/picard/quality/";
    private final File dict = new File(TEST_DATA_DIR+"chrM.reference.dict");
    private File tempSamFile;
    private File outfile;
    private File tsOutfile; // Theoretical sensitivity file
    private File perTargetOutfile;
    private static final int LENGTH = 99;
    private static final int RANDOM_SEED = 51;

    final String referenceFile = TEST_DATA_DIR + "chrM.reference.fasta";
    final String emptyIntervals = TEST_DATA_DIR + "chrM.empty.interval_list";
    final String singleIntervals = TEST_DATA_DIR + "chrM.single.interval_list";


    private static final String sample = "TestSample1";
    private static final String readGroupId = "TestReadGroup1";
    private static final String platform = "ILLUMINA";
    private static final String library = "TestLibrary1";
    private static final int numReads = 40000;

    @Override
    public String getCommandLineProgramName() {
        return CollectTargetedPcrMetrics.class.getSimpleName();
    }

    //create a samfile with 40000 reads for testing whether a cap is found.
    @BeforeTest
    void setupBuilder() throws IOException {
        final String readName = "TESTBARCODE";

        //Create Sam Files
        tempSamFile = VcfTestUtils.createTemporaryIndexedFile("CollectTargetedMetrics", ".bam");
        final File tempSamFileUnsorted = VcfTestUtils.createTemporaryIndexedFile("CollectTargetedMetrics", ".bam");

        final SAMFileHeader header = new SAMFileHeader();

        //Check that dictionary file is readable and then set header dictionary
        try {
            header.setSequenceDictionary(SAMSequenceDictionaryExtractor.extractDictionary(dict.toPath()));
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
        setBuilder.setRandomSeed(RANDOM_SEED);
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
        tsOutfile = File.createTempFile("test", ".TheoreticalSensitivityMetrics");
        outfile.deleteOnExit();
        perTargetOutfile.deleteOnExit();
        tsOutfile.deleteOnExit();
    }

    @DataProvider(name = "targetedIntervalDataProvider")
    public Object[][] targetedIntervalDataProvider() {

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

        final MetricsFile<TargetedPcrMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outfile));

        for (final TargetedPcrMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_READS, numReads * 2);
            Assert.assertEquals(metrics.HET_SNP_SENSITIVITY, .997972, .02);
        }
    }

    @DataProvider(name = "theoreticalSensitivityDataProvider")
    public Object[][] theoreticalSensitivityDataProvider() {

        return new Object[][] {
                // This test is primarily used as an integration test since theoretical sensitivity doesn't converge
                // well with a sample size of 2000 rather than 10,000.  The sample size is set so low as to prevent
                // the tests from taking too long to run.
                {tempSamFile, outfile, tsOutfile, perTargetOutfile, referenceFile, singleIntervals, 2000,
                        Arrays.asList(0.01, 0.05, 0.10,  0.30,  0.50), // Allele fraction
                        Arrays.asList(0.01, 0.52, 0.93,  0.99,  0.99)  // Expected sensitivity
                }
        };
    }

    @Test(dataProvider = "theoreticalSensitivityDataProvider")
    public void runCollectTargetedMetricsTheoreticalSensitivityTest(final File input, final File outfile, final File tsOutfile, final File perTargetOutfile, final String referenceFile,
                                              final String targetIntervals, final int sampleSize, final List<Double> alleleFractions, final List<Double> expectedSensitivities) throws IOException {

        final String[] args = new String[] {
                "TARGET_INTERVALS=" + targetIntervals,
                "INPUT=" + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "PER_TARGET_COVERAGE=" + perTargetOutfile.getAbsolutePath(),
                "LEVEL=ALL_READS",
                "AMPLICON_INTERVALS=" + targetIntervals,
                "THEORETICAL_SENSITIVITY_OUTPUT=" + tsOutfile.getAbsolutePath(),
                "SAMPLE_SIZE=" + sampleSize
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<TheoreticalSensitivityMetrics, Double> output = new MetricsFile<>();
        output.read(new FileReader(tsOutfile.getAbsolutePath()));

        // Create map of theoretical sensitivities read from file
        final Iterator<TheoreticalSensitivityMetrics> metrics = output.getMetrics().iterator();
        final Map<Double, Double> map = new HashMap<>();
        while (metrics.hasNext()) {
            TheoreticalSensitivityMetrics m = metrics.next();
            map.put(m.ALLELE_FRACTION, m.THEORETICAL_SENSITIVITY);
        }

        // Compare theoretical sensitivities created by CollectTargetedMetrics
        // with those in the test data provider.
        final Iterator<Double> alleleFraction = alleleFractions.iterator();
        final Iterator<Double> expectedSensitivity = expectedSensitivities.iterator();
        while (alleleFraction.hasNext() || expectedSensitivity.hasNext()) {
            // Provide wiggle room of 1% in the comparisons
            Assert.assertEquals(map.get(alleleFraction.next()), expectedSensitivity.next(), 0.01);
        }
    }

    @Test()
    public void testRawBqDistributionWithSoftClips() throws IOException {
        final String input = TEST_DATA_DIR + "chrMReadsWithClips.sam";

        final File outFile = File.createTempFile("test", ".TargetedMetrics_Coverage");
        outFile.deleteOnExit();

        final String[] args = new String[] {
                "TARGET_INTERVALS=" + singleIntervals,
                "INPUT=" + input,
                "OUTPUT=" + outFile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "LEVEL=ALL_READS",
                "AMPLICON_INTERVALS=" + singleIntervals,
                "SAMPLE_SIZE=" + 0
        };

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<TargetedPcrMetrics, Comparable<Integer>> output = new MetricsFile<>();
        output.read(new FileReader(outFile));

        Assert.assertEquals(output.getMetrics().size(), 1);

        for (final TargetedPcrMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.TOTAL_READS, 2);
        }
        Assert.assertEquals(output.getNumHistograms(), 2);
        final Histogram<Comparable<Integer>> histogram = output.getAllHistograms().get(1);
        Assert.assertTrue(TestNGUtil.compareDoubleWithAccuracy(histogram.getSumOfValues(), 62,0.01));
        Assert.assertTrue(TestNGUtil.compareDoubleWithAccuracy(histogram.get(32).getValue(), 52D, 0.01));
        Assert.assertTrue(TestNGUtil.compareDoubleWithAccuracy(histogram.get(33).getValue(), 10D, 0.01));
    }

    @Test
    public void testCoverageGetTotalOverflow() {
        final Interval interval = new Interval("chr1", 1, 2);
        final TargetMetricsCollector.Coverage coverage = new TargetMetricsCollector.Coverage(interval, 0);
        for (int offset = 0; offset <= interval.length(); offset++) {
            coverage.addBase(offset, Integer.MAX_VALUE - 1);
        }
        Assert.assertEquals((long)coverage.getTotal(), 2 * (long)(Integer.MAX_VALUE - 1));
    }
}
