package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.BufferedLineReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.analysis.CollectQualityYieldMetrics;
import picard.cmdline.CommandLineProgramTest;
import picard.nio.GATKBucketUtils;
import picard.nio.PicardHtsPath;
import picard.sam.util.SamTestUtil;
import htsjdk.samtools.DownsamplingIteratorFactory.Strategy;
import picard.util.GCloudTestUtils;
import picard.util.TestNGUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Random;

import static htsjdk.samtools.DownsamplingIteratorFactory.Strategy.Chained;
import static htsjdk.samtools.DownsamplingIteratorFactory.Strategy.ConstantMemory;
import static htsjdk.samtools.DownsamplingIteratorFactory.Strategy.HighAccuracy;


/**
 * Created by farjoun on 9/21/17.
 */
public class DownsampleSamTest extends CommandLineProgramTest {
    private static final String sample = "TestSample";
    private static final String readGroupId = "TestReadGroup";
    private static final String platform = "ILLUMINA";
    private static final String library = "TestLibrary";

    private static final File TEST_DIR = new File("testdata/picard/sam/DownsampleSam");
    private static final File dict = new File(TEST_DIR, "header.dict");

    private File tempDir;
    private File tempSamFile;

    private SAMRecordSetBuilder setBuilder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.coordinate);
    private SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(readGroupId);

    @Override
    public String getCommandLineProgramName() {
            return DownsampleSam.class.getSimpleName();
    }

    //create a samfile with one tile and randomly placed reads within it.
    @BeforeTest
    void setupBuilder() throws IOException {
        final int numReads = 10000;
        final String flowCellBarcode = "TESTBARCODE";


        final String separator = ":";
        final int lane = 1;
        final int tile = 2203;

        //use fixed random seed to eliminate possibility of intermittent failures
        final Random rg = new Random(31);
        setBuilder.setReadGroup(readGroupRecord);
        setBuilder.setUseNmFlag(true);

        for (int i = 0; i < numReads; i++) {

            final String readName = flowCellBarcode + separator + lane + separator + tile + separator + i;
            setBuilder.addPair(readName, 1, 1, 100);
        }

        tempDir = IOUtil.createTempDir("ds_test_Downsampling").toFile();
        tempSamFile = File.createTempFile("DownsampleSam", ".bam", tempDir);

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

    @DataProvider(name = "ValidArgumentsTestProvider")
    public Object[][] ValidArgumentsTestProvider() {
        final List<Object[]> objects = new ArrayList<>();
        for(final Strategy strategy : Strategy.values())
            for( final Integer seed : new Integer[]{1, null})
        for (double i = 0.3; i <= 1; i += .1) {
            final Object[] array = {i, strategy, seed};
            objects.add(array);
        }
        return objects.toArray(new Object[1][]);
    }

    // test removing some reads from a sparse, single tile bam
    @Test(dataProvider = "ValidArgumentsTestProvider")
    public void testDownsampleStrategies(final double fraction, final Strategy strategy, final Integer seed) throws IOException {
            testDownsampleWorker(new PicardHtsPath(tempSamFile), fraction, strategy.name(), seed);
    }

    private PicardHtsPath testDownsampleWorker(final PicardHtsPath samFile,  final double fraction, final String strategy, final Integer seed) throws IOException {
        return testDownsampleWorker(samFile, fraction, strategy, seed, Optional.empty(), Optional.empty(), Optional.empty());
    }

        // tsato: add 'final PicardHtsPath outputFile' as argument
    private PicardHtsPath testDownsampleWorker(final PicardHtsPath samFile,  final double fraction, final String strategy, final Integer seed,
                                               final Optional<PicardHtsPath> outputFile, final Optional<PicardHtsPath> metricsFile,
                                               final Optional<PicardHtsPath> referenceFile) throws IOException {
        final PicardHtsPath downsampled = outputFile.isEmpty() ? new PicardHtsPath(File.createTempFile("DownsampleSam", ".bam", tempDir)) : outputFile.get();

        final List<String> args = new ArrayList<>(Arrays.asList(
                "INPUT=" + samFile.getURIString(),
                "OUTPUT=" + downsampled.getURIString(),
                "PROBABILITY=" + fraction,
                "STRATEGY=" + strategy,
                "RANDOM_SEED=" + ((seed==null)?"null":seed.toString()),
                "CREATE_INDEX=true"));
        // consider detecting leading/trailing spaces...

        if (referenceFile.isPresent()){
            args.add("REFERENCE_SEQUENCE=" + referenceFile.get().getURIString());
        }

        if (metricsFile.isPresent()){
            args.add("METRICS_FILE=" + metricsFile.get().getURIString());
        }

        // make sure results is successful
        Assert.assertEquals(runPicardCommandLine(args), 0);

        // make sure that the resulting BAM is valid.
        final ValidateSamFile validateSamFile = new ValidateSamFile();

        validateSamFile.INPUT = downsampled;
        Assert.assertEquals(validateSamFile.doWork(), 0);

        // make sure that the total number of record in the resulting file in in the ballpark:
        // don't run this when the seed is null since that is non-deterministic and might (unlikely) fail to hit the bounds.
        if (seed!=null) {
            TestNGUtil.assertGreaterThan(SamTestUtil.countSamTotalRecord(downsampled.toPath()), fraction * .8 * SamTestUtil.countSamTotalRecord(samFile.toPath()));
            TestNGUtil.assertLessThan(SamTestUtil.countSamTotalRecord(downsampled.toPath()), fraction * 1.2 * SamTestUtil.countSamTotalRecord(samFile.toPath()));
        }

        if (metricsFile.isPresent()){
            final MetricsFile<CollectQualityYieldMetrics.QualityYieldMetrics, Integer> metrics = new MetricsFile<>();
            final Reader metricReader = IOUtil.openFileForBufferedReading(metricsFile.get().toPath());
            metrics.read(metricReader);
            CloserUtil.close(metricReader);
            Assert.assertEquals(metrics.getMetrics().get(0).READ_LENGTH, 101);
        }

        return downsampled;
    }

    @DataProvider(name = "RepeatedDownsamplingProvider")
    public Object[][] repeatedDownsamplingProvider() {
        final List<Object[]> rets = new ArrayList<>();
        rets.add(new Object[]{Arrays.asList(DownsamplingIteratorFactory.Strategy.ConstantMemory, ConstantMemory), Arrays.asList(2,1)});
        rets.add(new Object[]{Arrays.asList(DownsamplingIteratorFactory.Strategy.ConstantMemory, ConstantMemory), Arrays.asList(0,0)});
        rets.add(new Object[]{Arrays.asList(DownsamplingIteratorFactory.Strategy.ConstantMemory, ConstantMemory), Arrays.asList(Integer.MAX_VALUE,Integer.MAX_VALUE)});
        rets.add(new Object[]{Arrays.asList(DownsamplingIteratorFactory.Strategy.ConstantMemory, ConstantMemory), Arrays.asList(Integer.MIN_VALUE,Integer.MIN_VALUE)});
        rets.add(new Object[]{Arrays.asList(DownsamplingIteratorFactory.Strategy.ConstantMemory, ConstantMemory), Arrays.asList(Integer.MIN_VALUE,Integer.MAX_VALUE)});
        rets.add(new Object[]{Arrays.asList(DownsamplingIteratorFactory.Strategy.ConstantMemory, ConstantMemory), Arrays.asList(Integer.MAX_VALUE,0)});
        rets.add(new Object[]{Arrays.asList(ConstantMemory, ConstantMemory), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(Chained, ConstantMemory), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(Chained, ConstantMemory), Arrays.asList(1,3)});
        rets.add(new Object[]{Arrays.asList(ConstantMemory, Chained), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(HighAccuracy, ConstantMemory), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(ConstantMemory, HighAccuracy), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(Chained, Chained), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(HighAccuracy, Chained), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(HighAccuracy, HighAccuracy), Arrays.asList(1,1)});
        rets.add(new Object[]{Arrays.asList(Chained, HighAccuracy), Arrays.asList(1,1)});

        //randomly generate some sequences to test out
        final Strategy[] availableStratagies = Strategy.values();
        final Random random = new Random(12345);
        for (int i=0; i<20; i++) {
            final List<Strategy> strategies = new ArrayList<>();
            final List<Integer> seeds = new ArrayList<>();

            while (strategies.size() < 5) {
                final int seed = random.nextInt(3);
                final Strategy strategy = availableStratagies[random.nextInt(availableStratagies.length)];

                seeds.add(seed);
                strategies.add(strategy);
            }
            rets.add(new Object[]{strategies, seeds});
        }

        return rets.toArray(new Object[0][]);
    }

    @Test(dataProvider = "RepeatedDownsamplingProvider")
    public void testRepeatedDownsampling(List<Strategy> strategies, List<Integer> seeds) throws IOException {
        PicardHtsPath input = new PicardHtsPath(tempSamFile);
        final long nReadsOriginal = SamTestUtil.countSamTotalRecord(input.toPath());
        double totalFraction = 1;
        for (int i = 0 ; i < strategies.size(); i++) {
            // The downsampled file will be the input for the next iteration
            input = testDownsampleWorker(input, 0.5, strategies.get(i).toString(), seeds.get(i));
            totalFraction *= 0.5;

            final long nReadsNow = SamTestUtil.countSamTotalRecord(input.toPath());
            Assert.assertTrue(nReadsNow > 0.8 * totalFraction * nReadsOriginal);
            Assert.assertTrue(nReadsNow < 1.2 * totalFraction * nReadsOriginal);
        }
    }

    private static final PicardHtsPath NA12878_MINI = new PicardHtsPath(GCloudTestUtils.getTestInputPath() + "picard/bam/CEUTrio.HiSeq.WGS.b37.NA12878.20.21_n100.bam");
    private static final PicardHtsPath NA12878_MINI_CRAM = new PicardHtsPath(GCloudTestUtils.getTestInputPath() + "picard/bam/CEUTrio.HiSeq.WGS.b37.NA12878.20.21_n100.cram");
    final int DEFAULT_RANDOM_SEED = 42;
    private static final PicardHtsPath HG19_REFERENCE = new PicardHtsPath("gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta");


    @Test(groups = "cloud")
    public void testCloud() throws IOException {
        // Local test files for testing, to delete before merge
//        final String bamOutput = GATKBucketUtils.getTempFilePath("downsample", "bam");
//        final String metricsOutput = GATKBucketUtils.getTempFilePath("metrics", "txt");

        // Test bam (input/output)
        // Local output
        testDownsampleWorker(NA12878_MINI, 0.5, ConstantMemory.toString(), DEFAULT_RANDOM_SEED);
        // Output is in the cloud
        final String bamOutput1 = GATKBucketUtils.getTempFilePath(GCloudTestUtils.TEST_OUTPUT_DEFAULT + "downsample", "bam");
        testDownsampleWorker(NA12878_MINI, 0.5, ConstantMemory.toString(), DEFAULT_RANDOM_SEED,
                Optional.of(new PicardHtsPath(bamOutput1)), Optional.empty(), Optional.empty());

        // Test generating a metrics file in the cloud
        final String metricsOutput = GATKBucketUtils.getTempFilePath(GCloudTestUtils.TEST_OUTPUT_DEFAULT + "metrics", "txt");
        final String bamOutput2 = GATKBucketUtils.getTempFilePath(GCloudTestUtils.TEST_OUTPUT_DEFAULT + "downsample", "bam");
        testDownsampleWorker(NA12878_MINI, 0.5, ConstantMemory.toString(), DEFAULT_RANDOM_SEED,
                Optional.of(new PicardHtsPath(bamOutput2)), Optional.of(new PicardHtsPath(metricsOutput)), Optional.empty());

        // Test cram (input/output). Temporarily turned off as we investigate why it hangs
        boolean testCram = false;
        if (testCram){
            final String cramOutputInCloud = GATKBucketUtils.getTempFilePath("downsample", "cram");
            testDownsampleWorker(NA12878_MINI_CRAM, 0.5, ConstantMemory.toString(), 42, Optional.of(new PicardHtsPath(cramOutputInCloud)), Optional.empty(), Optional.of(HG19_REFERENCE));
        }


    }
}