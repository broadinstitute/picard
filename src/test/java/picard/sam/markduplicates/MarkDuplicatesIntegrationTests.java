package picard.sam.markduplicates;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.Optional;
import java.util.stream.Collectors;

public class MarkDuplicatesIntegrationTests extends CommandLineProgramTest {

    protected static String TEST_BASE_NAME = null;
    protected static File TEST_DATA_DIR = null;

    @BeforeClass
    public void setUp() {
        TEST_BASE_NAME = "MD_IT";
        TEST_DATA_DIR = new File("testdata/picard/sam/MarkDuplicates/IntTest");
    }

    @DataProvider(name = "MarkDuplicatesIntegrationTest")
    public Object[][] createMarkDuplicatesIntegrationTest() {
        return new Object[][]{
                {"duplicates_test.sam"},
                {"duplicates_big_test.sam"},
                {"optical_duplicates_test.sam"},
        };
    }

    @DataProvider(name = "SecondarySupplementaryUnmappedTest")
    public Object[][] createSecondarySupplementaryUnmappedTest() {
        return new Object[][]{
                {"secondary_supplementary_test.sam"},
                {"secondary_supplementary_small_test.sam"},
                {"secondary_supplementary_unmapped_big_test.bam"}
        };
    }

    @Test(dataProvider = "MarkDuplicatesIntegrationTest")
    public void markDuplicatesWithoutRemovingDuplicatesIntegrationTest(final String input) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final ArrayList<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR,input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineDuplicationMetrics(input, output, metrics, false);
    }

    @Test(dataProvider = "MarkDuplicatesIntegrationTest")
    public void markDuplicatesWithRemovingDuplicatesIntegrationTest(final String input) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR,input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());
        args.add("REMOVE_DUPLICATES=true");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineDuplicationMetrics(input, output, metrics, true);
    }


    @Test(dataProvider = "SecondarySupplementaryUnmappedTest")
    public void secondarySupplementaryUnmappedWithRemovingDuplicatesIntegrationTest(final String input) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR,input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());
        args.add("REMOVE_DUPLICATES=true");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineDuplicationMetrics(input, output, metrics, true);
    }


    @Test(dataProvider = "SecondarySupplementaryUnmappedTest")
    public void secondarySupplementaryUnmappedWithoutRemovingDuplicatesIntegrationTest(final String input) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR,input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineDuplicationMetrics(input, output, metrics, false);
    }

    @Test
    public void BarcodesWithoutRemovingDuplicatesIntegrationTest() {
        final String input = "barcodes_test.sam";
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR,input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineDuplicationMetrics(input, output, metrics, false);
    }


    @Test
    public void BarcodesWithRemovingDuplicatesIntegrationTest() {
        final String input = "barcodes_test.sam";
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + new File(TEST_DATA_DIR,input).getAbsolutePath());
        final File output = new File(outputDir, TEST_BASE_NAME + ".sam");
        args.add("OUTPUT=" + output.getAbsolutePath());
        final File metrics = new File(outputDir, TEST_BASE_NAME + ".integration_metrics");
        args.add("METRICS_FILE=" + metrics.getAbsolutePath());
        args.add("REMOVE_DUPLICATES=true");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineDuplicationMetrics(input, output, metrics, true);
    }

    private void examineDuplicationMetrics(String input, File output, File metrics, boolean removeDuplicates) {
        final List<DuplicationMetrics> metricList = MetricsFile.readBeans(metrics);

        DuplicationMetrics metric = metricList.get(0);
        for(DuplicationMetrics currentMetric:metricList) {
            if (currentMetric.LIBRARY.equals("Unknown Library")) {
                metric = currentMetric;
            }
        }

        switch (input) {
            case "duplicates_test.sam" : examineDuplicatesFile(metric, output, removeDuplicates);
                break;
            case "duplicates_big_test.sam" : examineDuplicatesBigFile(metric, output, removeDuplicates);
                break;
            case "optical_duplicates_test.sam" : examineOpticalDuplicatesFile(metric, output, removeDuplicates);
                break;
            case "barcodes_test.sam" : examineBarcodesFile(metric, output, removeDuplicates);
                break;
            case "secondary_supplementary_test.sam" : examineSecondarySupplementaryFile(metric, output, removeDuplicates);
                break;
            case "secondary_supplementary_small_test.sam" : examineSecondarySupplementarySmallFile(metric, output, removeDuplicates);
                break;
            case "secondary_supplementary_unmapped_big_test.bam" : examineSecondarySupplementaryUnmappedBigFile(metric, output, removeDuplicates);
                break;
            default: // if file from data provider is not in case, it must fall
                Assert.assertEquals(true,false);
        }
    }

    private void examineDuplicatesFile(DuplicationMetrics metric, File output, boolean removeDuplicates) {
        HashMap<String, Boolean> duplicateFlagsMap = getSamRecordBooleanHashMap(output);
        if (removeDuplicates) {
            Assert.assertEquals(duplicateFlagsMap.getOrDefault("read_2_too_many_gaps", false).booleanValue(), false); ///////
        }
        else {
            Assert.assertEquals(duplicateFlagsMap.getOrDefault("read_2_too_many_gaps", false).booleanValue(), true);
        }

        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 1);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 3);
        Assert.assertEquals(metric.UNMAPPED_READS, 1);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 1);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 1);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 0);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.428571);
    }

    private void examineDuplicatesBigFile(final DuplicationMetrics metric, final File output, boolean removeDuplicates) {
        if (removeDuplicates) {
            Assert.assertEquals(getSamRecordsList(output), getSamRecordsList(new File(TEST_DATA_DIR, "expected_output_after_removing_duplicates" + ".sam")));
        } else {
            Assert.assertEquals(getSamRecordBooleanHashMap(output), getSamRecordBooleanHashMap(new File(TEST_DATA_DIR, "expected_dupes_flags" + ".sam")));
        }
        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 439);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 500);
        Assert.assertEquals(metric.UNMAPPED_READS, 0);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 0);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.878);
    }

    private void examineOpticalDuplicatesFile(DuplicationMetrics metric, File output, boolean removeDuplicates) {
        HashMap<String, Boolean> duplicateFlagsMap = getSamRecordBooleanHashMap(output);
        if (removeDuplicates) {
            Assert.assertEquals(duplicateFlagsMap.getOrDefault("C4N4WACXX140821:8:1112:2344:1985", false).booleanValue(), false);
        }
        else {
            Assert.assertEquals(duplicateFlagsMap.getOrDefault("C4N4WACXX140821:8:1112:2344:1985", false).booleanValue(), true);
        }

        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 1);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 1);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 2);
        Assert.assertEquals(metric.UNMAPPED_READS, 0);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 0);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.5);
    }

    private void examineSecondarySupplementaryUnmappedBigFile(DuplicationMetrics metric, File output, boolean removeDuplicates) {
        if (removeDuplicates) {
            Assert.assertEquals(getSamRecordsNameList(output), getSamRecordsNameList(new File(TEST_DATA_DIR, "expected_sec_sup_rd_true.bam")));
        } else {
            Assert.assertEquals(getSamRecordsNameList(output), getSamRecordsNameList(new File(TEST_DATA_DIR, "expected_sec_sup_rd_false.bam")));
        }

        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 18254);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 9601);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 136562);
        Assert.assertEquals(metric.UNMAPPED_READS, 1211);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 255);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 1211);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 1424);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.134008);
    }

    private void examineSecondarySupplementarySmallFile(DuplicationMetrics metric, File output, boolean removeDuplicates) {
        HashMap<String, Boolean> duplicateFlagsMap = getSamRecordBooleanHashMap(output);
        if (removeDuplicates) {
            Assert.assertEquals(duplicateFlagsMap.getOrDefault("ST-E00297:149016593:H3GVWCCXX:5:2214:10145:57038", false).booleanValue(), false);
            Assert.assertEquals(getSamRecordsList(output).stream()
                                                         .filter(rec -> rec.getReadName().equals("ST-E00297:149016593:H3GVWCCXX:5:2214:10145:57038"))
                                                         .count(),2);
            Assert.assertEquals(getSamRecordsList(output).stream()
                                                         .filter(rec -> rec.getReadName().equals("ST-E00297:149016593:H3GVWCCXX:3:1218:20812:27591"))
                                                         .count(),2);
        }
        else {
            Assert.assertEquals(duplicateFlagsMap.getOrDefault("ST-E00297:149016593:H3GVWCCXX:5:2214:10145:57038", false).booleanValue(), false);
            Assert.assertEquals(getSamRecordsList(output).stream()
                                                         .filter(rec -> rec.getReadName().equals("ST-E00297:149016593:H3GVWCCXX:5:2214:10145:57038"))
                                                         .count(),2);
            Assert.assertEquals(getSamRecordsList(output).stream()
                                                         .filter(rec -> rec.getReadName().equals("ST-E00297:149016593:H3GVWCCXX:3:1218:20812:27591"))
                                                         .count(),2);
        }

        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 0);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 1);
        Assert.assertEquals(metric.UNMAPPED_READS, 0);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 1);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.0);
    }

    private void examineSecondarySupplementaryFile(DuplicationMetrics metric, File output, boolean removeDuplicates) {
        if (removeDuplicates) {
            Assert.assertEquals(getSamRecordsList(output).stream()
                                                         .filter(rec -> rec.getReadName().equals("ST-E00297:149016593:H3GVWCCXX:5:2214:10145:57038"))
                                                         .count(),0);
        }
        else {
            Assert.assertEquals(getSamRecordsList(output).stream()
                                                         .filter(rec -> rec.getReadName().equals("ST-E00297:149016593:H3GVWCCXX:5:2214:10145:57038"))
                                                         .count(),102);
        }

        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 1);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 2);
        Assert.assertEquals(metric.UNMAPPED_READS, 0);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 100);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.5);
    }

    private void examineBarcodesFile(DuplicationMetrics metric, File output, boolean removeDuplicates) {
        if (removeDuplicates) {
            Assert.assertEquals(getBarcodeFlagsHashMap(output), getBarcodeFlagsHashMap(new File("testdata/picard/sam/MarkDuplicates/IntTest/expected_barcodes_with_removing_duplicates.sam")));
        } else {
            Assert.assertEquals(getBarcodeFlagsHashMap(output), getBarcodeFlagsHashMap(new File("testdata/picard/sam/MarkDuplicates/IntTest/expected_barcodes_without_removing_duplicates.sam")));
        }

        Assert.assertEquals(metric.READ_PAIR_DUPLICATES, 2);
        Assert.assertEquals(metric.READ_PAIR_OPTICAL_DUPLICATES, 0);
        Assert.assertEquals(metric.READ_PAIRS_EXAMINED, 3);
        Assert.assertEquals(metric.UNMAPPED_READS, 0);
        Assert.assertEquals(metric.UNPAIRED_READ_DUPLICATES, 0);
        Assert.assertEquals(metric.UNPAIRED_READS_EXAMINED, 0);
        Assert.assertEquals(metric.SECONDARY_OR_SUPPLEMENTARY_RDS, 0);
        Assert.assertEquals(metric.PERCENT_DUPLICATION, 0.666667);

    }

    private HashMap<String, String> getBarcodeFlagsHashMap(File output) {
        SamReader reader = SamReaderFactory.make().open(output);
        HashMap<String, String> barcodeFlagsMap = new HashMap<>();
        for(SAMRecord record : reader) {
            barcodeFlagsMap.put(record.getReadName(), getBarcodeFlag(record));
        }
        return barcodeFlagsMap;
    }

    private String getBarcodeFlag(SAMRecord record) {
        if (EstimateLibraryComplexity.getReadBarcodeValue(record, "BC") != 0 ) {
            return "BARCODE_TAG";
        }
        if (EstimateLibraryComplexity.getReadBarcodeValue(record, "RX") != 0) {
            return "READ_ONE_BARCODE_TAG";
        }
        if (EstimateLibraryComplexity.getReadBarcodeValue(record, "RX") != 0) {
            return "READ_TWO_BARCODE_TAG";
        }
        return null;
    }

    private HashMap<String, Boolean> getSamRecordBooleanHashMap(File output) {
        SamReader reader = SamReaderFactory.make().open(output);
        HashMap<String,Boolean> duplicateFlagsMap = new HashMap<>();
        for(SAMRecord record : reader) {
            Set<SAMFlag> flags = record.getSAMFlags();
            Optional<SAMFlag> duplicateFlag = flags.stream().filter(flag -> flag.intValue()==1024).findFirst();
            duplicateFlagsMap.put(record.getReadName(), duplicateFlag.isPresent());
        }
        return duplicateFlagsMap;
    }


    private ArrayList<SAMRecord> getSamRecordsList(File output) {
        SamReader reader = SamReaderFactory.make().open(output);
        ArrayList <SAMRecord> SamRecordList = new ArrayList<>();
        for(SAMRecord record : reader) {
            SamRecordList.add(record);
        }
        return SamRecordList;
    }

    private ArrayList<String> getSamRecordsNameList(File output) {
        SamReader reader = SamReaderFactory.make().open(output);
        ArrayList <String> SamRecordList = new ArrayList<>();
        for(SAMRecord record : reader) {
            SamRecordList.add(record.getReadName());
        }
        return SamRecordList;
    }
    @Override
    public String getCommandLineProgramName() { return MarkDuplicates.class.getSimpleName();
    }
}

