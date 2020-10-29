/*
 * The MIT License
 *
 * Pierre Lindenbaum - Institut du Thorax - 2016
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
package picard.sam;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.vcf.VcfTestUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.StreamSupport;


public class FilterSamReadsTest extends CommandLineProgramTest {
    @Override
    public String getCommandLineProgramName() {
        return FilterSamReads.class.getSimpleName();
    }

    private static final int READ_LENGTH = 151;
    private final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
    private static final String TEST_DIR = "testdata/picard/sam/";


    @BeforeTest
    public void setUp() {
        builder.setReadLength(READ_LENGTH);

        builder.addPair("mapped_pair_chr1", 0, 1, 151); //should be kept in first test, filtered out in third
        builder.addPair("mapped_pair_chr2", 1, 1, 151); //should be filtered out for first test, and kept in third
        builder.addPair("prove_one_of_pair", 0, 1000, 1000); //neither of these will be kept in any test
        builder.addPair("one_of_pair", 0, 1, 1000); //first read should pass, second should not, but both will be kept in first test
    }

    @DataProvider(name = "dataTestReadFilter")
    public Object[][] dataTestReadFilter() {
        List<String> reads = Arrays.asList(
                "mapped_pair_chr1",
                "prove_one_of_pair",
                "one_of_pair");

        return new Object[][]{
                {FilterSamReads.Filter.includeReadList, reads, 3 * 2},
                {FilterSamReads.Filter.excludeReadList, reads, 1 * 2}
        };
    }

    /**
     * filters a SAM using a reads file
     */
    @Test(dataProvider = "dataTestReadFilter")
    public void testReadFilters(final FilterSamReads.Filter filterType, final List<String> readList, final int expectNumber) throws Exception {

        final File inputSam = File.createTempFile("testSam", ".sam", new File(TEST_DIR));
        inputSam.deleteOnExit();
        final File sortedSamIdx = new File(TEST_DIR, inputSam.getName() + ".idx");
        sortedSamIdx.deleteOnExit();

        try (final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeSAMWriter(builder.getHeader(), false, inputSam)) {

            for (final SAMRecord record : builder) {
                writer.addAlignment(record);
            }
        }

        final File reads = File.createTempFile(TEST_DIR, "reads");
        reads.deleteOnExit();

        try (final FileWriter writer = new FileWriter(reads)) {
            writer.write(String.join("\n", readList));
        } catch (IOException e) {
            e.printStackTrace();
        }
        FilterSamReads filterTest = setupProgram(reads, inputSam, filterType);
        Assert.assertEquals(filterTest.doWork(), 0);

        long count = getReadCount(filterTest);

        Assert.assertEquals(count, expectNumber);
    }

    @DataProvider(name = "badArgumentCombinationsdata")
    public Object[][] badArgumentCombinationsdata() {
        return new Object[][]{
                {FilterSamReads.Filter.includeJavascript, "READ_LIST_FILE"},
                {FilterSamReads.Filter.excludeAligned, "READ_LIST_FILE"},
                {FilterSamReads.Filter.includeAligned, "READ_LIST_FILE"},
                {FilterSamReads.Filter.includePairedIntervals, "READ_LIST_FILE"},
                {FilterSamReads.Filter.includeTagValues, "READ_LIST_FILE"},
                {FilterSamReads.Filter.excludeTagValues, "READ_LIST_FILE"},

                {FilterSamReads.Filter.includeJavascript, "INTERVAL_LIST"},
                {FilterSamReads.Filter.excludeReadList, "INTERVAL_LIST"},
                {FilterSamReads.Filter.includeReadList, "INTERVAL_LIST"},
                {FilterSamReads.Filter.excludeAligned, "INTERVAL_LIST"},
                {FilterSamReads.Filter.includeAligned, "INTERVAL_LIST"},
                {FilterSamReads.Filter.includeTagValues, "INTERVAL_LIST"},
                {FilterSamReads.Filter.excludeTagValues, "INTERVAL_LIST"},

                {FilterSamReads.Filter.excludeReadList, "JAVASCRIPT_FILE"},
                {FilterSamReads.Filter.includeReadList, "JAVASCRIPT_FILE"},
                {FilterSamReads.Filter.excludeAligned, "JAVASCRIPT_FILE"},
                {FilterSamReads.Filter.includeAligned, "JAVASCRIPT_FILE"},
                {FilterSamReads.Filter.includePairedIntervals, "JAVASCRIPT_FILE"},
                {FilterSamReads.Filter.includeTagValues, "JAVASCRIPT_FILE"},
                {FilterSamReads.Filter.excludeTagValues, "JAVASCRIPT_FILE"},

                {FilterSamReads.Filter.excludeReadList, "TAG"},
                {FilterSamReads.Filter.includeReadList, "TAG"},
                {FilterSamReads.Filter.excludeAligned, "TAG"},
                {FilterSamReads.Filter.includeAligned, "TAG"},
                {FilterSamReads.Filter.includePairedIntervals, "TAG"},
                {FilterSamReads.Filter.excludeAligned, "TAG"},
                {FilterSamReads.Filter.includeAligned, "TAG"},

                {FilterSamReads.Filter.excludeReadList, "TAG_VALUE"},
                {FilterSamReads.Filter.includeReadList, "TAG_VALUE"},
                {FilterSamReads.Filter.excludeAligned, "TAG_VALUE"},
                {FilterSamReads.Filter.includeAligned, "TAG_VALUE"},
                {FilterSamReads.Filter.includePairedIntervals, "TAG_VALUE"},
                {FilterSamReads.Filter.excludeAligned, "TAG_VALUE"},
                {FilterSamReads.Filter.includeAligned, "TAG_VALUE"},

        };
    }

    @Test(dataProvider = "badArgumentCombinationsdata")
    public void testBadArgumentCombinations(final FilterSamReads.Filter filter, final String fileArgument) throws IOException {

        final File dummyFile = File.createTempFile(TEST_DIR, "dummy");
        dummyFile.deleteOnExit();
        try (final FileWriter writer = new FileWriter(dummyFile)) {
            writer.write("\n");
        } catch (IOException e) {
            e.printStackTrace();
        }

        final File temp = File.createTempFile("FilterSamReads.output.", ".sam");
        temp.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + TEST_DIR + "aligned.sam",
                String.format("OUTPUT=%s", temp.getAbsolutePath()),
                String.format("FILTER=%s", filter.toString()),
                String.format("%s=%s", fileArgument, dummyFile.getAbsoluteFile()),
        };

        // make sure program invocation failed - inputs validation error
        Assert.assertEquals(runPicardCommandLine(args), 1);
    }

    @DataProvider(name = "dataTestJsFilter")
    public Object[][] dataTestJsFilter() {
        return new Object[][]{
                {TEST_DIR + "aligned.sam", TEST_DIR + "FilterSamReads/filterOddStarts.js", 3},
                {TEST_DIR + "aligned.sam", TEST_DIR + "FilterSamReads/filterReadsWithout5primeSoftClip.js", 0}
        };
    }

    @DataProvider(name = "dataTestPairedIntervalFilter")
    public Object[][] dataTestPairedIntervalFilter() {
        return new Object[][]{
                {"FilterSamReads/filter1.interval_list", 4},
                {"FilterSamReads/filter2.interval_list", 0}
        };
    }

    @DataProvider(name = "dataTestTagFilter")
    public Object[][] dataTestTagFilter() {
        return new Object[][]{
                {TEST_DIR + "aligned.sam", "RG", "0", true, 8},
                {TEST_DIR + "aligned.sam", "RG", "0", false, 0},
                {TEST_DIR + "aligned.sam", "CB", "ACG", false, 3},
                {TEST_DIR + "aligned.sam", "CB", "ACG", true, 5}
        };
    }

    /**
     * filters a SAM using a javascript filter
     */
    @Test(dataProvider = "dataTestJsFilter")
    public void testJavaScriptFilters(final String samFilename, final String javascriptFilename, final int expectNumber) throws Exception {
        // input as SAM file 
        final File inputSam = new File(samFilename);
        final File javascriptFile = new File(javascriptFilename);

        FilterSamReads filterTest = setupProgram(javascriptFile, inputSam, FilterSamReads.Filter.includeJavascript);
        Assert.assertEquals(filterTest.doWork(), 0);

        long count = getReadCount(filterTest);

        Assert.assertEquals(count, expectNumber);
    }

    /**
     * filters a SAM using an interval filter
     */
    @Test(dataProvider = "dataTestPairedIntervalFilter")
    public void testPairedIntervalFilter(final String intervalFilename, final int expectNumber) throws Exception {
        // Build a sam file for testing
        final File inputSam = VcfTestUtils.createTemporaryIndexedFile("testSam", ".bam");

        final SAMFileWriter writer = new SAMFileWriterFactory()
                .setCreateIndex(true).makeSAMWriter(builder.getHeader(), false, inputSam);

        for (final SAMRecord record : builder) {
            writer.addAlignment(record);
        }
        writer.close();

        final File intervalFile = new File(TEST_DIR,intervalFilename);

        FilterSamReads filterTest = setupProgram(intervalFile, inputSam, FilterSamReads.Filter.includePairedIntervals);

        Assert.assertEquals(filterTest.doWork(), 0);

        long count = getReadCount(filterTest);

        Assert.assertEquals(count, expectNumber);
    }

    @DataProvider(name = "dataTestDebugOption")
    public Object[][] dataTestDebugOption() {
        return new Object[][]{
                {null, false},
                {true, true},
                {false, false}
        };
    }

    /**
     * makes sure debug files are created properly
     */
    @Test(dataProvider = "dataTestDebugOption")
    public void testDebugOption(Boolean writeDebugReads, boolean isDebugFileExpected) throws Exception {
        // input as SAM file
        final File inputSam = new File(TEST_DIR + "aligned.sam");
        final File javascriptFile = new File(TEST_DIR + "FilterSamReads/filterOddStarts.js");

        FilterSamReads filterTest = setupProgram(javascriptFile, inputSam, FilterSamReads.Filter.includeJavascript, writeDebugReads);
        Assert.assertEquals(filterTest.doWork(), 0);

        final File inputReadsFile = new File(filterTest.OUTPUT.getParentFile(), IOUtil.basename(filterTest.INPUT) + ".reads");
        Assert.assertEquals(inputReadsFile.exists(), isDebugFileExpected);

        final File outputReadsFile = new File(filterTest.OUTPUT.getParentFile(), IOUtil.basename(filterTest.OUTPUT) + ".reads");
        outputReadsFile.deleteOnExit();
        Assert.assertEquals(outputReadsFile.exists(), isDebugFileExpected);

        // We have to clean up the debug files after each test is run to make sure a clean state is preserved in between tests
        // This mostly affects the input *.reads file because it will always be called "aligned.reads" and will cause future
        // tests to fail if it sticks around and we dont expect it to be written
        Files.deleteIfExists(inputReadsFile.toPath());
    }

    private FilterSamReads setupProgram(final File inputFile, final File inputSam, final FilterSamReads.Filter filter, final Boolean writeDebugReads) throws Exception {
        final FilterSamReads program = new FilterSamReads();
        program.INPUT = inputSam;
        program.OUTPUT = File.createTempFile("FilterSamReads.output.", ".sam");
        program.OUTPUT.deleteOnExit();
        program.FILTER = filter;
        switch (filter) {
            case includePairedIntervals:
                program.INTERVAL_LIST = inputFile;
                break;
            case includeJavascript:
                program.JAVASCRIPT_FILE = inputFile;
                break;
            case includeReadList:
            case excludeReadList:
                program.READ_LIST_FILE = inputFile;
                break;
            default:
                throw new IllegalArgumentException("Not configured for filter=" + filter);
        }
        if (writeDebugReads != null) {
            program.WRITE_READS_FILES = writeDebugReads;
        }

        return program;
    }

    private FilterSamReads setupProgram(final File inputFile, final File inputSam, final FilterSamReads.Filter filter) throws Exception {
        return setupProgram(inputFile, inputSam, filter, null);
    }

       /**
        * filters a SAM using Tag Values
        */
       @Test(dataProvider = "dataTestTagFilter")
       public void testTagFilter(final String samFilename, final String tag, final String tagValue, final boolean includeReads, final int expectNumber) throws Exception {
           // input as SAM file
           final File inputSam = new File(samFilename);
           final FilterSamReads filterTest = new FilterSamReads();
           filterTest.INPUT = inputSam;
           filterTest.OUTPUT = File.createTempFile("FilterSamReads.output.", ".sam");
           filterTest.OUTPUT.deleteOnExit();
           filterTest.FILTER = includeReads ? FilterSamReads.Filter.includeTagValues : FilterSamReads.Filter.excludeTagValues;
           filterTest.TAG = tag;
           filterTest.TAG_VALUE = Arrays.asList(tagValue);
           Assert.assertEquals(filterTest.doWork(),0);
           long count = getReadCount(filterTest);
           Assert.assertEquals(count, expectNumber);
           }

    private long getReadCount(FilterSamReads filterTest) throws Exception {
        final SamReader samReader = SamReaderFactory.makeDefault().open(filterTest.OUTPUT);

        long count = StreamSupport.stream(samReader.spliterator(), false)
                .count();

        samReader.close();
        return count;
    }
}
