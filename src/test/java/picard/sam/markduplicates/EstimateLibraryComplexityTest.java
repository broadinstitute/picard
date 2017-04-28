/*
 * The MIT License
 *
 * Copyright (c) 2015 Nils Homer
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

package picard.sam.markduplicates;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class EstimateLibraryComplexityTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/EstimateLibraryComplexity");

    @DataProvider(name = "testSimpleDuplicate")
    public Object[][] createDataTestSimpleDuplicate() {
        return new Object[][]{
                {"dupes.sam", 2, 2},
                {"big_dupes.sam", 12, 500}
        };
    }

    @DataProvider(name = "testMaxDiffRate")
    public Object[][] createDataTestMaxDiffRate() {
        return new Object[][]{
                {"dupes.sam", 0, 2},
                {"big_dupes.sam", 8, 500}
        };
    }

    @DataProvider(name = "testSimpleDuplicateWithMaxReadLength")
    public Object[][] createDataTestSimpleDuplicateWithMaxReadLength() {
        return new Object[][]{
                {"dupes.sam", 2, 2},
                {"big_dupes.sam", 512, 500}
        };
    }

    @DataProvider(name = "testDefaultMinGroupCount")
    public Object[][] createDataTestDefaultMinGroupCount() {
        return new Object[][]{
                {"dupes.sam", 0, 0},
                {"big_dupes.sam", 8, 497}
        };
    }

    @DataProvider(name = "testSimpleDuplicatesWithSecondaryAndSupplementaryRecords")
    public Object[][] createDataTestSimpleDuplicatesWithSecondaryAndSupplementaryRecords() {
        return new Object[][]{
                {"dupes_with_sos.sam", 2, 2},
                {"big_dupes_with_sos.sam", 12, 500}
        };
    }

    public String getCommandLineProgramName() {
        return EstimateLibraryComplexity.class.getSimpleName();
    }

    private void examineMetricsFile(final File output, final int numDuplicates, final int numReadPairsExamined) {
        final List<DuplicationMetrics> metricsList = MetricsFile.readBeans(output);
        Assert.assertEquals(metricsList.size(), 1);
        final DuplicationMetrics metrics = metricsList.get(0);
        Assert.assertEquals(metrics.READ_PAIR_DUPLICATES * 2 + metrics.UNPAIRED_READ_DUPLICATES, numDuplicates);
        Assert.assertEquals(metrics.READ_PAIRS_EXAMINED, numReadPairsExamined);
    }

    /**
     * Finds duplicates as expected.
     */
    @Test(dataProvider = "testSimpleDuplicate")
    public void testSimpleDuplicate(final String testName, final int numDuplicates, final int numReadPairsExamined) throws IOException {
        final File input = new File(TEST_DATA_DIR, testName);
        final File output = File.createTempFile("estimateLibraryComplexity", ".els_metrics");
        output.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("MIN_GROUP_COUNT=1");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineMetricsFile(output, numDuplicates, numReadPairsExamined);
    }

    /**
     * Finds duplicates as expected ignoring secondary and supplementary records.
     */
    @Test(dataProvider = "testSimpleDuplicatesWithSecondaryAndSupplementaryRecords")
    public void testSimpleDuplicatesWithSecondaryAndSupplementaryRecords(final String testName, final int numDuplicates, final int numReadPairsExamined) throws IOException {
        final File input = new File(TEST_DATA_DIR, testName);
        final File output = File.createTempFile("estimateLibraryComplexity", ".els_metrics");
        output.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("MIN_GROUP_COUNT=1");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineMetricsFile(output, numDuplicates, numReadPairsExamined);
    }

    /**
     * Does not find duplicates since the difference rate was too high across the entire read
     */
    @Test(dataProvider = "testMaxDiffRate")
    public void testMaxDiffRate(final String testName, final int numDuplicates, final int numReadPairsExamined) throws IOException {
        final File input = new File(TEST_DATA_DIR, testName);
        final File output = File.createTempFile("estimateLibraryComplexity", ".els_metrics");
        output.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("MAX_DIFF_RATE=0.0");
        args.add("MIN_GROUP_COUNT=1");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineMetricsFile(output, numDuplicates, numReadPairsExamined);
    }

    /**
     * Finds duplicates since the we examine only the fist ten bases.
     */
    @Test(dataProvider = "testSimpleDuplicateWithMaxReadLength")
    public void testSimpleDuplicateWithMaxReadLength(final String testName, final int numDuplicates, final int numReadPairsExamined) throws IOException {
        final File input = new File(TEST_DATA_DIR, testName);
        final File output = File.createTempFile("estimateLibraryComplexity", ".els_metrics");
        output.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("OUTPUT=" + output.getAbsolutePath());
        args.add("MAX_DIFF_RATE=0.0");
        args.add("MIN_GROUP_COUNT=1");
        args.add("MAX_READ_LENGTH=10");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineMetricsFile(output, numDuplicates, numReadPairsExamined);
    }

    /**
     * Does not find any duplicates since there was only one group of duplicates of size one.  Also
     * there are no reads examined due to this filtering step.
     */
    @Test(dataProvider = "testDefaultMinGroupCount")
    public void testDefaultMinGroupCount(final String testName, final int numDuplicates, final int numReadPairsExamined) throws IOException {
        final File input = new File(TEST_DATA_DIR, testName);
        final File output = File.createTempFile("estimateLibraryComplexity", ".els_metrics");
        output.deleteOnExit();

        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("OUTPUT=" + output.getAbsolutePath());

        Assert.assertEquals(runPicardCommandLine(args), 0);
        examineMetricsFile(output, numDuplicates, numReadPairsExamined); // no read pairs examined!!!
    }
}
