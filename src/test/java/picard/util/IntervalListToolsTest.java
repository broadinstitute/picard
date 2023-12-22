/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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
package picard.util;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;
import picard.nio.PicardBucketUtils;
import picard.nio.PicardIOUtils;
import picard.nio.PicardHtsPath;
import picard.util.IntervalList.IntervalListScatterMode;
import picard.util.IntervalList.IntervalListScatterer;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.Scanner;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by farjoun on 10/22/17.
 */

public class IntervalListToolsTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIR = "testdata/picard/util/";
    private static final String CLOUD_DATA_DIR = GCloudTestUtils.getTestInputPath() + "picard/intervals/";
    private static final String CLOUD_OUTPUT_DIR = GCloudTestUtils.getTestStaging() + "picard/";
    private final Path scatterable = Paths.get(TEST_DATA_DIR, "scatterable.interval_list");
    private final PicardHtsPath scatterableCloud = new PicardHtsPath(CLOUD_DATA_DIR + "scatterable.interval_list");
    private final Path scatterableStdin = Paths.get(TEST_DATA_DIR, "scatterable_stdin");
    private final Path secondInput = Paths.get(TEST_DATA_DIR, "secondInput.interval_list");
    private final PicardHtsPath secondInputCloud = new PicardHtsPath(CLOUD_DATA_DIR + "secondInput.interval_list");
    private final Path largeScatterable = Paths.get(TEST_DATA_DIR, "large_scatterable.interval_list");
    private final Path abutting = Paths.get(TEST_DATA_DIR, "abutting.interval_list");
    private final Path abuttingCombined = Paths.get(TEST_DATA_DIR, "abutting_combined.interval_list");
    private final Path abuttingNotCombined = Paths.get(TEST_DATA_DIR, "abutting_notcombined.interval_list");
    private static final List<File> LARGER_EXPECTED_WITH_REMAINDER_FILES = Arrays.asList(new File(TEST_DATA_DIR, "largeScattersWithRemainder").listFiles());
    private static final List<IntervalList> LARGER_EXPECTED_WITH_REMAINDER_LISTS = LARGER_EXPECTED_WITH_REMAINDER_FILES.stream()
            .sorted().flatMap(l -> Arrays.asList(l.listFiles())
                    .stream().map(f -> IntervalList.fromFile(f))).collect(Collectors.toList());
    private static final boolean DEFAULT_INVERT_ARG = false;
    private static final boolean DEFAULT_UNIQUE_ARG = false;
    private static final IntervalListTools.Action DEFAULT_ACTION = IntervalListTools.Action.INTERSECT;
    private static final boolean DEFAULT_DONT_MERGE_ABUTTING = false;
    private static final String INTERVAL_LIST_EXTENSION = ".interval_list"; // Note that '.' is included

    @Test
    public void testSecondInputValidation() {
        IntervalListTools intervalListTools = new IntervalListTools();
        String[] errors = intervalListTools.customCommandLineValidation();
        Assert.assertNull(errors);

        intervalListTools.SECOND_INPUT = new ArrayList<>();
        errors = intervalListTools.customCommandLineValidation();
        Assert.assertNull(errors);

        intervalListTools.SECOND_INPUT.add(new PicardHtsPath("fakefile"));
        errors = intervalListTools.customCommandLineValidation();
        Assert.assertEquals(errors.length, 1);
    }

    @Test()
    public void testCountOutputValidation() {
        final IntervalListTools intervalListTools = new IntervalListTools();

        for (IntervalListTools.Output output_value : IntervalListTools.Output.values()) {
            intervalListTools.OUTPUT_VALUE = output_value;
            intervalListTools.COUNT_OUTPUT = null;
            String[] errors = intervalListTools.customCommandLineValidation();
            Assert.assertNull(errors);

            intervalListTools.COUNT_OUTPUT = new PicardHtsPath("fakefile");
            errors = intervalListTools.customCommandLineValidation();
            if (output_value == IntervalListTools.Output.NONE) {
                Assert.assertEquals(errors.length, 1);
            } else {
                Assert.assertNull(errors);
            }
        }
    }

    @Override
    public String getCommandLineProgramName() {
        return IntervalListTools.class.getSimpleName();
    }

    @DataProvider
    public Iterator<Object[]> ActionsTest() {
        final Stream stream1 = Arrays.stream(IntervalListTools.Action.values()).map(a -> new Object[]{a, scatterable});
        final Stream stream2 = Arrays.stream(IntervalListTools.Action.values()).map(a -> new Object[]{a, scatterableStdin});
        return Stream.concat(stream1, stream2).iterator();
    }

    // test that all actions work. but not test output at all.
    @Test(dataProvider = "ActionsTest")
    public void testAllActions(final IntervalListTools.Action action, final Path path) throws IOException {
        final File ilOut = File.createTempFile("IntervalListTools", "interval_list");
        ilOut.deleteOnExit();

        final List<String> args = new ArrayList<>();

        args.add("ACTION=" + action.toString());
        args.add("INPUT=" + path.toUri());

        if (action.takesSecondInput) {
            args.add("SECOND_INPUT=" + secondInput);
        }
        args.add("OUTPUT=" + ilOut);

        Assert.assertEquals(runPicardCommandLine(args), 0);
    }

    @DataProvider
    public Object[][] actionAndTotalBasesData() {
        return new Object[][]{
                {IntervalListTools.Action.CONCAT, 341, 4},
                {IntervalListTools.Action.UNION, 201, 2},
                {IntervalListTools.Action.INTERSECT, 140, 2},
                {IntervalListTools.Action.SUBTRACT, 60, 2},
                {IntervalListTools.Action.SYMDIFF, 61, 3},
                {IntervalListTools.Action.OVERLAPS, 150, 2}
        };
    }

    @Test(dataProvider = "actionAndTotalBasesData")
    public void testActions(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final IntervalList il = tester(action);
        Assert.assertEquals(il.getBaseCount(), bases, "unexpected number of bases found.");
        Assert.assertEquals(il.getIntervals().size(), intervals, "unexpected number of intervals found.");

        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.BASES), bases, "unexpected number of bases written to count_output file.");
        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.INTERVALS), intervals, "unexpected number of intervals written to count_output file.");
    }

    @DataProvider
    public Object[][] actionAndTotalBasesWithInvertData() {
        final long totalBasesInDict = IntervalList.fromPath(secondInput).getHeader().getSequenceDictionary().getReferenceLength();
        final int totalContigsInDict = IntervalList.fromPath(secondInput).getHeader().getSequenceDictionary().size();
        return new Object[][]{
                {IntervalListTools.Action.CONCAT, totalBasesInDict - 201, 2 + totalContigsInDict},
                {IntervalListTools.Action.UNION, totalBasesInDict - 201, 2 + totalContigsInDict},
                {IntervalListTools.Action.INTERSECT, totalBasesInDict - 140, 2 + totalContigsInDict},
                {IntervalListTools.Action.SUBTRACT, totalBasesInDict - 60, 2 + totalContigsInDict},
                {IntervalListTools.Action.SYMDIFF, totalBasesInDict - 61, 3 + totalContigsInDict},
                {IntervalListTools.Action.OVERLAPS, totalBasesInDict - 150, 2 + totalContigsInDict}
        };
    }

    @Test(dataProvider = "actionAndTotalBasesWithInvertData")
    public void testActionsWithInvert(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final IntervalList il = tester(action, true, false, false, scatterable, secondInput, false);
        Assert.assertEquals(il.getBaseCount(), bases, "unexpected number of bases found.");
        Assert.assertEquals(il.getIntervals().size(), intervals, "unexpected number of intervals found.");

        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.BASES, true, false, false, scatterable, secondInput), bases, "unexpected number of bases written to count_output file.");
        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.INTERVALS, true, false, false, scatterable, secondInput), intervals, "unexpected number of intervals written to count_output file.");
    }

    @DataProvider
    public Object[][] actionAndTotalBasesWithUniqueData() {
        return new Object[][]{
                {IntervalListTools.Action.CONCAT, 201, 2},
                {IntervalListTools.Action.UNION, 201, 2},
                {IntervalListTools.Action.INTERSECT, 140, 2},
                {IntervalListTools.Action.SUBTRACT, 60, 2},
                {IntervalListTools.Action.SYMDIFF, 61, 3},
                {IntervalListTools.Action.OVERLAPS, 150, 2}
        };
    }

    @Test(dataProvider = "actionAndTotalBasesWithUniqueData")
    public void testActionsWithUnique(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final IntervalList il = tester(action, false, true, false, scatterable, secondInput, false);
        Assert.assertEquals(il.getBaseCount(), bases, "unexpected number of bases found.");
        Assert.assertEquals(il.getIntervals().size(), intervals, "unexpected number of intervals found.");

        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.BASES, false, true, false, scatterable, secondInput), bases, "unexpected number of bases written to count_output file.");
        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.INTERVALS, false, true, false, scatterable, secondInput), intervals, "unexpected number of intervals written to count_output file.");
    }

    @DataProvider
    public Object[][] actionAndTotalBasesWithDontMergeAbuttingData() {
        return new Object[][]{
                {IntervalListTools.Action.CONCAT, 8, 3},
                {IntervalListTools.Action.UNION, 8, 3},
                {IntervalListTools.Action.INTERSECT, 8, 2},
                {IntervalListTools.Action.SUBTRACT, 0, 0},
                {IntervalListTools.Action.SYMDIFF, 0, 0},
                {IntervalListTools.Action.OVERLAPS, 8, 3}
        };
    }

    @Test(dataProvider = "actionAndTotalBasesWithDontMergeAbuttingData")
    public void testActionsWithDontMergeAbutting(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final IntervalList il = tester(action, false, true, true, abutting, abutting, false);
        Assert.assertEquals(il.getBaseCount(), bases, "unexpected number of bases found.");
        Assert.assertEquals(il.getIntervals().size(), intervals, "unexpected number of intervals found.");

        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.BASES, false, true, true, abutting, abutting), bases, "unexpected number of bases written to count_output file.");
        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.INTERVALS, false, true, true, abutting, abutting), intervals, "unexpected number of intervals written to count_output file.");
    }

    private IntervalList tester(final IntervalListTools.Action action) {
        return tester(action, DEFAULT_INVERT_ARG, DEFAULT_UNIQUE_ARG, DEFAULT_DONT_MERGE_ABUTTING, scatterable, secondInput, false);
    }

    private IntervalList tester(final IntervalListTools.Action action, final boolean invert, final boolean unique,
                                final boolean dontMergeAbutting, final Path input1, final Path input2, final boolean cloudOutput) {
        final String outputDirFullName = cloudOutput ? CLOUD_OUTPUT_DIR : null;
        final String prefix = "IntervalListTools";
        final PicardHtsPath output = PicardBucketUtils.getTempFilePath(outputDirFullName, prefix, INTERVAL_LIST_EXTENSION);
        final List<String> args = buildStandardTesterArguments(action, invert, unique, dontMergeAbutting, input1, input2);
        args.add("OUTPUT=" + output);
        Assert.assertEquals(runPicardCommandLine(args), 0);

        return IntervalList.fromPath(output.toPath());
    }

    private long testerCountOutput(IntervalListTools.Action action, IntervalListTools.Output outputValue) throws IOException {
        return testerCountOutput(action, outputValue, false, false, false, scatterable, secondInput, false);
    }

    private long testerCountOutput(IntervalListTools.Action action, IntervalListTools.Output outputValue, boolean invert,
                                   boolean unique, boolean dontMergeAbutting, Path input1, Path input2) throws IOException{
        return testerCountOutput(action, outputValue, invert, unique, dontMergeAbutting, input1, input2, false);
    }
    private long testerCountOutput(IntervalListTools.Action action, IntervalListTools.Output outputValue, boolean invert,
                                   boolean unique, boolean dontMergeAbutting, Path input1, Path input2,
                                   final boolean cloudOutput) throws IOException {
        final String outputDirFullName = cloudOutput ? CLOUD_OUTPUT_DIR : null;
        final String prefix = "IntervalListTools";
        final PicardHtsPath countOutput = PicardBucketUtils.getTempFilePath(outputDirFullName, prefix, ".txt");

        final List<String> args = buildStandardTesterArguments(action, invert, unique, dontMergeAbutting, input1, input2);
        args.add("OUTPUT_VALUE=" + outputValue);
        args.add("COUNT_OUTPUT=" + countOutput);

        Assert.assertEquals(runPicardCommandLine(args), 0);
        final Scanner reader = new Scanner(countOutput.toPath());
        return reader.nextLong();
    }

    // The caller is responsible for supplying the output argument(s)
    private List<String> buildStandardTesterArguments(final IntervalListTools.Action action, final boolean invert, final boolean unique,
                                                      final boolean dontMergeAbutting, final Path input1, final Path input2){
        final List<String> args = new ArrayList<>();

        args.add("ACTION=" + action.toString());
        args.add("INPUT=" + input1.toUri());

        if (action.takesSecondInput) {
            args.add("SECOND_INPUT=" + input2.toUri());
        } else {
            args.add("INPUT=" + input2.toUri()  );
        }

        if (invert) {
            args.add("INVERT=true");
        }

        if (unique) {
            args.add("UNIQUE=true");
        }

        if (dontMergeAbutting) {
            args.add("DONT_MERGE_ABUTTING=true");
        }

        return args;
    }

    /**
     * These are split out separately and private because the treatment of remainders causes different behavior whether
     * scatter count is specified (as similated here) versus whether scatter *content* is specified, as in in the IntervalListTools
     * integration tests
     * @return test cases with expected behavior for distributed remainder mode given a specified number of output lists
     */
    private static List<IntervalListScattererTest.Testcase> getRemainderTestcases() {
        final List<IntervalListScattererTest.Testcase> testCases = new ArrayList<>();
        testCases.add(new IntervalListScattererTest.Testcase(
                IntervalListScattererTest.LARGER_INTERVAL_FILE, 67, IntervalListScatterMode.INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER,
                LARGER_EXPECTED_WITH_REMAINDER_LISTS
        ));

        testCases.add(new IntervalListScattererTest.Testcase(
                IntervalListScattererTest.LARGER_INTERVAL_FILE, 20, IntervalListScatterMode.INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER,
                IntervalListScattererTest.LARGER_NO_REMAINDER_EXPECTED_LISTS
        ));
        return testCases;
    }

    // test scatter with different kinds of balancing.
    @DataProvider
    public static Iterator<Object[]> testScatterTestcases() {
        final List<IntervalListScattererTest.Testcase> testcases = new ArrayList<>(IntervalListScattererTest.getScatterTestcases());
        testcases.addAll(getRemainderTestcases());
        return testcases.stream().map(tc -> new Object[]{tc}).iterator();
    }

    private final List<File> dirsToDelete = new ArrayList<>();

    @AfterTest
    void deleteTempDirs() {
        for (File file : dirsToDelete) {
            try {
                IOUtil.recursiveDelete(file.toPath());
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    @Test(dataProvider = "testScatterTestcases")
    public void testScatter(final IntervalListScattererTest.Testcase tc) {

        final IntervalListScatterer scatterer = tc.mode.make();
        final List<IntervalList> scatter = scatterer.scatter(tc.source, tc.scatterWidth);
        Assert.assertEquals(scatter.size(), tc.expectedScatter.size());
        for (int i = 0; i < scatter.size(); i++) {
            Assert.assertEquals(scatter.get(i).getIntervals(), tc.expectedScatter.get(i).getIntervals(), "Problem with the " + i + " scatter");
        }

        final List<String> args = new ArrayList<>();

        args.add("ACTION=CONCAT");
        args.add("INPUT=" + tc.file.getAbsolutePath());

        args.add("SUBDIVISION_MODE=" + tc.mode);

        final File ilOutDir = IOUtil.createTempDir("IntervalListTools_lists").toFile();
        dirsToDelete.add(ilOutDir);

        if (tc.scatterWidth == 1) {
            final File subDir = new File(ilOutDir, "temp_1_of_1");
            Assert.assertTrue(subDir.mkdir(), "was unable to create directory");
            args.add("OUTPUT=" + new File(subDir, "scattered.interval_list"));
        } else {
            args.add("OUTPUT=" + ilOutDir);
        }

        args.add("SCATTER_COUNT=" + tc.scatterWidth);

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final String[] files = ilOutDir.list();
        Assert.assertNotNull(files);
        Arrays.sort(files);

        Assert.assertEquals(files.length, tc.expectedScatter.size());

        final Iterator<IntervalList> intervals = tc.expectedScatter.iterator();

        for (final String fileName : files) {
            final IntervalList intervalList = IntervalList.fromFile(new File(new File(ilOutDir, fileName), "scattered.interval_list"));
            final IntervalList expected = intervals.next();

            Assert.assertEquals(intervalList.getIntervals(), expected);
        }
    }

    @Test(dataProvider = "testScatterTestcases")
    public void testScatterByContent(final IntervalListScattererTest.Testcase tc) throws IOException {

        final List<String> args = new ArrayList<>();

        args.add("ACTION=CONCAT");
        args.add("INPUT=" + tc.file.getAbsolutePath());

        args.add("SUBDIVISION_MODE=" + tc.mode);

        final File ilOutDir = IOUtil.createTempDir("IntervalListTools_lists").toFile();
        dirsToDelete.add(ilOutDir);

        if (tc.scatterWidth == 1) {
            final File subDir = new File(ilOutDir, "temp_1_of_1");
            Assert.assertTrue(subDir.mkdir(), "was unable to create directory");
            args.add("OUTPUT=" + new File(subDir, "scattered.interval_list"));
        } else {
            args.add("OUTPUT=" + ilOutDir);
        }

        args.add("SCATTER_CONTENT=" + (int)Math.round((double)tc.mode.make().listWeight(tc.source) / tc.expectedScatter.size()));

        Assert.assertEquals(runPicardCommandLine(args), 0);

        final String[] files = ilOutDir.list();
        Assert.assertNotNull(files);
        Arrays.sort(files);

        Assert.assertEquals(files.length, tc.expectedScatter.size());

        final Iterator<IntervalList> intervals = tc.expectedScatter.iterator();

        for (final String fileName : files) {
            final IntervalList intervalList = IntervalList.fromFile(new File(new File(ilOutDir, fileName), "scattered.interval_list"));
            final IntervalList expected = intervals.next();

            Assert.assertEquals(intervalList.getIntervals(), expected);
        }
    }

    //should only take about 3 seconds now.
    @Test(timeOut = 40_000)
    public void testLargeScatters() throws IOException {
        final int scatterCount=1_000;
        final File ilOutDir = IOUtil.createTempDir("IntervalListTools_lists").toFile();
        dirsToDelete.add(ilOutDir);

        //scatter
        {
            final List<String> args = new ArrayList<>();

            args.add("SCATTER_COUNT=" + scatterCount);
            args.add("INPUT=" + largeScatterable);
            args.add("OUTPUT=" + ilOutDir);

            Assert.assertEquals(runPicardCommandLine(args), 0);
        }

        final List<File> files = Arrays.asList(Objects.requireNonNull(ilOutDir.listFiles()));
        Assert.assertEquals(files.size(), scatterCount);

        //gather
        final File ilOut = File.createTempFile("IntervalListTools", ".interval_list");
        ilOut.deleteOnExit();
        {
            final List<String> args = new ArrayList<>();
            files.forEach(f-> args.add("INPUT=" + f.toPath().resolve("scattered.interval_list").toAbsolutePath()));

            args.add("OUTPUT=" + ilOut);
            args.add("ACTION=UNION");

            Assert.assertEquals(runPicardCommandLine(args), 0);
        }
        final IntervalList gather = IntervalList.fromFile(ilOut);
        final IntervalList original = IntervalList.fromPath(largeScatterable).uniqued();

        Assert.assertEquals(gather, original);
    }

    @DataProvider
    public Object[][] combineAbuttingIntervals() {
        return new Object[][] {
                {false, abuttingCombined},
                {true, abuttingNotCombined},
        };
    }

    @Test(dataProvider = "combineAbuttingIntervals")
    public void testCombineAbuttingIntervals(boolean dontMergeAbutting, Path outputPath) throws IOException {
        // Test the default behavior of UNION, which is to combine abutting and overlapping intervals.
        //gather
        final File ilOut = File.createTempFile("IntervalListTools", ".interval_list");
        ilOut.deleteOnExit();
        final List<String> args = new ArrayList<>();
        args.add("INPUT=" + abutting.toUri());
        args.add("OUTPUT=" + ilOut);
        args.add("ACTION=UNION");
        args.add("DONT_MERGE_ABUTTING="+dontMergeAbutting);
        Assert.assertEquals(runPicardCommandLine(args), 0);
        final IntervalList gather = IntervalList.fromFile(ilOut);
        final IntervalList original = IntervalList.fromPath(outputPath);

        Assert.assertEquals(gather, original); // equal to expected output
    }

    // Check that we get the same output using the cloud (input/output) files as when we do everything locally.
    @Test(groups = "cloud")
    public void testCloudInputAgainstLocalInput() {
        // expected = local input files
        final IntervalList expectedIntervals = tester(DEFAULT_ACTION, DEFAULT_INVERT_ARG, DEFAULT_UNIQUE_ARG, DEFAULT_DONT_MERGE_ABUTTING, scatterable, secondInput, false);
        for (boolean outputInCloud : Arrays.asList(true, false)) {
            final IntervalList intervalListCloudInput = tester(DEFAULT_ACTION, DEFAULT_INVERT_ARG, DEFAULT_UNIQUE_ARG, DEFAULT_DONT_MERGE_ABUTTING, scatterableCloud.toPath(), secondInputCloud.toPath(), outputInCloud);

            Assert.assertEquals(expectedIntervals.size(), intervalListCloudInput.size());

            for (int i = 0; i < expectedIntervals.size(); i++) {
                final Interval local = expectedIntervals.getIntervals().get(i);
                final Interval cloud = intervalListCloudInput.getIntervals().get(i);
                Assert.assertEquals(local, cloud);
            }
        }
    }

    // Scatter intervals and write output file in the local filesystem (vs cloud) to be used as the "expected" output
    private String getScatteredIntervalListLocally(final int scatterCount) {
        final List<String> args = buildStandardTesterArguments(DEFAULT_ACTION, DEFAULT_INVERT_ARG, DEFAULT_UNIQUE_ARG, DEFAULT_DONT_MERGE_ABUTTING, scatterableCloud.toPath(), secondInputCloud.toPath());
        try {
            final Path expectedOutput = Files.createTempDirectory("scatter_test");
            args.add("SCATTER_COUNT=" + scatterCount);
            args.add("OUTPUT=" + expectedOutput.toUri());
            Assert.assertEquals(runPicardCommandLine(args), 0);
            return expectedOutput.toUri().getPath();
        } catch (IOException e){
            throw new PicardException("Encountered an IO Exception while running IntevalListTools", e);
        }
    }

    @Test(groups = "cloud")
    public void testScatterOutputInCloud() {
        final List<String> args = buildStandardTesterArguments(DEFAULT_ACTION, DEFAULT_INVERT_ARG, DEFAULT_UNIQUE_ARG, DEFAULT_DONT_MERGE_ABUTTING, scatterableCloud.toPath(), secondInputCloud.toPath());
        final int scatterCount = 3;

        // Files.createTempDirectory() seems to throw an error when used with a gs:// path
        final PicardHtsPath outputDir = new PicardHtsPath(CLOUD_OUTPUT_DIR + "scatter/");
        PicardIOUtils.deleteOnExit(outputDir.toPath());

        args.add("SCATTER_COUNT=" + scatterCount);
        args.add("OUTPUT=" + outputDir.toPath().toUri());
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final String expectedOutputDir = getScatteredIntervalListLocally(scatterCount);
        try {
            Assert.assertEquals(Files.list(outputDir.toPath()).count(), scatterCount);
            // All this is clunky because we create intermediate output directories e.g. temp_0001_of_3. Consider removing them.
            final List<Path> outputPaths = new ArrayList<>(scatterCount);
            for (int i = 1; i <= scatterCount; i++){
                outputPaths.add(IntervalListTools.createSubDirectoryAndGetScatterFile(outputDir, scatterCount, new DecimalFormat("0000").format(i)));
            }

            final List<Path> expectedOutput = new ArrayList<>(scatterCount);
            for (int i = 1; i <= scatterCount; i++){
                expectedOutput.add(IntervalListTools.createSubDirectoryAndGetScatterFile(new PicardHtsPath(expectedOutputDir), scatterCount, new DecimalFormat("0000").format(i)));
            }

            for (int j = 0; j < scatterCount; j++){
                final IntervalList output = IntervalList.fromPath(outputPaths.get(j));
                final IntervalList expected = IntervalList.fromPath(expectedOutput.get(j));
                Assert.assertEquals(output, expected);
            }
        } catch (IOException e){
            throw new PicardException("Encountered an IOException", e);
        }
    }

    /** Test {@link IntervalListTools.COUNT_OUTPUT} file in the cloud. Code copied from testActionsWithInvert **/
    @Test(dataProvider = "actionAndTotalBasesWithInvertData", groups = "cloud")
    public void testCountOutputInCloud(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final boolean cloudOutput = true;
        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.BASES, true, false, false, scatterable, secondInput, cloudOutput), bases, "unexpected number of bases written to count_output file.");
        Assert.assertEquals(testerCountOutput(action, IntervalListTools.Output.INTERVALS, true, false, false, scatterable, secondInput, cloudOutput), intervals, "unexpected number of intervals written to count_output file.");
    }
}
