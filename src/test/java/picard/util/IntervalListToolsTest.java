package picard.util;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.*;
import picard.cmdline.CommandLineProgramTest;
import picard.cmdline.PicardCommandLineTest;
import picard.util.IntervalList.IntervalListScatterMode;
import picard.util.IntervalList.IntervalListScatterer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * Created by farjoun on 10/22/17.
 */

public class IntervalListToolsTest extends CommandLineProgramTest {
    private final File TEST_DATA_DIR = new File("testdata/picard/util/");
    private final File scatterable = new File(TEST_DATA_DIR, "scatterable.interval_list");
    private final File secondInput = new File(TEST_DATA_DIR, "secondInput.interval_list");

    @Test
    public void testSecondInputValidation() {
        IntervalListTools intervalListTools = new IntervalListTools();
        String[] errors = intervalListTools.customCommandLineValidation();
        Assert.assertNull(errors);

        intervalListTools.SECOND_INPUT = new ArrayList<>();
        errors = intervalListTools.customCommandLineValidation();
        Assert.assertNull(errors);

        intervalListTools.SECOND_INPUT.add(new File("fakefile"));
        errors = intervalListTools.customCommandLineValidation();
        Assert.assertTrue(errors.length == 1);
    }

    @Override
    public String getCommandLineProgramName() {
        return IntervalListTools.class.getSimpleName();
    }

    @DataProvider
    public Iterator<Object[]> ActionsTest() {
        return Arrays.stream(IntervalListTools.Action.values()).map(a -> new Object[]{a}).iterator();
    }

    // test that all actions work. but not test output at all.
    @Test(dataProvider = "ActionsTest")
    public void testAllActions(final IntervalListTools.Action action) throws IOException {
        final File ilOut = File.createTempFile("IntervalListTools", "interval_list");
        ilOut.deleteOnExit();

        final List<String> args = new ArrayList<>();

        args.add("ACTION=" + action.toString());
        args.add("INPUT=" + scatterable);

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
                {IntervalListTools.Action.OVERLAPS, 150, 2},
        };
    }

    @Test(dataProvider = "actionAndTotalBasesData")
    public void testActions(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final IntervalList il = tester(action);
        Assert.assertEquals(il.getBaseCount(), bases, "unexpected number of bases found.");
        Assert.assertEquals(il.getIntervals().size(), intervals, "unexpected number of intervals found.");
    }

    @DataProvider
    public Object[][] actionAndTotalBasesWithInvertData() {
        final long totalBasesInDict = IntervalList.fromFile(secondInput).getHeader().getSequenceDictionary().getReferenceLength();
        final int totalContigsInDict = IntervalList.fromFile(secondInput).getHeader().getSequenceDictionary().size();
        return new Object[][]{
                {IntervalListTools.Action.CONCAT, totalBasesInDict - 201, 2 + totalContigsInDict},
                {IntervalListTools.Action.UNION, totalBasesInDict - 201, 2 + totalContigsInDict},
                {IntervalListTools.Action.INTERSECT, totalBasesInDict - 140, 2 + totalContigsInDict},
                {IntervalListTools.Action.SUBTRACT, totalBasesInDict - 60, 2 + totalContigsInDict},
                {IntervalListTools.Action.SYMDIFF, totalBasesInDict - 61, 3 + totalContigsInDict},
                {IntervalListTools.Action.OVERLAPS, totalBasesInDict - 150, 2 + totalContigsInDict},
        };
    }

    @Test(dataProvider = "actionAndTotalBasesWithInvertData")
    public void testActionsWithInvert(final IntervalListTools.Action action, final long bases, final int intervals) throws IOException {
        final IntervalList il = tester(action, true);
        Assert.assertEquals(il.getBaseCount(), bases, "unexpected number of bases found.");
        Assert.assertEquals(il.getIntervals().size(), intervals, "unexpected number of intervals found.");
    }

    private IntervalList tester(IntervalListTools.Action action) throws IOException {
        return tester(action, false);
    }

    private IntervalList tester(IntervalListTools.Action action, boolean invert) throws IOException {
        final File ilOut = File.createTempFile("IntervalListTools", "interval_list");
        ilOut.deleteOnExit();

        final List<String> args = new ArrayList<>();

        args.add("ACTION=" + action.toString());
        args.add("INPUT=" + scatterable);

        if (action.takesSecondInput) {
            args.add("SECOND_INPUT=" + secondInput);
        } else {
            args.add("INPUT=" + secondInput);
        }

        if (invert) {
            args.add("INVERT=true");
        }
        args.add("OUTPUT=" + ilOut);

        Assert.assertEquals(runPicardCommandLine(args), 0);

        return IntervalList.fromFile(ilOut);
    }

    // test scatter with different kinds of balancing.
    @DataProvider
    public static Iterator<Object[]> testScatterTestcases() {
        return IntervalListScattererTest.testScatterTestcases();
    }

    private final List<File> dirsToDelete = new ArrayList<>();

    @AfterTest
    void deleteTempDirs() {
        dirsToDelete.forEach(TestUtil::recursiveDelete);
    }

    @Test(dataProvider = "testScatterTestcases")
    public void testScatter(final IntervalListScattererTest.Testcase tc) throws IOException {

        final IntervalListScatterer scatterer = tc.mode.make();
        final List<IntervalList> scatter = scatterer.scatter(tc.source, tc.scatterWidth);
        Assert.assertEquals(scatter.size(), tc.expectedScatter.size());
        for (int i = 0; i < scatter.size(); i++) {
            Assert.assertEquals(scatter.get(i).uniqued().getIntervals(), tc.expectedScatter.get(i).getIntervals(), "Problem with the " + i + " scatter");
        }

        final List<String> args = new ArrayList<>();

        args.add("ACTION=UNION");
        args.add("INPUT=" + tc.file.getAbsolutePath());

        args.add("SUBDIVISION_MODE=" + tc.mode);

        final File ilOutDir = IOUtil.createTempDir("IntervalListTools", "lists");
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
}
