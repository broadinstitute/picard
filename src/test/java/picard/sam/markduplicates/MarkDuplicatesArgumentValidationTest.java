/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.sam.markduplicates.MarkDuplicatesEndToEndRunner.Result;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Tests how {@link MarkDuplicates} rejects arguments and inputs it cannot use: arguments the parser refuses
 * (mutually exclusive ASSUME_SORTED and ASSUME_SORT_ORDER, missing required arguments, values that do not parse),
 * which make it exit with status 1 without writing anything, and the errors raised once it runs (an input that is
 * neither coordinate nor queryname sorted, files it cannot read or write).
 */
public class MarkDuplicatesArgumentValidationTest {

    private static MarkDuplicatesEndToEndRunner input(final SAMFileHeader.SortOrder sortOrder) {
        final MarkDuplicatesEndToEndRunner input = new MarkDuplicatesEndToEndRunner(sortOrder).readGroup("rg", "lib");
        input.pair("a", "rg", 0, 100, false, "10M", 0, 291, true, "10M", 30);
        input.pair("b", "rg", 0, 100, false, "10M", 0, 291, true, "10M", 20);
        return input;
    }

    @DataProvider(name = "rejectedArguments")
    public Object[][] rejectedArguments() {
        return new Object[][]{
                {new String[]{"ASSUME_SORTED=true", "ASSUME_SORT_ORDER=coordinate"}},
                {new String[]{"ASSUME_SORTED=false", "ASSUME_SORT_ORDER=coordinate"}},
                {new String[]{"TAGGING_POLICY=Sometimes"}},
                {new String[]{"DUPLICATE_SCORING_STRATEGY=BEST"}},
                {new String[]{"OPTICAL_DUPLICATE_PIXEL_DISTANCE=near"}},
                {new String[]{"MAX_OPTICAL_DUPLICATE_SET_SIZE=1.5"}},
                {new String[]{"SORTING_COLLECTION_SIZE_RATIO=small"}},
                {new String[]{"NO_SUCH_ARGUMENT=1"}},
        };
    }

    @Test(dataProvider = "rejectedArguments")
    public void testRejectedArgumentsExitWithOne(final String[] args) {
        final File dir = IOUtil.createTempDir("MarkDuplicatesRejected").toFile();
        try {
            final File in = input(SAMFileHeader.SortOrder.coordinate).write(new File(dir, "input.sam"));
            final Result result = MarkDuplicatesEndToEndRunner.runOn(dir, Collections.singletonList(in), args);
            Assert.assertEquals(result.exitCode, 1);
            Assert.assertFalse(new File(dir, "output.sam").exists());
            Assert.assertFalse(new File(dir, "metrics.txt").exists());
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
    }

    @DataProvider(name = "missingRequiredArguments")
    public Object[][] missingRequiredArguments() {
        return new Object[][]{{"OUTPUT"}, {"METRICS_FILE"}};
    }

    @Test(dataProvider = "missingRequiredArguments")
    public void testMissingRequiredArgumentExitsWithOne(final String missing) {
        final File dir = IOUtil.createTempDir("MarkDuplicatesMissing").toFile();
        try {
            final File in = input(SAMFileHeader.SortOrder.coordinate).write(new File(dir, "input.sam"));
            final List<String> args = new ArrayList<>();
            args.add("INPUT=" + in.getAbsolutePath());
            if (!missing.equals("OUTPUT")) {
                args.add("OUTPUT=" + new File(dir, "output.sam").getAbsolutePath());
            }
            if (!missing.equals("METRICS_FILE")) {
                args.add("METRICS_FILE=" + new File(dir, "metrics.txt").getAbsolutePath());
            }
            Assert.assertEquals(new MarkDuplicates().instanceMain(args.toArray(new String[0])), 1);
            Assert.assertFalse(new File(dir, "output.sam").exists());
            Assert.assertFalse(new File(dir, "metrics.txt").exists());
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
    }

    /** The header's sort order, extra arguments, and the unusable sort order the error must name. */
    @DataProvider(name = "unusableSortOrders")
    public Object[][] unusableSortOrders() {
        return new Object[][]{
                {SAMFileHeader.SortOrder.unsorted, new String[0], SAMFileHeader.SortOrder.unsorted},
                {SAMFileHeader.SortOrder.unknown, new String[0], SAMFileHeader.SortOrder.unknown},
                {SAMFileHeader.SortOrder.coordinate, new String[]{"ASSUME_SORT_ORDER=unsorted"}, SAMFileHeader.SortOrder.unsorted},
                {SAMFileHeader.SortOrder.coordinate, new String[]{"ASSUME_SORT_ORDER=unknown"}, SAMFileHeader.SortOrder.unknown},
                {SAMFileHeader.SortOrder.coordinate, new String[]{"ASSUME_SORT_ORDER=duplicate"}, SAMFileHeader.SortOrder.duplicate},
        };
    }

    @Test(dataProvider = "unusableSortOrders")
    public void testInputMustBeCoordinateOrQuerynameSorted(final SAMFileHeader.SortOrder headerOrder, final String[] args,
                                                           final SAMFileHeader.SortOrder unusableOrder) {
        final File dir = IOUtil.createTempDir("MarkDuplicatesSortOrder").toFile();
        try {
            final File in = input(headerOrder).write(new File(dir, "input.sam"));
            final PicardException e = Assert.expectThrows(PicardException.class,
                    () -> MarkDuplicatesEndToEndRunner.runOn(dir, Collections.singletonList(in), args));
            Assert.assertTrue(e.getMessage().contains(unusableOrder.name()), e.getMessage());
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
    }

    @DataProvider(name = "unusableFiles")
    public Object[][] unusableFiles() {
        return new Object[][]{{"INPUT"}, {"OUTPUT"}, {"METRICS_FILE"}};
    }

    /** An input that does not exist, or an output or metrics file in a directory that does not exist. */
    @Test(dataProvider = "unusableFiles")
    public void testUnreadableInputOrUnwritableOutput(final String which) {
        final File dir = IOUtil.createTempDir("MarkDuplicatesFiles").toFile();
        try {
            final File missingDir = new File(dir, "missing");
            final File in = which.equals("INPUT") ? new File(missingDir, "input.sam") :
                    input(SAMFileHeader.SortOrder.coordinate).write(new File(dir, "input.sam"));
            final File out = new File(which.equals("OUTPUT") ? missingDir : dir, "output.sam");
            final File metrics = new File(which.equals("METRICS_FILE") ? missingDir : dir, "metrics.txt");
            final List<String> args = Arrays.asList("INPUT=" + in.getAbsolutePath(), "OUTPUT=" + out.getAbsolutePath(),
                    "METRICS_FILE=" + metrics.getAbsolutePath());
            try {
                new MarkDuplicates().instanceMain(args.toArray(new String[0]));
                Assert.fail("expected a SAMException");
            } catch (final SAMException e) {
                Assert.assertTrue(e.getMessage().contains(missingDir.getName()), e.getMessage());
            }
            Assert.assertFalse(out.exists());
            Assert.assertFalse(metrics.exists());
        } finally {
            MarkDuplicatesEndToEndRunner.deleteTempDir(dir);
        }
    }
}
