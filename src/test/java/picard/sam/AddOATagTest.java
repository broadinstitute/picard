package picard.sam;

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

import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class AddOATagTest {
    private static final Path TEST_DIR = Paths.get("testdata/picard/sam/");
    private static final Path OA_TEST_DIR = TEST_DIR.resolve("AddOATag/");

    @DataProvider()
    Object[][] testOAData() {
        return new Object[][]{
                new Object[]{OA_TEST_DIR.resolve("aligned_withoutOA.sam").toFile(), OA_TEST_DIR.resolve("aligned_withOA.sam").toFile()},
        };
    }

    @DataProvider()
    Object[][] testIntervalListData() {
        return new Object[][]{
                new Object[]{TEST_DIR.resolve("aligned.sam").toFile(), null, OA_TEST_DIR.resolve("aligned_no_intervals.sam").toFile()},
                new Object[]{TEST_DIR.resolve("aligned.sam").toFile(), new File("testdata/picard/analysis/directed/CollectHsMetrics/chrM.interval_list"), TEST_DIR.resolve("aligned.sam").toFile()},
                new Object[]{TEST_DIR.resolve("aligned.sam").toFile(), TEST_DIR.resolve("contiguous.interval_list").toFile(), OA_TEST_DIR.resolve("aligned_with_intervals.sam").toFile()}
        };
    }

    @Test(dataProvider = "testOAData")
    public void testWritingOATag(final File testSam, final File truthSam) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(testSam, clpOutput, null);
        validateOATag(clpOutput, truthSam);
    }

    @Test(dataProvider = "testIntervalListData")
    public void testIntervalList(final File inputSam, final File intervalList, final File truthSam) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(inputSam, clpOutput, intervalList);
        validateOATag(clpOutput, truthSam);
    }

    private static void runAddOATag(final File inputSam, final File output, final File intervalList) {
        final List<String> args = new ArrayList<>(Arrays.asList(
          "INPUT=" + inputSam,
          "OUTPUT=" + output
        ));

        if (intervalList != null) {
            args.add("INTERVAL_LIST=" + intervalList);
        }

        AddOATag addOATag = new AddOATag();
        Assert.assertEquals(addOATag.instanceMain(args.toArray(new String[0])), 0, "Running addOATag did not succeed");
    }

    private static void validateOATag(final File testSam, final File truthSam) throws IOException {
        final ArrayList<String> truthOAValues = new ArrayList<>();
        final ArrayList<String> testOAValues = new ArrayList<>();

        try(SamReader readerPre = SamReaderFactory.makeDefault().open(truthSam)) {
            readerPre.forEach(rec -> truthOAValues.add(rec.getStringAttribute(SAMTag.OA.name())));
        }

        try (SamReader readerPost = SamReaderFactory.makeDefault().open(testSam)) {
            readerPost.forEach(rec -> testOAValues.add(rec.getStringAttribute(SAMTag.OA.name())));
        }

        Assert.assertEquals(testOAValues, truthOAValues);
    }
}
