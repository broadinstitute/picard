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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class AddOATagTest {
    private static final String TEST_DIR = "testdata/picard/sam/";
    private static final String OA_TEST_DIR = TEST_DIR + "AddOATag/";

    @DataProvider(name="testOAData")
    Object[][] testOAData() {
        return new Object[][]{
                new Object[]{new File(OA_TEST_DIR + "aligned_withoutOA.sam"), new File(OA_TEST_DIR + "aligned_withOA.sam")},
        };
    }

    @DataProvider(name="testIntervalListData")
    Object[][] testIntervalListData() {
        return new Object[][]{
                new Object[]{new File(TEST_DIR + "aligned.sam"), null, new File(OA_TEST_DIR + "aligned_no_intervals.sam")},
                new Object[]{new File(TEST_DIR + "aligned.sam"), new File("testdata/picard/analysis/directed/CollectHsMetrics/chrM.interval_list"), new File(TEST_DIR + "aligned.sam")},
                new Object[]{new File(TEST_DIR + "aligned.sam"), new File(TEST_DIR + "contiguous.interval_list"), new File(OA_TEST_DIR + "aligned_with_intervals.sam")}
        };
    }

    @DataProvider(name="testOverwriteData")
    Object[][] overwriteTestData() {
        return new Object[][]{
                new Object[]{new File(TEST_DIR + "aligned.sam"), new File(TEST_DIR + "contiguous.interval_list"), true, new File(OA_TEST_DIR + "aligned_with_intervals.sam")},
                new Object[]{new File(TEST_DIR + "aligned.sam"), null, false, new File(OA_TEST_DIR + "aligned_two_oa_tags.sam")}
        };
    }

    @Test(dataProvider = "testOAData")
    public void testWritingOATag(final File testSam, final File truthSam) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(testSam, clpOutput, null, null);

        validateOATag(clpOutput, truthSam);
    }

    @Test(dataProvider = "testIntervalListData")
    public void testIntervalList(final File inputSam, final File intervalList, final File truthSam) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(inputSam, clpOutput, intervalList, null);

        validateOATag(clpOutput, truthSam);
    }

    @Test(dataProvider = "testOverwriteData")
    public void testOverWrite(final File inputSam, final File interval_list, final Boolean overWriteTag, final File truthSam) throws IOException {
        final File firstPassOutput = File.createTempFile("FirstPassAddOATag", ".bam");
        firstPassOutput.deleteOnExit();
        final File secondPassOutput = File.createTempFile("SecondPassAddOATag", ".bam");
        secondPassOutput.deleteOnExit();

        // make first pass bam that only has one value per OA tag
        runAddOATag(inputSam, firstPassOutput, interval_list, true);

        // make bam we want to test the overwrite option with
        runAddOATag(firstPassOutput, secondPassOutput, interval_list, overWriteTag);

        validateOATag(secondPassOutput, truthSam);
    }

    private void runAddOATag(final File inputSam, final File output, final File intervalList, Boolean overwriteTag) throws IOException {
        final ArrayList<String> args = new ArrayList<String>(){
            {
                add("INPUT=" + inputSam);
                add("OUTPUT=" + output);
            }
        };
        if (intervalList != null) {
            args.add("INTERVAL_LIST=" + intervalList);
        }
        if (overwriteTag != null) {
            args.add("OVERWRITE_TAG=" + overwriteTag);
        }
        AddOATag addOATag = new AddOATag();
        Assert.assertEquals(addOATag.instanceMain(args.toArray(new String[args.size()])), 0, "Running addOATag did not succeed");
    }

    private void validateOATag(final File testSam, final File truthSam) {
        final ArrayList<String> truthOAValues = new ArrayList<>();
        final ArrayList<String> testOAValues = new ArrayList<>();

        SAMRecordIterator iterator = SamReaderFactory.makeDefault().open(truthSam).iterator();
        while (iterator.hasNext()){
            SAMRecord rec = iterator.next();
            truthOAValues.add(rec.getStringAttribute("OA"));
        }

        iterator = SamReaderFactory.makeDefault().open(testSam).iterator();
        while (iterator.hasNext()){
            SAMRecord rec = iterator.next();
            testOAValues.add(rec.getStringAttribute("OA"));
        }

        Assert.assertEquals(testOAValues, truthOAValues);
    }
}
