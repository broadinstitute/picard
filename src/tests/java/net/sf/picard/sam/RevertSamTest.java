/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.samtools.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: ktibbett
 * Date: Jul 20, 2010
 * Time: 10:27:58 AM
 * To change this template use File | Settings | File Templates.
 */
public class RevertSamTest {

    public static final String basicSamToRevert = "testdata/net/sf/picard/sam/revert_sam_basic.sam";
    public static final String negativeTestSamToRevert = "testdata/net/sf/picard/sam/revert_sam_negative.sam";

    private static final String revertedQualities  =
        "11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111";

    private static final String unmappedRead = "both_reads_present_only_first_aligns/2";

    @Test(dataProvider="positiveTestData")
    public void basicPositiveTests(SAMFileHeader.SortOrder so, boolean removeDuplicates, boolean removeAlignmentInfo,
                                   boolean restoreOriginalQualities, String sample, String library,
                                   List<String> attributesToClear) throws Exception {

        File output = File.createTempFile("reverted", ".sam");
        RevertSam reverter = new RevertSam();
        String args[] = new String[5 + (so != null ? 1 : 0) + attributesToClear.size() + (sample != null ? 1 : 0) + (library != null ? 1 : 0)];
        int index = 0;
        args[index++] = "INPUT=" + basicSamToRevert;
        args[index++] = "OUTPUT=" + output.getAbsolutePath();
        if (so != null) {
            args[index++] = "SORT_ORDER=" + so.name();
        }
        args[index++] = "REMOVE_DUPLICATE_INFORMATION=" + removeDuplicates;
        args[index++] = "REMOVE_ALIGNMENT_INFORMATION=" + removeAlignmentInfo;
        args[index++] = "RESTORE_ORIGINAL_QUALITIES=" + restoreOriginalQualities;
        if (sample != null) {
            args[index++] = "SAMPLE_ALIAS=" + sample;
        }
        if (library != null) {
            args[index++] = "LIBRARY_NAME=" + library;
        }
        for (String attr : attributesToClear) {
            args[index++] = "ATTRIBUTE_TO_CLEAR=" + attr;
        }
        reverter.instanceMain(args);

        SAMFileReader reader = new SAMFileReader(output);
        SAMFileHeader header = reader.getFileHeader();
        Assert.assertEquals(header.getSortOrder(), SAMFileHeader.SortOrder.queryname);
        Assert.assertEquals(header.getProgramRecords().size(), removeAlignmentInfo ? 0 : 1);
        for (SAMReadGroupRecord rg : header.getReadGroups()) {
            if (sample != null) {
                Assert.assertEquals(rg.getSample(), sample);
            }
            else {
                Assert.assertEquals(rg.getSample(), "Hi,Mom!");
            }
            if (library != null) {
                Assert.assertEquals(rg.getLibrary(), library);
            }
            else {
                Assert.assertEquals(rg.getLibrary(), "my-library");
            }
        }
        SAMRecordIterator it = reader.iterator();
        while (it.hasNext()) {
            SAMRecord rec = it.next();

            if (removeDuplicates) {
                Assert.assertFalse(rec.getDuplicateReadFlag(),
                    "Duplicates should have been removed: " + rec.getReadName());
            }

            if (removeAlignmentInfo) {
                Assert.assertTrue(rec.getReadUnmappedFlag(),
                     "Alignment info should have been removed: " + rec.getReadName());
            }

            if (restoreOriginalQualities && !unmappedRead.equals(
                    rec.getReadName() + "/" + (rec.getFirstOfPairFlag() ? "1" : "2"))) {

                Assert.assertEquals(rec.getBaseQualityString(), revertedQualities);
            }
            else {
                Assert.assertNotSame(rec.getBaseQualityString(), revertedQualities);
            }

            for (SAMRecord.SAMTagAndValue attr : rec.getAttributes()) {
                if (removeAlignmentInfo || (!attr.tag.equals("PG") && !attr.tag.equals("NM")
                    && !attr.tag.equals("MQ"))) {
                    Assert.assertFalse(reverter.ATTRIBUTE_TO_CLEAR.contains(attr.tag),
                            attr.tag + " should have been cleared.");
                }
            }
        }
    }


    @DataProvider(name="positiveTestData")
    public Object[][] getPostitiveTestData() {
        return new Object[][] {
                {null, true, true, true, null, null, Collections.EMPTY_LIST},
                {SAMFileHeader.SortOrder.queryname, true, true, true, "Hey,Dad!", null, Arrays.asList("XT")},
                {null, false, true, false, "Hey,Dad!", "NewLibraryName", Arrays.asList("XT")},
                {null, false, false, false, null, null, Collections.EMPTY_LIST}
        };
    }


    @Test(dataProvider="negativeTestData", expectedExceptions = {PicardException.class})
    public void basicNegativeTest(String sample, String library) throws Exception {

        File output = File.createTempFile("bad", ".sam");
        RevertSam reverter = new RevertSam();
        String args[] = new String[2 + (sample != null ? 1 : 0) + (library != null ? 1 : 0)];
        int index = 0;
        args[index++] = "INPUT=" + negativeTestSamToRevert;
        args[index++] = "OUTPUT=" + output.getAbsolutePath();
        if (sample != null) {
            args[index++] = "SAMPLE_ALIAS=" + sample;
        }
        if (library != null) {
            args[index++] = "LIBRARY_NAME=" + library;
        }
        reverter.instanceMain(args);
        Assert.fail("Negative test should have thrown an exception and didn't");
    }

    @DataProvider(name="negativeTestData")
    public Object[][] getNegativeTestData() {
        return new Object[][] {
                {"NewSample", null},
                {null, "NewLibrary"},
                {"NewSample", "NewLibrary"}
        };
    }
}
