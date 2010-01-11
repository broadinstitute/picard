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

import net.sf.samtools.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

/**
 * Test new functionality that truncates sequence names at first whitespace in order to deal
 * with older BAMs that had spaces in sequence names.
 * 
 * @author alecw@broadinstitute.org
 */
public class SequenceNameTruncationAndValidationTest {
    private static File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam");
    @Test(expectedExceptions = {SAMException.class}, dataProvider = "badSequenceNames")
    public void testSequenceRecordThrowsWhenInvalid(final String sequenceName) {
        new SAMSequenceRecord(sequenceName, 123);
        Assert.fail("Should not reach here.");
    }
    @DataProvider(name="badSequenceNames")
    public Object[][] badSequenceNames() {
        return new Object[][] {
                {" "},
                {"\t"},
                {"\n"},
                {"="},
                {"Hi, Mom!"}
        };
    }

    @Test(dataProvider = "goodSequenceNames")
    public void testSequenceRecordPositiveTest(final String sequenceName) {
        new SAMSequenceRecord(sequenceName, 123);
    }
    @DataProvider(name="goodSequenceNames")
    public Object[][] goodSequenceNames() {
        return new Object[][] {
                {"Hi,@Mom!"}
        };
    }

    @Test(dataProvider = "samFilesWithSpaceInSequenceName")
    public void testSamSequenceTruncation(final String filename) {
        final SAMFileReader reader = new SAMFileReader(new File(TEST_DATA_DIR, filename));
        for (final SAMSequenceRecord sequence : reader.getFileHeader().getSequenceDictionary().getSequences()) {
            Assert.assertFalse(sequence.getSequenceName().contains(" "), sequence.getSequenceName());
        }
        for (final SAMRecord rec: reader) {
            Assert.assertFalse(rec.getReferenceName().contains(" "));
        }
    }
    @DataProvider(name="samFilesWithSpaceInSequenceName")
    public Object[][] samFilesWithSpaceInSequenceName() {
        return new Object[][] {
                {"sequenceWithSpace.sam"},
                {"sequenceWithSpace.bam"}
        };
    }

    @Test(expectedExceptions = {SAMFormatException.class})
    public void testBadRname() {
        final SAMFileReader reader = new SAMFileReader(new File(TEST_DATA_DIR, "readWithBadRname.sam"));
        for (final SAMRecord rec: reader) {
        }
        Assert.fail("Should not reach here.");
    }
}
