/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class MergeSamFilesTest {
    private static File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam/MergeSamFiles");

    /**
     * Confirm that unsorted input can result in coordinate sorted output, with index created.
     */
    @Test
    public void unsortedInputSortedOutputTest() throws Exception {
        File unsortedInputTestDataDir = new File(TEST_DATA_DIR, "unsorted_input");
        File mergedOutput = File.createTempFile("unsortedInputSortedOutputTest.", ".bam");
        mergedOutput.deleteOnExit();
        String[] args = {
                "I=" + new File(unsortedInputTestDataDir, "1.sam").getAbsolutePath(),
                "I=" + new File(unsortedInputTestDataDir, "2.sam").getAbsolutePath(),
                "O=" + mergedOutput.getAbsolutePath(),
                "SO=coordinate"
        };
        final int mergeExitStatus = new MergeSamFiles().instanceMain(args);
        Assert.assertEquals(mergeExitStatus, 0);
        final SAMFileReader reader = new SAMFileReader(mergedOutput);
        Assert.assertEquals(reader.getFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);

        final int validateExitStatus = new ValidateSamFile().instanceMain(new String[] {"I=" + mergedOutput.getAbsolutePath()});
        Assert.assertEquals(validateExitStatus, 0);
    }
}
