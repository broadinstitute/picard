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
package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.testers.ValidateSamTester;

import java.io.File;

public class MergeSamFilesTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/MergeSamFiles");

    public String getCommandLineProgramName() {
        return MergeSamFiles.class.getSimpleName();
    }
    /**
     * Confirm that unsorted input can result in coordinate sorted output, with index created.
     */
    @Test
    public void unsortedInputSortedOutputTest() throws Exception {
        final File unsortedInputTestDataDir = new File(TEST_DATA_DIR, "unsorted_input");
        final File mergedOutput = File.createTempFile("unsortedInputSortedOutputTest.", BamFileIoUtils.BAM_FILE_EXTENSION);
        mergedOutput.deleteOnExit();
        final String[] args = {
                "I=" + new File(unsortedInputTestDataDir, "1.sam").getAbsolutePath(),
                "I=" + new File(unsortedInputTestDataDir, "2.sam").getAbsolutePath(),
                "O=" + mergedOutput.getAbsolutePath(),
                "SO=coordinate"
        };
        final int mergeExitStatus = runPicardCommandLine(args);
        Assert.assertEquals(mergeExitStatus, 0);
        final SamReader reader = SamReaderFactory.makeDefault().open(mergedOutput);
        Assert.assertEquals(reader.getFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);

        new ValidateSamTester().assertSamValid(mergedOutput);
    }
}
