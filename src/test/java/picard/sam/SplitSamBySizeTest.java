/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.PicardException;
import picard.sam.testers.ValidateSamTester;

import java.io.File;
import java.io.IOException;

/**
 * Tests for SAMSplitter
 */
public class SplitSamBySizeTest extends CommandLineProgramTest {
    private static final String TEST_DATA_DIR = "testdata/picard/sam/bam2fastq";
    private static final File PAIRED_FILE = new File(TEST_DATA_DIR + "/paired/ok/sorted-pair.sam");
    private static final ValidateSamTester VALIDATE_SAM_TESTER = new ValidateSamTester();

    public String getCommandLineProgramName() {
        return SplitSamBySize.class.getSimpleName();
    }

    @Test
    public void testOkFile() throws IOException {
        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";

        final String [] args = new String[]{
                "INPUT=" + PAIRED_FILE.getAbsolutePath(),
                "TOTAL_READS_IN_INPUT=10",
                "SPLIT_TO_N_TEMPLATES=3",
                "OUTPUT_DIR=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out1 = new File(tmpDir + "1.bam");
        final File out2 = new File(tmpDir + "2.bam");
        final SamReader reader1 = SamReaderFactory.makeDefault().open(out1);
        final SamReader reader2 = SamReaderFactory.makeDefault().open(out2);

        Assert.assertEquals(reader1.getFileHeader(), reader2.getFileHeader());
        VALIDATE_SAM_TESTER.assertSamValid(out1);
        VALIDATE_SAM_TESTER.assertSamValid(out2);
        Assert.assertEquals(reader1.iterator().toList().size(), 6);
        Assert.assertEquals(reader2.iterator().toList().size(), 4);
    }

    @Test
    public void testOneOutput() throws IOException {
        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";

        final String [] args = new String[]{
                "INPUT=" + PAIRED_FILE.getAbsolutePath(),
                "TOTAL_READS_IN_INPUT=10",
                "SPLIT_TO_N_TEMPLATES=5",
                "OUTPUT_DIR=" + tmpDir
        };

        runPicardCommandLine(args);

        final File out = new File(tmpDir + "1.bam");
        final SamReader outReader = SamReaderFactory.makeDefault().open(out);
        final SamReader inReader = SamReaderFactory.makeDefault().open(PAIRED_FILE);
        VALIDATE_SAM_TESTER.assertSamValid(out);
        Assert.assertEquals(inReader.iterator().toList().size(), outReader.iterator().toList().size());
    }

    @Test(expectedExceptions = PicardException.class)
    public void testZeroOutput() throws IOException {
        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";

        final String [] args = new String[]{
                "INPUT=" + PAIRED_FILE.getAbsolutePath(),
                "TOTAL_READS_IN_INPUT=10",
                "SPLIT_TO_N_TEMPLATES=100",
                "OUTPUT_DIR=" + tmpDir
        };

        runPicardCommandLine(args);
    }
}