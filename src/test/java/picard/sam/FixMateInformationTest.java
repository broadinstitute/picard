/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.samtools.SAMException;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.IOException;

public class FixMateInformationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/FixMateInformation");
    private static final String MISSING_MATE_TEST = "missingMate.sam";

    public String getCommandLineProgramName() {
        return FixMateInformation.class.getSimpleName();
    }

    public int missingMateTestHelper(final boolean ignoreMissingMates) throws IOException {
        final File inSamFile = new File(TEST_DATA_DIR, MISSING_MATE_TEST);
        final File outSamFile = File.createTempFile("outMissingMateTest", "sam");
        outSamFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + inSamFile.getAbsolutePath(),
                "OUTPUT=" + outSamFile.getAbsolutePath(),
                "IGNORE_MISSING_MATES=" + ignoreMissingMates
        };

        return new FixMateInformation().instanceMain(args);
    }

    @Test
    public void ignoreMissingMateTest() throws IOException {
        Assert.assertEquals(missingMateTestHelper(true), 0);
    }

    @Test(expectedExceptions = SAMException.class)
    public void ignoreMissingMateExceptionTest() throws IOException {
        missingMateTestHelper(false);
    }
}
