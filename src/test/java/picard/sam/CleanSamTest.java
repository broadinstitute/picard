/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.sam.testers.CleanSamTester;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

public class CleanSamTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/CleanSam");
    private static final String qualityScore = "&/,&-.1/6/&&)&).)/,&0768)&/.,/874,&.4137572)&/&&,&1-&.0/&&*,&&&&&&&&&&18775799,&16:8775-56256/69::;0";

    public String getCommandLineProgramName() {
        return CleanSam.class.getSimpleName();
    }

    @Test(dataProvider = "testCleanSamDataProvider")
    public void testCleanSam(final String samFile, final String expectedCigar) throws IOException {
        final File cleanedFile = File.createTempFile(samFile + ".", ".sam");
        cleanedFile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + new File(TEST_DATA_DIR, samFile).getAbsolutePath(),
                "OUTPUT=" + cleanedFile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
        validator.setIgnoreWarnings(true);
        validator.setVerbose(true, 1000);
        validator.setErrorsToIgnore(Arrays.asList(SAMValidationError.Type.MISSING_READ_GROUP));
        SamReader samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(cleanedFile);
        final SAMRecord rec = samReader.iterator().next();
        samReader.close();
        Assert.assertEquals(rec.getCigarString(), expectedCigar);
        samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT).open(cleanedFile);
        final boolean validated = validator.validateSamFileVerbose(samReader, null);
        samReader.close();
        Assert.assertTrue(validated, "ValidateSamFile failed");
    }

    @DataProvider(name = "testCleanSamDataProvider")
    public Object[][] testCleanSamDataProvider() {
        return new Object[][]{
                {"simple_fits.sam", "100M"},
                {"simple_overhang.sam", "99M1S"},
                {"fits_with_deletion.sam", "91M2D9M"},
                {"overhang_with_deletion.sam", "91M2D8M1S"},
                {"trailing_insertion.sam", "99M1I"},
                {"long_trailing_insertion.sam", "90M10I"},
        };
    }

    //identical test case using the SamFileTester to generate that SAM file on the fly
    @Test(dataProvider = "testCleanSamTesterDataProvider")
    public void testCleanSamTester(final String originalCigar, final String expectedCigar, final int defaultChromosomeLength, final int alignStart) throws IOException {
        final CleanSamTester cleanSamTester = new CleanSamTester(expectedCigar, 100, defaultChromosomeLength);
        // NB: this will add in the mate cigar, when enabled in SamPairUtil, for additional validation
        cleanSamTester.addMappedPair(0, alignStart, alignStart, false, false, originalCigar, originalCigar, false, 50);
        cleanSamTester.runTest();
    }

    @DataProvider(name = "testCleanSamTesterDataProvider")
    public Object[][] testCleanSamTesterDataProvider() {
        return new Object[][]{
                {"100M", "100M", 101, 2}, // simple_filts.sam
                {"100M", "99M1S", 101, 3}, // simple_overhang.sam
                {"91M2D9M", "91M2D9M", 102, 1}, // fits_with_deletion.sam
                {"91M2D9M", "91M2D8M1S", 101, 1}, // overhang_with_deletion.sam
                {"99M1I", "99M1I", 101, 3}, // trailing_insertion.sam
                {"90M10I", "90M10I", 101, 3} // long_trailing_insertion.sam
        };
    }
}
