/*
 * The MIT License
 *
 * Copyright (c) 2016 Nils Homer
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
 *
 */

package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;

/**
 * Tests for FixVcfHeader.
 *
 * @author Nils Homer
 */
public class FixVcfHeaderTest {
    private static final File OUTPUT_DATA_PATH             = IOUtil.createTempDir("FixVcfHeaderTest", null);
    private static final File TEST_DATA_PATH               = new File("testdata/picard/vcf/FixVcfHeaderTest/");
    private static final File INPUT_VCF                    = new File(TEST_DATA_PATH, "input.vcf");
    private static final File OUTPUT_VCF                   = new File(TEST_DATA_PATH, "output.vcf");
    private static final File HEADER_VCF                   = new File(TEST_DATA_PATH, "header.vcf");
    private static final File HEADER_VCF_WITH_EXTRA_SAMPLE = new File(TEST_DATA_PATH, "header_with_extra_sample.vcf");

    private void runFixVcfHeader(final int checkFirstNRecords,
                                 final File replacementHeader,
                                 final boolean enforceSampleSamples) {
        final FixVcfHeader program = new FixVcfHeader();
        final File outputVcf = new File(OUTPUT_DATA_PATH, "output.vcf");
        outputVcf.deleteOnExit();

        program.INPUT      = INPUT_VCF;
        program.OUTPUT     = outputVcf;
        if (replacementHeader == null) {
            program.CHECK_FIRST_N_RECORDS = checkFirstNRecords;
        }
        else {
            program.HEADER = replacementHeader;
            program.ENFORCE_SAME_SAMPLES = enforceSampleSamples;
        }

        Assert.assertEquals(program.instanceMain(new String[0]), 0);

        IOUtil.assertFilesEqual(OUTPUT_VCF, outputVcf);
    }

    @Test(dataProvider = "testFixVcfHeaderDataProvider")
    public void testFixVcfHeader(final int checkFirstNRecords,
                                 final File replacementHeader,
                                 final boolean enforceSampleSamples) {
        runFixVcfHeader(checkFirstNRecords, replacementHeader, enforceSampleSamples);
    }

    @DataProvider(name="testFixVcfHeaderDataProvider")
    public Object[][] testFixVcfHeaderDataProvider() {
        return new Object[][] {
                {-1,                HEADER_VCF,                   true},
                {-1,                HEADER_VCF,                   false},
                {-1,                HEADER_VCF_WITH_EXTRA_SAMPLE, false},
                {-1,                null,                         true},
                {1,                 null,                         true},
                {Integer.MAX_VALUE, null,                         true}
        };
    }

    @Test(expectedExceptions=PicardException.class)
    public void testReplaceHeaderWithDifferentSamplesError() {
        runFixVcfHeader(-1, HEADER_VCF_WITH_EXTRA_SAMPLE, true);
    }
}