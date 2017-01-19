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
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;
import java.io.IOException;

/**
 * Tests for FixVcfHeader.
 *
 * @author Nils Homer
 */
public class FixVcfHeaderTest {
    private File INPUT_VCF;
    private File OUTPUT_VCF;
    private File HEADER_VCF;
    private File HEADER_VCF_WITH_EXTRA_SAMPLE;

    @BeforeTest
    void setup() throws IOException {
        File OUTPUT_DATA_PATH = IOUtil.createTempDir("FixVcfHeaderTest", null);
        OUTPUT_DATA_PATH.deleteOnExit();
        final File testDataPath      = new File("testdata/picard/vcf/FixVcfHeaderTest/");
        INPUT_VCF                    = new File(testDataPath, "input.vcf");
        OUTPUT_VCF                   = new File(testDataPath, "output.vcf");
        HEADER_VCF                   = new File(testDataPath, "header.vcf");
        HEADER_VCF_WITH_EXTRA_SAMPLE = new File(testDataPath, "header_with_extra_sample.vcf");
    }

    private void runFixVcfHeader(final int checkFirstNRecords,
                                 final File replacementHeader,
                                 final boolean enforceSampleSamples) throws IOException {
        final FixVcfHeader program = new FixVcfHeader();
        final File outputVcf = VcfTestUtils.createTemporaryIndexedVcfFile("output.", ".vcf");

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
                                 final boolean enforceSampleSamples) throws IOException {
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
    public void testReplaceHeaderWithDifferentSamplesError() throws IOException {
        runFixVcfHeader(-1, HEADER_VCF_WITH_EXTRA_SAMPLE, true);
    }
}