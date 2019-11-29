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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.testng.Assert;
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
    void setup() {
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
        final File outputVcf = VcfTestUtils.createTemporaryIndexedFile("output.", ".vcf");

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

        final VCFFileReader actualReader   = new VCFFileReader(OUTPUT_VCF, false);
        final VCFFileReader expectedReader = new VCFFileReader(outputVcf, false);

        // Check that the headers match (order does not matter
        final VCFHeader actualHeader   = actualReader.getFileHeader();
        final VCFHeader expectedHeader = expectedReader.getFileHeader();
        Assert.assertEquals(actualHeader.getFilterLines().size(), expectedHeader.getFilterLines().size());
        Assert.assertEquals(actualHeader.getInfoHeaderLines().size(), expectedHeader.getInfoHeaderLines().size());
        Assert.assertEquals(actualHeader.getFormatHeaderLines().size(), expectedHeader.getFormatHeaderLines().size());

        // Check the number of records match, since we don't touch them
        Assert.assertEquals(actualReader.iterator().stream().count(), expectedReader.iterator().stream().count(), "The wrong number of variants was found.");

        CloserUtil.close(actualReader);
        CloserUtil.close(expectedReader);
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
                {-1,                HEADER_VCF,                   true},  // tests replacing with a new header, enforcing the same samples
                {-1,                HEADER_VCF,                   false}, // tests replacing with a new header, not enforcing the same samples
                {-1,                HEADER_VCF_WITH_EXTRA_SAMPLE, false}, // tests replacing with a new header with an extra sample, not enforcing the same samples
                {-1,                null,                         true},  // tests discovering the header lines from the VCF using all records
                {2,                 null,                         true},  // tests that the MissingFilter FILTER is found (in the second record)
                {Integer.MAX_VALUE, null,                         true}   // tests discovering the header lines from the VCF using all records (with CHECK_FIRST_N_RECORDS set to infinity)
        };
    }

    @Test(expectedExceptions=java.lang.IllegalStateException.class)
    public void testMissingFilterNotInFirstNRecords() throws IOException {
        // NB: the MissingFilter FILTER is in the second record, so not found in the first pass of the VCF
        runFixVcfHeader(1, null, true);
    }

    @Test(expectedExceptions=PicardException.class)
    public void testReplaceHeaderWithDifferentSamplesError() throws IOException {
        runFixVcfHeader(-1, HEADER_VCF_WITH_EXTRA_SAMPLE, true);
    }
}