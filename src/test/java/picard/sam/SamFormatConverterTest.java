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

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.util.SamComparison;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class SamFormatConverterTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/SamFormatConverterTest");
    private static final File unmappedSam = new File(TEST_DATA_DIR, "unmapped.sam");
    private static final File unmappedBam = new File(TEST_DATA_DIR, "unmapped.bam");
    private static final File unmappedCram = new File(TEST_DATA_DIR, "unmapped.cram");

    private static final File ESSENTIALLY_EMPTY_REFERENCE_TO_USE_WITH_UNMAPPED_CRAM = new File(TEST_DATA_DIR, "basicallyEmpty.fasta");

    @DataProvider
    public Object[][] conversionCases() {
        return new Object[][]{
                {unmappedSam, unmappedBam, ".bam"},
                {unmappedSam, unmappedCram, ".cram"},
                {unmappedBam, unmappedCram, ".cram"},
                {unmappedBam, unmappedSam, ".sam"},
                {unmappedCram, unmappedBam, ".bam"},
                {unmappedCram, unmappedSam, ".sam"},

        };
    }

    @Test(dataProvider = "conversionCases")
    public void testConvert(File input, File expected, String extension) throws IOException {
        convertFile(input, expected, extension);
    }


    private void convertFile(final File inputFile, final File fileToCompare, final String extension) throws IOException {
        final List<File> samFiles = new ArrayList<>();
        final ValidateSamFile validateSamFile = new ValidateSamFile();

        final File output = File.createTempFile("SamFormatConverterTest." + inputFile.getName(), extension);
        output.deleteOnExit();
        SamFormatConverter.convert(inputFile, output, ESSENTIALLY_EMPTY_REFERENCE_TO_USE_WITH_UNMAPPED_CRAM, Defaults.CREATE_INDEX);

        validateSamFile.INPUT = output;
        assertEquals(validateSamFile.doWork(), 0);
        final SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().referenceSequence(ESSENTIALLY_EMPTY_REFERENCE_TO_USE_WITH_UNMAPPED_CRAM);
        try (final SamReader samReader1 = samReaderFactory.open(output);
             final SamReader samReader2 = samReaderFactory.open(fileToCompare)) {
            final SamComparison samComparison = new SamComparison(samReader1, samReader2);
            Assert.assertTrue(samComparison.areEqual());
        }

    }
}
