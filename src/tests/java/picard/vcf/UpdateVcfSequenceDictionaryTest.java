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
package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;
import picard.PicardException;

import java.io.File;

/**
 * @author George Grant
 */
public class UpdateVcfSequenceDictionaryTest {
    private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");
    private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("UpdateVcfSequenceDictionaryTest", null);

    @AfterClass
    public void teardown() {
        IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
    }

    @Test
    public void testUpdateVcfSequenceDictionary() {
        final File input = new File(TEST_DATA_PATH, "vcfFormatTest.vcf");
        // vcfFormatTest.bad_dict.vcf is a vcf with two (2) ##contig lines deleted
        final File samSequenceDictionaryVcf = new File(TEST_DATA_PATH, "vcfFormatTest.bad_dict.vcf");
        final File outputFile = new File(OUTPUT_DATA_PATH, "updateVcfSequenceDictionaryTest-delete-me.vcf");

        outputFile.deleteOnExit();

        final UpdateVcfSequenceDictionary updateVcfSequenceDictionary = new UpdateVcfSequenceDictionary();
        updateVcfSequenceDictionary.INPUT = input;
        updateVcfSequenceDictionary.SEQUENCE_DICTIONARY = samSequenceDictionaryVcf;
        updateVcfSequenceDictionary.OUTPUT = outputFile;

        Assert.assertEquals(updateVcfSequenceDictionary.instanceMain(new String[0]), 0);

        IOUtil.assertFilesEqual(samSequenceDictionaryVcf, outputFile);

        // A little extra checking.
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(input).size(), 84);
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(samSequenceDictionaryVcf).size(), 82);
        Assert.assertEquals(SAMSequenceDictionaryExtractor.extractDictionary(outputFile).size(), 82);
    }
}
