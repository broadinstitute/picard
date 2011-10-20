/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.samtools;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;


import java.io.File;

public class SAMFileReaderTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/samtools");

    @Test(dataProvider = "variousFormatReaderTestCases")
    public void variousFormatReaderTest(final String inputFile) {
        final File input = new File(TEST_DATA_DIR, inputFile);
        final SAMFileReader reader = new SAMFileReader(input);
        for (final SAMRecord rec: reader) {}
        reader.close();
    }

    @DataProvider(name = "variousFormatReaderTestCases")
    public Object[][] variousFormatReaderTestCases() {
        final Object[][] scenarios = new Object[][] {
                { "block_compressed.sam.gz"},
                { "uncompressed.sam"},
                { "compressed.sam.gz"},
                { "compressed.bam"},
        };
        return scenarios;
    }

    // Tests for the SAMRecordFactory usage
    class SAMRecordFactoryTester extends DefaultSAMRecordFactory {
        int samRecordsCreated;
        int bamRecordsCreated;

        public SAMRecord createSAMRecord(final SAMFileHeader header) {
            ++samRecordsCreated;
            return super.createSAMRecord(header);
        }

        public BAMRecord createBAMRecord(final SAMFileHeader header, final int referenceSequenceIndex, final int alignmentStart, final short readNameLength, final short mappingQuality, final int indexingBin, final int cigarLen, final int flags, final int readLen, final int mateReferenceSequenceIndex, final int mateAlignmentStart, final int insertSize, final byte[] variableLengthBlock) {
            ++bamRecordsCreated;
            return super.createBAMRecord(header, referenceSequenceIndex, alignmentStart, readNameLength, mappingQuality, indexingBin, cigarLen, flags, readLen, mateReferenceSequenceIndex, mateAlignmentStart, insertSize, variableLengthBlock);
        }
    }

    @Test(dataProvider = "variousFormatReaderTestCases")
    public void samRecordFactoryTest(final String inputFile) {
        final File input = new File(TEST_DATA_DIR, inputFile);
        final SAMFileReader reader = new SAMFileReader(input);
        final SAMRecordFactoryTester factory = new SAMRecordFactoryTester();
        reader.setSAMRecordFactory(factory);

        int i=0;
        for (final SAMRecord rec: reader) {
            ++i;
        }
        reader.close();

        Assert.assertTrue(i > 0);
        if (inputFile.endsWith(".sam") || inputFile.endsWith(".sam.gz")) Assert.assertEquals(factory.samRecordsCreated, i);
        else if (inputFile.endsWith(".bam")) Assert.assertEquals(factory.bamRecordsCreated, i);
    }

}
