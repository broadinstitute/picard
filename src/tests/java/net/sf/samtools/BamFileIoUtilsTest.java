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

package net.sf.samtools;

import net.sf.picard.PicardException;
import net.sf.picard.sam.CompareSAMs;
import net.sf.samtools.BamFileIoUtils;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class BamFileIoUtilsTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam");
    private static final File HEADER_PROVIDER = new File(TEST_DATA_DIR, "samHeaderProvider.sam");
    private static final File INPUT_FILE = new File(TEST_DATA_DIR, "aligned_queryname_sorted.bam");

    @Test
    public void testReheader() throws Exception {
        final SAMFileHeader samFileHeader = new SAMFileReader(HEADER_PROVIDER).getFileHeader();
        final File outputFile = createOutputFile();
        BamFileIoUtils.reheaderBamFile(samFileHeader, INPUT_FILE, outputFile);
        assertReheader(outputFile, HEADER_PROVIDER);
    }

    @Test
    public void testReapplySameHeader() throws Exception {
        final SAMFileHeader samFileHeader = new SAMFileReader(INPUT_FILE).getFileHeader();
        final File outputFile = createOutputFile();
        BamFileIoUtils.reheaderBamFile(samFileHeader, INPUT_FILE, outputFile);

        assertReheader(outputFile, INPUT_FILE);
        final CompareSAMs compareSAMs = new CompareSAMs();
        compareSAMs.instanceMain(new String[]{INPUT_FILE.getAbsolutePath(), outputFile.getAbsolutePath()});
        Assert.assertTrue(compareSAMs.areEqual());
    }

    @Test(expectedExceptions = SAMException.class)
    public void mismatchSortOrderTest() throws Exception {
        final SAMFileHeader samFileHeader = new SAMFileReader(HEADER_PROVIDER).getFileHeader();
        samFileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        BamFileIoUtils.reheaderBamFile(samFileHeader, INPUT_FILE, createOutputFile());
        throw new IllegalStateException("We shouldn't be here");
    }

    private void assertReheader(final File outputFile, final File headerProvider) {
        final SAMFileHeader origHeader = new SAMFileReader(headerProvider).getFileHeader();
        final SAMFileHeader newHeader = new SAMFileReader(outputFile).getFileHeader();
        Assert.assertEquals(origHeader, newHeader);
    }

    private File createOutputFile() throws IOException {
        return File.createTempFile("reheaderTest.", BamFileIoUtils.BAM_FILE_EXTENSION);
    }

}
