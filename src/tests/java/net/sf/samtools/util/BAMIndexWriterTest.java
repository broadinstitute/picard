/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.samtools.util;

import net.sf.samtools.SAMFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import net.sf.samtools.BAMFileIndexWriter;
import net.sf.samtools.BAMIndexTextWriter;
import net.sf.picard.io.IoUtil;

import java.io.File;

import static org.testng.Assert.*;

/**
 * Test BAM file index creation
 */
public class BAMIndexWriterTest
{
    // Two input files
    private final String BAM_FILE_LOCATION = "testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam";
    private final String BAI_FILE_LOCATION = "testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam.bai";
    private final File BAM_FILE = new File(BAM_FILE_LOCATION);
    private final File BAI_FILE = new File(BAI_FILE_LOCATION);

    private final boolean mVerbose = false;

    @Test
    public void testWriteText() throws Exception {
        // Compare the text form of the c-generated bai file and a java-generated one
        // final String cBaiTxtFileName = BAM_FILE_LOCATION + ".bai.txt";
        final File cBaiTxtFile = File.createTempFile("cBai.", ".bai.txt");
        cBaiTxtFile.deleteOnExit();
        final BAMIndexTextWriter bfi = new BAMIndexTextWriter(BAI_FILE, cBaiTxtFile);
        bfi.writeText(true);
        verbose ("Wrote Textual BAM Index file " + cBaiTxtFile);

        // Text compare of javaBaiTxtFileName.txt(bfi2)
        //             and cBaiTxtFileName(bfi)should be the same
        //final String javaBaiTxtFileName = BAM_FILE_LOCATION + ".java.bai.txt";  // java-generated
        final File javaBaiTxtFile = File.createTempFile("javaBai.", "java.bai.txt");
        javaBaiTxtFile.deleteOnExit();
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        bam.enableFileSource(true);
        final BAMFileIndexWriter javaBai = new BAMFileIndexWriter(javaBaiTxtFile,
                    bam.getFileHeader().getSequenceDictionary().size());
        int n_records = javaBai.createIndex(BAM_FILE, true, true);
        Assert.assertEquals(9721, n_records);
        // diff index_test.bam.java.bai.txt index_text.bam.bai.txt
        verbose ("diff " + javaBaiTxtFile + " " + cBaiTxtFile);
        IoUtil.assertFilesEqual(javaBaiTxtFile, cBaiTxtFile);
    }

    @Test
    public void testWriteBinary() throws Exception {
        // Compare c-generated and sorted bai file with a java-generated bai file
        // final String javaBaiFileName = BAM_FILE.getPath() + ".java.bai";  // java-generated
        final File javaBaiFile = File.createTempFile("javaBai.", ".bai");
        javaBaiFile.deleteOnExit();
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        bam.enableFileSource(true);
        final BAMFileIndexWriter javaBai = new BAMFileIndexWriter(javaBaiFile,
                    bam.getFileHeader().getSequenceDictionary().size());
        int n_records = javaBai.createIndex(BAM_FILE, false, true);
        Assert.assertEquals(9721, n_records);
        verbose ("Wrote Binary BAM Index file " + javaBaiFile);

        // final String cRegeneratedBaiFileName = BAM_FILE.getPath() + ".generated.bai";  // java-generated
        final File cRegeneratedBaiFile = File.createTempFile("cBai.", ".bai.txt");
        cRegeneratedBaiFile.deleteOnExit();
        final BAMIndexTextWriter bfi2 = new BAMIndexTextWriter(BAI_FILE, cRegeneratedBaiFile);
        bfi2.writeBinary(true, 0);
        verbose ("Wrote Binary BAM Index file " + javaBaiFile);
        // Binary compare of javaBaiFileName and cBaiFileName.sorted should be the same
        // diff index_test.bam.java.bai index_test.bam.generated.bai
        verbose ("diff " + javaBaiFile + " " + cRegeneratedBaiFile);
        IoUtil.assertFilesEqual(javaBaiFile, cRegeneratedBaiFile);

        // todo - write both of these as text format and compare
    }

    private void verbose(final String text) {
        if (mVerbose) {
            System.out.println("# " + text);
        }
    }
}
