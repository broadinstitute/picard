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

import net.sf.samtools.BAMIndexer;
import net.sf.samtools.BamIndexerForExistingBai;
import net.sf.samtools.SAMFileReader;
import org.testng.annotations.Test;
import net.sf.picard.io.IoUtil;

import java.io.File;


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

    private final boolean mVerbose = true;

    @Test
    public void testWriteText() throws Exception {
        // Compare the text form of the c-generated bai file and a java-generated one
        final File cBaiTxtFile = File.createTempFile("cBai.", ".bai.txt");
        final BamIndexerForExistingBai bfi = new BamIndexerForExistingBai(BAI_FILE, true);
        bfi.createIndex(cBaiTxtFile, true);
        verbose ("Wrote Textual version of C BAM Index file " + cBaiTxtFile);

        // Text compare of javaBaiTxtFileName.txt(bfi2)
        //             and cBaiTxtFileName(bfi)should be the same
        final File javaBaiTxtFile = File.createTempFile("javaBai.", "java.bai.txt");
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        bam.enableFileSource(true);
        final BAMIndexer javaBai = new BAMIndexer(BAM_FILE, javaBaiTxtFile,
                    bam.getFileHeader().getSequenceDictionary().size(), true, true);
        javaBai.createIndex();
        verbose ("Wrote Textual version of Java BAM Index file " + javaBai);
        // diff index_test.bam.java.bai.txt index_text.bam.bai.txt
        verbose ("diff " + javaBaiTxtFile + " " + cBaiTxtFile);
        IoUtil.assertFilesEqual(javaBaiTxtFile, cBaiTxtFile);
        cBaiTxtFile.deleteOnExit();
        javaBaiTxtFile.deleteOnExit();
    }

    @Test
    public void testWriteBinary() throws Exception {
        // Compare java-generated bai file with c-generated and sorted bai file
        final File javaBaiFile = File.createTempFile("javaBai.", ".bai");
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        bam.enableFileSource(true);
        final BAMIndexer javaBai = new BAMIndexer(BAM_FILE, javaBaiFile,
                    bam.getFileHeader().getSequenceDictionary().size(), true, false);
        javaBai.createIndex();
        verbose ("Wrote sorted C Binary BAM Index file " + javaBaiFile);

        final File cRegeneratedBaiFile = File.createTempFile("cBai.", ".bai");
        final BamIndexerForExistingBai bfi2 = new BamIndexerForExistingBai(BAI_FILE, true);
        bfi2.createIndex(cRegeneratedBaiFile, false);
        verbose ("Wrote Java-generated Binary BAM Index file " + cRegeneratedBaiFile);
        // Binary compare of javaBaiFile and cRegeneratedBaiFile should be the same
        // diff index_test.bam.java.bai index_test.bam.generated.bai
        verbose ("diff " + javaBaiFile + " " + cRegeneratedBaiFile);
        IoUtil.assertFilesEqual(javaBaiFile, cRegeneratedBaiFile);
        javaBaiFile.deleteOnExit();
        cRegeneratedBaiFile.deleteOnExit();

    }

    private void verbose(final String text) {
        if (mVerbose) {
            System.out.println("# " + text);
        }
    }
}
