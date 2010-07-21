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

import net.sf.samtools.*;
import org.testng.annotations.Test;
import net.sf.picard.io.IoUtil;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertNotNull;
import static org.testng.Assert.assertTrue;


/**
 * Test BAM file index creation
 */
public class BAMIndexWriterTest
{
    // Two input files for basic test
    private final String BAM_FILE_LOCATION = "testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam";
    private final String BAI_FILE_LOCATION = "testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam.bai";
    private final File BAM_FILE = new File(BAM_FILE_LOCATION);
    private final File BAI_FILE = new File(BAI_FILE_LOCATION);

    private final String LARGE_BAM_URL_STRING = "http://picard.sourceforge.net/testdata/test_human.bam";
    private final File BAM_Index_File = new File("testdata/net/sf/samtools/BAMFileIndexTest/test_human.bai");
    private final boolean mVerbose = true;

    @Test
    public void testWriteText() throws Exception {
        // Compare the text form of the c-generated bai file and a java-generated one
        final File cBaiTxtFile = File.createTempFile("cBai.", ".bai.txt");
        final BamIndexerForExistingBai bfi = new BamIndexerForExistingBai(BAI_FILE);
        bfi.createIndex(cBaiTxtFile, true, true);
        verbose ("Wrote textual C BAM Index file " + cBaiTxtFile);

        final File javaBaiFile = File.createTempFile("javaBai.", "java.bai");
        final File javaBaiTxtFile = new File(javaBaiFile.getAbsolutePath() + ".txt");
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        bam.enableFileSource(true);
        final BAMIndexer javaBai = new BAMIndexer(BAM_FILE, javaBaiFile,
                    bam.getFileHeader().getSequenceDictionary().size(), true);
        javaBai.createIndex();
        verbose ("Wrote binary Java BAM Index file " + javaBaiFile);
        
        // now, turn the bai file into text
        new BamIndexerForExistingBai(javaBaiFile).createIndex(javaBaiTxtFile, true, true);
        // and compare them
        verbose ("diff " + javaBaiTxtFile + " " + cBaiTxtFile);
        IoUtil.assertFilesEqual(javaBaiTxtFile, cBaiTxtFile);
        cBaiTxtFile.deleteOnExit();
        javaBaiFile.deleteOnExit();
        javaBaiTxtFile.deleteOnExit();
    }

    @Test
    public void testWriteBinary() throws Exception {
        // Compare java-generated bai file with c-generated and sorted bai file
        final File javaBaiFile = File.createTempFile("javaBai.", ".bai");
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        bam.enableFileSource(true);
        final BAMIndexer javaBai = new BAMIndexer(BAM_FILE, javaBaiFile,
                    bam.getFileHeader().getSequenceDictionary().size(), true);
        javaBai.createIndex();
        verbose ("Wrote binary java BAM Index file " + javaBaiFile);

        final File cRegeneratedBaiFile = File.createTempFile("cBai.", ".bai");
        final BamIndexerForExistingBai bfi2 = new BamIndexerForExistingBai(BAI_FILE);
        bfi2.createIndex(cRegeneratedBaiFile, false, true);
        verbose ("Wrote sorted C binary BAM Index file " + cRegeneratedBaiFile);

        // Binary compare of javaBaiFile and cRegeneratedBaiFile should be the same
        verbose ("diff " + javaBaiFile + " " + cRegeneratedBaiFile);
        IoUtil.assertFilesEqual(javaBaiFile, cRegeneratedBaiFile);
        javaBaiFile.deleteOnExit();
        cRegeneratedBaiFile.deleteOnExit();

    }

    @Test
    public void testBadIndex() throws Exception {
        
        final URL bamURL = new URL(LARGE_BAM_URL_STRING);
        final String sequence = "chr22";
        final int startPos = 624553;
        final int endPos = 5892114;
        int startWindow = LinearIndex.convertToLinearIndexOffset(startPos);

        verbose("Testing query " + sequence + ":" + startPos + "-" + endPos + " in window " + startWindow + " ...");

        final SAMFileReader reader1 = new SAMFileReader(bamURL, BAM_Index_File, false);
        final Iterator<SAMRecord> iter = reader1.queryOverlapping(sequence, startPos, endPos);
        int count = 0;
        while (iter.hasNext()){
            count ++;
        }
        verbose("Found " + count + " records in " + sequence + " at position " + startPos + " in window " + startWindow);
        assertEquals(count, 0);
    }

    private void verbose(final String text) {
        if (mVerbose) {
            System.out.println("#BAMIndexWriterTest " + text);
        }
    }
}
