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

import net.sf.picard.io.IoUtil;
import net.sf.picard.sam.BuildBamIndex;
import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Iterator;

import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

/**
 * Test BAM file index creation
 */
public class BAMIndexWriterTest {
    // Two input files for basic test
    private final String BAM_FILE_LOCATION = "testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam";
    private final String BAI_FILE_LOCATION = "testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam.bai";
    private final File BAM_FILE = new File(BAM_FILE_LOCATION);
    private final File BAI_FILE = new File(BAI_FILE_LOCATION);

    // another larger test case
    private final URL bamURL;
    private final String LARGE_BAM_URL_STRING = "http://picard.sourceforge.net/testdata/test_human.bam";
    private final File BAM_Index_File = new File("testdata/net/sf/samtools/BAMFileIndexTest/test_human.bai");
    
    private final boolean mVerbose = true;

    public BAMIndexWriterTest() throws Exception {
        bamURL = new URL(LARGE_BAM_URL_STRING);
    }

    @Test(enabled=true)
    public void testWriteText() throws Exception {
        // Compare the text form of the c-generated bai file and a java-generated one
        final File cBaiTxtFile = File.createTempFile("cBai.", ".bai.txt");
        BAMIndexer.createAndWriteIndex(BAI_FILE, cBaiTxtFile, true);
        verbose ("Wrote textual C BAM Index file " + cBaiTxtFile);

        final File javaBaiFile = File.createTempFile("javaBai.", "java.bai");
        final File javaBaiTxtFile = new File(javaBaiFile.getAbsolutePath() + ".txt");
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        BuildBamIndex.createIndex(bam, javaBaiFile);
        verbose ("Wrote binary Java BAM Index file " + javaBaiFile);
        
        // now, turn the bai file into text
        BAMIndexer.createAndWriteIndex(javaBaiFile, javaBaiTxtFile, true);
        // and compare them
        verbose ("diff " + javaBaiTxtFile + " " + cBaiTxtFile);
        IoUtil.assertFilesEqual(javaBaiTxtFile, cBaiTxtFile);
        cBaiTxtFile.deleteOnExit();
        javaBaiFile.deleteOnExit();
        javaBaiTxtFile.deleteOnExit();
    }

    @Test(enabled=true)
    public void testWriteBinary() throws Exception {
        // Compare java-generated bai file with c-generated and sorted bai file
        final File javaBaiFile = File.createTempFile("javaBai.", ".bai");
        final SAMFileReader bam = new SAMFileReader(BAM_FILE);
        BuildBamIndex.createIndex(bam, javaBaiFile);
        verbose ("Wrote binary java BAM Index file " + javaBaiFile);
                                                    
        final File cRegeneratedBaiFile = File.createTempFile("cBai.", ".bai");
        BAMIndexer.createAndWriteIndex(BAI_FILE, cRegeneratedBaiFile, false);
        verbose ("Wrote sorted C binary BAM Index file " + cRegeneratedBaiFile);

        // Binary compare of javaBaiFile and cRegeneratedBaiFile should be the same
        verbose ("diff " + javaBaiFile + " " + cRegeneratedBaiFile);
        IoUtil.assertFilesEqual(javaBaiFile, cRegeneratedBaiFile);
        javaBaiFile.deleteOnExit();
        cRegeneratedBaiFile.deleteOnExit();

    }

    @Test(enabled=true)
    public void testBadIndex() throws Exception {
        // This test case is one that fails when bin meta data is incorrectly used as data

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

    @Test(enabled = false, dataProvider = "linearIndexTestData")
    /** Test linear index at specific references and windows */
    public void testLinearIndex(String testName, String filepath, int problemReference, int problemWindowStart, int problemWindowEnd, int expectedCount) {
        final SAMFileReader sfr = new SAMFileReader(new File(filepath));
        for (int problemWindow = problemWindowStart; problemWindow <= problemWindowEnd; problemWindow++) {
            int count = countAlignmentsInWindow(problemReference, problemWindow, sfr, expectedCount);
            if (expectedCount != -1)
                assertEquals(expectedCount, count);
        }
    }

    @DataProvider(name = "linearIndexTestData")
    public Object[][] getLinearIndexTestData() {
        // Add data here for test cases, reference, and windows where linear index needs testing
        return new Object[][]{
                new Object[]{"index_test", BAM_FILE_LOCATION, 1, 29, 66, -1},  // 29-66
                new Object[]{"index_test", BAM_FILE_LOCATION, 1, 68, 118, -1},  // 29-66

        };
    }

    private int countAlignmentsInWindow(int reference, int window, SAMFileReader reader, int expectedCount) {
        final int SIXTEEN_K = 1 << 14;       // 1 << LinearIndex.BAM_LIDX_SHIFT
        final int start = window >> 14;             // window * SIXTEEN_K;
        final int stop = ((window + 1) >> 14) - 1; // (window + 1 * SIXTEEN_K) - 1;

        final String chr = reader.getFileHeader().getSequence(reference).getSequenceName();

        // get records for the entire linear index window
        SAMRecordIterator iter = reader.queryOverlapping(chr, start, stop);
        SAMRecord rec;
        int count = 0;
        while (iter.hasNext()) {
            rec = iter.next();
            count++;
            if (expectedCount == -1)
                System.err.println(rec.getReadName());
        }
        iter.close();
        return count;
    }


    @Test(enabled=false, dataProvider = "indexComparisonData")
    /** Test linear index at all references and windows, comparing with existing index */
    public void compareLinearIndex(String testName, String bamFile, String bamIndexFile) throws IOException{
        // compare index generated from bamFile with existing bamIndex file
        // by testing all the references' windows and comparing the counts

        // 1. generate bai file
        // 2. count its references
        // 3. count bamIndex references comparing counts

        // 1. generate bai file
        File bam = new File(bamFile);
        assertTrue(bam.exists(), testName + " input bam file doesn't exist: " + bamFile);

        File indexFile1 = createIndexFile(bam);
        assertTrue(indexFile1.exists(), testName + " generated bam file's index doesn't exist: " + indexFile1);

        // 2. count its references
        File indexFile2 = new File (bamIndexFile);
        assertTrue(indexFile2.exists(), testName + " input index file doesn't exist: " + indexFile2);

        final CachingBAMFileIndex existingIndex1 = new CachingBAMFileIndex(indexFile1, null); // todo null sequence dictionary?
        final CachingBAMFileIndex existingIndex2 = new CachingBAMFileIndex(indexFile2, null);
        final int n_ref = existingIndex1.getNumberOfReferences();
        assertEquals(n_ref, existingIndex2.getNumberOfReferences());

        final SAMFileReader reader1 = new SAMFileReader(bam, indexFile1, false);
        final SAMFileReader reader2 = new SAMFileReader(bam, indexFile2, false );

        System.out.println("Comparing " + n_ref + " references in " + indexFile1 + " and " + indexFile2);

        for (int i = 0; i < n_ref; i++) {
            final BAMIndexContent content1 = existingIndex1.getQueryResults(i);
            final BAMIndexContent content2 = existingIndex2.getQueryResults(i);
            if (content1 == null){
                assertTrue(content2 == null, "No content for 1st bam index, but content for second at reference" + i);
                continue;
            }
            int[] counts1 = new int[LinearIndex.MAX_LINEAR_INDEX_SIZE];
            int[] counts2 = new int[LinearIndex.MAX_LINEAR_INDEX_SIZE];
            LinearIndex li1 = content1.getLinearIndex();
            LinearIndex li2 = content2.getLinearIndex();
            // todo not li1 and li2 sizes may differ. Implies 0's in the smaller index windows
            // 3. count bamIndex references comparing counts
            int baiSize = Math.max(li1.size(), li2.size());
            for (int win = 0; win < baiSize; win++) {
                counts1[win] = countAlignmentsInWindow(i, win, reader1, 0);
                counts2[win] = countAlignmentsInWindow(i, win, reader2, counts1[win]);
                assertEquals(counts2[win], counts1[win], "Counts don't match for reference " + i  +
                        " window " + win );
            }
        }

        indexFile1.deleteOnExit();

    }

    @DataProvider(name="indexComparisonData")
    public Object[][] getIndexComparisonData() {
        // enter bam file and alternate index file to be tested against generated bam index
        return new Object[][]{
           new Object[]{"index_test", BAM_FILE_LOCATION, BAI_FILE_LOCATION} ,
        };
    }

    /**  generates the index file using the latest java index generating code  */
    private File createIndexFile(File bamFile) throws IOException {
        final File bamIndexFile = File.createTempFile("Bai.", ".bai");
        final SAMFileReader bam = new SAMFileReader(bamFile);
        BuildBamIndex.createIndex(bam, bamIndexFile);
        verbose ("Wrote BAM Index file " + bamIndexFile);
        return bamIndexFile;
    }

    private void verbose(final String text) {
        if (mVerbose) {
            System.out.println("#BAMIndexWriterTest " + text);
        }
    }
}
