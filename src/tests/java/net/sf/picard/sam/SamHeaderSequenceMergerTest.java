/**
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
 **/


package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import static org.testng.Assert.assertEquals;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author aaron
 * @version 1.0
 * @date May 20, 2009
 * <p/>
 * Class SamFileHeaderMergerTest
 * <p/>
 * Tests the ability of the SamFileHeaderMerger class to merge sequence dictionaries.
 */
public class SamHeaderSequenceMergerTest {
    private static File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam");

    /** tests that if we've set the merging to false, we get a PicardException for bam's with different dictionaries. */
    @Test(expectedExceptions = PicardException.class)
    public void testMergedException() {
        File INPUT[] = {new File(TEST_DATA_DIR, "Chromosome1to10.bam"),
                        new File(TEST_DATA_DIR, "Chromosome5to9.bam")};
        final List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        for (final File inFile : INPUT) {
            IoUtil.assertFileIsReadable(inFile);
            final SAMFileReader in = new SAMFileReader(inFile);
            readers.add(in);
        }
        final MergingSamRecordIterator iterator;
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SAMFileHeader.SortOrder.unsorted, false);

    }

    /** Tests that we can successfully merge two files with */
    @Test
    public void testMerging() {
        File INPUT[] = {new File(TEST_DATA_DIR, "Chromosome1to10.bam"),
                        new File(TEST_DATA_DIR, "Chromosome5to9.bam")};
        final List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        for (final File inFile : INPUT) {
            IoUtil.assertFileIsReadable(inFile);
            final SAMFileReader in = new SAMFileReader(inFile);
            readers.add(in);
        }
        final MergingSamRecordIterator iterator;
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SAMFileHeader.SortOrder.unsorted, true);
        iterator = new MergingSamRecordIterator(headerMerger, false);
        final SAMFileHeader header = headerMerger.getMergedHeader();

        // count the total reads, and record read counts for each sequence
        Map<Integer, Integer> seqCounts = new HashMap<Integer, Integer>();
        int totalCount = 0;

        while (iterator.hasNext()) {
            SAMRecord r = iterator.next();
            if (seqCounts.containsKey(r.getReferenceIndex())) {
                seqCounts.put(r.getReferenceIndex(), seqCounts.get(r.getReferenceIndex()) + 1);
            } else {
                seqCounts.put(r.getReferenceIndex(), 1);
            }
            ++totalCount;
        }
        assertEquals(totalCount, 1500);
        for (Integer i : seqCounts.keySet()) {
            if (i < 4 || i > 8) {
                // seqeunce 5 - 9 should have 200 reads (indices 4 - 8)
                assertEquals(seqCounts.get(i).intValue(), 100);
            } else {
                // the others should have 100
                assertEquals(seqCounts.get(i).intValue(), 200);
            }
        }
    }

}
