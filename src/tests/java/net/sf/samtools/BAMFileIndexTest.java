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
package net.sf.samtools;

import org.testng.annotations.Test;
import static org.testng.Assert.*;

import java.io.*;
import java.util.*;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.BAMFileIndex;

/**
 * Test BAM file indexing.
 */
public class BAMFileIndexTest
{
    private File BAM_FILE = new File("etc/testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam");
    private boolean mVerbose = false;

    @Test
    public void testGetSearchBins()
        throws Exception {
        BAMFileIndex bfi = new BAMFileIndex(new File(BAM_FILE.getPath() + ".bai"));
        long[] bins = bfi.getSearchBins(1, 0, 0);
        /***
        if (bins == null) {
            System.out.println("Search bins: " + bins);
            return;
        }
        System.out.println("Search bins:");
        for (int i = 0; i < bins.length; i++) {
            System.out.println(" " + Long.toHexString(bins[i]));
        }
        ***/
        assertNotNull(bins);
        assertEquals(bins.length, 2);
    }

    @Test
    public void testSpecificQueries()
        throws Exception {
        assertEquals(runQueryTest(BAM_FILE, "chrM", 10400, 10600, true), 1);
        assertEquals(runQueryTest(BAM_FILE, "chrM", 10400, 10600, false), 2);
    }

    @Test(enabled = false)
    public void testRandomQueries()
        throws Exception {
        runRandomTest(BAM_FILE, 1000, new Random());
    }

    @Test
    public void testWholeChromosomes() {
        checkChromosome("chrM", 23);
        checkChromosome("chr1", 885);
        checkChromosome("chr2", 837);
        /***
        checkChromosome("chr3", 683);
        checkChromosome("chr4", 633);
        checkChromosome("chr5", 611);
        checkChromosome("chr6", 585);
        checkChromosome("chr7", 521);
        checkChromosome("chr8", 507);
        checkChromosome("chr9", 388);
        checkChromosome("chr10", 477);
        checkChromosome("chr11", 467);
        checkChromosome("chr12", 459);
        checkChromosome("chr13", 327);
        checkChromosome("chr14", 310);
        checkChromosome("chr15", 280);
        checkChromosome("chr16", 278);
        checkChromosome("chr17", 269);
        checkChromosome("chr18", 265);
        checkChromosome("chr19", 178);
        checkChromosome("chr20", 228);
        checkChromosome("chr21", 123);
        checkChromosome("chr22", 121);
        checkChromosome("chrX", 237);
        checkChromosome("chrY", 29);
        ***/
    }

    private void checkChromosome(String name, int expectedCount) {
        int count = runQueryTest(BAM_FILE, name, 0, 0, true);
        assertEquals(count, expectedCount);
        count = runQueryTest(BAM_FILE, name, 0, 0, false);
        assertEquals(count, expectedCount);
    }

    private void runRandomTest(File bamFile, int count, Random generator) {
        int maxCoordinate = 10000000;
        List<String> referenceNames = getReferenceNames(bamFile);
        for (int i = 0; i < count; i++) {
            String refName = referenceNames.get(generator.nextInt(referenceNames.size()));
            int coord1 = generator.nextInt(maxCoordinate+1);
            int coord2 = generator.nextInt(maxCoordinate+1);
            int startPos = Math.min(coord1, coord2);
            int endPos = Math.max(coord1, coord2);
            System.out.println("Testing query " + refName + ":" + startPos + "-" + endPos + " ...");
            try {
                runQueryTest(bamFile, refName, startPos, endPos, true);
                runQueryTest(bamFile, refName, startPos, endPos, false);
            } catch (Throwable exc) {
                String message = "Query test failed: " + refName + ":" + startPos + "-" + endPos;
                message += ": " + exc.getMessage();
                throw new RuntimeException(message, exc);
            }
        }
    }

    private List<String> getReferenceNames(File bamFile) {
        SAMFileReader reader = new SAMFileReader(bamFile);
        List<String> result = new ArrayList<String>();
        List<SAMSequenceRecord> seqRecords = reader.getFileHeader().getSequenceDictionary().getSequences();
        for (SAMSequenceRecord seqRecord : seqRecords) {
            if (seqRecord.getSequenceName() != null) {
                result.add(seqRecord.getSequenceName());
            }
        }
        reader.close();
        return result;
    }

    private int runQueryTest(File bamFile, String sequence, int startPos, int endPos, boolean contained) {
        verbose("Testing query " + sequence + ":" + startPos + "-" + endPos + " ...");
        SAMFileReader reader1 = new SAMFileReader(bamFile);
        SAMFileReader reader2 = new SAMFileReader(bamFile);
        Iterator<SAMRecord> iter1 = reader1.query(sequence, startPos, endPos, contained);
        Iterator<SAMRecord> iter2 = reader2.iterator();
        // Compare ordered iterators.
        // Confirm that iter1 is a subset of iter2 that properly filters.
        SAMRecord record1 = null;
        SAMRecord record2 = null;
        int count1 = 0;
        int count2 = 0;
        int beforeCount = 0;
        int afterCount = 0;
        while (true) {
            if (record1 == null && iter1.hasNext()) {
                record1 = iter1.next();
                count1++;
            }
            if (record2 == null && iter2.hasNext()) {
                record2 = iter2.next();
                count2++;
            }
            // System.out.println("Iteration:");
            // System.out.println(" Record1 = " + ((record1 == null) ? "null" : record1.format()));
            // System.out.println(" Record2 = " + ((record2 == null) ? "null" : record2.format()));
            if (record1 == null && record2 == null) {
                break;
            }
            if (record1 == null) {
                checkPassesFilter(false, record2, sequence, startPos, endPos, contained);
                record2 = null;
                afterCount++;
                continue;
            }
            assertNotNull(record2);
            int ordering = compareCoordinates(record1, record2);
            if (ordering > 0) {
                checkPassesFilter(false, record2, sequence, startPos, endPos, contained);
                record2 = null;
                beforeCount++;
                continue;
            }
            assertTrue(ordering == 0);
            checkPassesFilter(true, record1, sequence, startPos, endPos, contained);
            checkPassesFilter(true, record2, sequence, startPos, endPos, contained);
            assertEquals(record1.getReadName(), record2.getReadName());
            assertEquals(record1.getReadString(), record2.getReadString());
            record1 = null;
            record2 = null;
        }
        reader1.close();
        reader2.close();
        verbose("Checked " + count1 + " records against " + count2 + " records.");
        verbose("Found " + (count2 - beforeCount - afterCount) + " records matching.");
        verbose("Found " + beforeCount + " records before.");
        verbose("Found " + afterCount + " records after.");
        return count1;
    }

    private void checkPassesFilter(boolean expected, SAMRecord record, String sequence, int startPos, int endPos, boolean contained) {
        boolean passes = passesFilter(record, sequence, startPos, endPos, contained);
        if (passes != expected) {
            System.out.println("Error: Record erroneously " +
                               (passes ? "passed" : "failed") +
                               " filter.");
            System.out.println(" Record: " + record.format());
            System.out.println(" Filter: " + sequence + ":" +
                               startPos + "-" + endPos +
                               " (" + (contained ? "contained" : "overlapping") + ")");
            assertEquals(passes, expected);
        }
    }

    private boolean passesFilter(SAMRecord record, String sequence, int startPos, int endPos, boolean contained) {
        if (record == null) {
            return false;
        }
        if (!safeEquals(record.getReferenceName(), sequence)) {
            return false;
        }
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();
        if (alignmentStart <= 0) {
            assertTrue(record.getReadUnmappedFlag());
            return false;
        }
        if (alignmentEnd <= 0) {
            // For indexing-only records, treat as single base alignment.
            assertTrue(record.getReadUnmappedFlag());
            alignmentEnd = alignmentStart;
        }
        if (contained) {
            if (startPos != 0 && alignmentStart < startPos) {
                return false;
            }
            if (endPos != 0 && alignmentEnd > endPos) {
                return false;
            }
        } else {
            if (startPos != 0 && alignmentEnd < startPos) {
                return false;
            }
            if (endPos != 0 && alignmentStart > endPos) {
                return false;
            }
        }
        return true;
    }

    private int compareCoordinates(SAMRecord record1, SAMRecord record2) {
        int seqIndex1 = record1.getReferenceIndex();
        int seqIndex2 = record2.getReferenceIndex();
        if (seqIndex1 == -1) {
            return ((seqIndex2 == -1) ? 0 : -1);
        } else if (seqIndex2 == -1) {
            return 1;
        }
        int result = seqIndex1 - seqIndex2;
        if (result != 0) {
            return result;
        }
        result = record1.getAlignmentStart() - record2.getAlignmentStart();
        return result;
    }

    private boolean safeEquals(Object o1, Object o2) {
        if (o1 == o2) {
            return true;
        } else if (o1 == null || o2 == null) {
            return false;
        } else {
            return o1.equals(o2);
        }
    }

    private void verbose(String text) {
        if (mVerbose) {
            System.out.println("# " + text);
        }
    }
}
