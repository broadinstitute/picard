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

import static org.testng.Assert.*;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

/**
 * Test BAM file indexing.
 */
public class BAMRemoteFileTest {
    private final File BAM_INDEX_FILE = new File("testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam.bai");
    private final File BAM_FILE = new File("testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam");
    private final String BAM_URL_STRING = "http://picard.sourceforge.net/testdata/index_test.bam";
    private final URL bamURL;

    private final boolean mVerbose = false;

    public BAMRemoteFileTest() throws Exception {
        bamURL = new URL(BAM_URL_STRING);
    }


    @Test
    void testRemoteLocal()
            throws Exception {
        runLocalRemoteTest(bamURL, BAM_FILE, "chrM", 10400, 10600, false);
    }

    @Test
    public void testSpecificQueries()
            throws Exception {
        assertEquals(runQueryTest(bamURL, "chrM", 10400, 10600, true), 1);
        assertEquals(runQueryTest(bamURL, "chrM", 10400, 10600, false), 2);
    }

    @Test(enabled = true)
    public void testRandomQueries()
            throws Exception {
        runRandomTest(bamURL, 20, new Random());
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


    private void checkChromosome(final String name, final int expectedCount) {
        int count = runQueryTest(bamURL, name, 0, 0, true);
        assertEquals(count, expectedCount);
        count = runQueryTest(bamURL, name, 0, 0, false);
        assertEquals(count, expectedCount);
    }

    private void runRandomTest(final URL bamFile, final int count, final Random generator) throws IOException {
        final int maxCoordinate = 10000000;
        final List<String> referenceNames = getReferenceNames(bamFile);
        for (int i = 0; i < count; i++) {
            final String refName = referenceNames.get(generator.nextInt(referenceNames.size()));
            final int coord1 = generator.nextInt(maxCoordinate + 1);
            final int coord2 = generator.nextInt(maxCoordinate + 1);
            final int startPos = Math.min(coord1, coord2);
            final int endPos = Math.max(coord1, coord2);
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

    private List<String> getReferenceNames(final URL bamFile) throws IOException {


        final SAMFileReader reader = new SAMFileReader(bamFile.openStream());

        final List<String> result = new ArrayList<String>();
        final List<SAMSequenceRecord> seqRecords = reader.getFileHeader().getSequenceDictionary().getSequences();
        for (final SAMSequenceRecord seqRecord : seqRecords) {
            if (seqRecord.getSequenceName() != null) {
                result.add(seqRecord.getSequenceName());
            }
        }
        reader.close();
        return result;
    }

    private void runLocalRemoteTest(final URL bamURL, final File bamFile, final String sequence, final int startPos, final int endPos, final boolean contained) {
        verbose("Testing query " + sequence + ":" + startPos + "-" + endPos + " ...");
        final SAMFileReader reader1 = new SAMFileReader(bamFile, BAM_INDEX_FILE, false);
        final SAMFileReader reader2 = new SAMFileReader(bamURL, BAM_INDEX_FILE, false);
        final Iterator<SAMRecord> iter1 = reader1.query(sequence, startPos, endPos, contained);
        final Iterator<SAMRecord> iter2 = reader2.query(sequence, startPos, endPos, contained);

        final List<SAMRecord> records1 = new ArrayList<SAMRecord>();
        final List<SAMRecord> records2 = new ArrayList<SAMRecord>();

        while (iter1.hasNext()) {
            records1.add(iter1.next());
        }
        while (iter2.hasNext()) {
            records2.add(iter2.next());
        }

        assertTrue(records1.size() > 0);
        assertEquals(records1.size(), records2.size());
        for (int i = 0; i < records1.size(); i++) {
            //System.out.println(records1.get(i).format());
            assertEquals(records1.get(i).format(), records2.get(i).format());
        }


    }

    private int runQueryTest(final URL bamURL, final String sequence, final int startPos, final int endPos, final boolean contained) {
        verbose("Testing query " + sequence + ":" + startPos + "-" + endPos + " ...");
        final SAMFileReader reader1 = new SAMFileReader(bamURL, BAM_INDEX_FILE, false);
        final SAMFileReader reader2 = new SAMFileReader(bamURL, BAM_INDEX_FILE, false);
        final Iterator<SAMRecord> iter1 = reader1.query(sequence, startPos, endPos, contained);
        final Iterator<SAMRecord> iter2 = reader2.iterator();
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
            final int ordering = compareCoordinates(record1, record2);
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

    private void checkPassesFilter(final boolean expected, final SAMRecord record, final String sequence, final int startPos, final int endPos, final boolean contained) {
        final boolean passes = passesFilter(record, sequence, startPos, endPos, contained);
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

    private boolean passesFilter(final SAMRecord record, final String sequence, final int startPos, final int endPos, final boolean contained) {
        if (record == null) {
            return false;
        }
        if (!safeEquals(record.getReferenceName(), sequence)) {
            return false;
        }
        final int alignmentStart = record.getAlignmentStart();
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

    private int compareCoordinates(final SAMRecord record1, final SAMRecord record2) {
        final int seqIndex1 = record1.getReferenceIndex();
        final int seqIndex2 = record2.getReferenceIndex();
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

    private boolean safeEquals(final Object o1, final Object o2) {
        if (o1 == o2) {
            return true;
        } else if (o1 == null || o2 == null) {
            return false;
        } else {
            return o1.equals(o2);
        }
    }

    private void verbose(final String text) {
        if (mVerbose) {
            System.out.println("# " + text);
        }
    }
}