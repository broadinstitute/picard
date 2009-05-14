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

import net.sf.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;

import net.sf.samtools.*;

/**
 * Test that BAM writing doesn't blow up.  For presorted writing, the resulting BAM file is read and contents are
 * compared with the original SAM file. 
 */
public class BAMFileWriterTest {

    private SAMRecordSetBuilder getSAMReader(final boolean sortForMe, final SAMFileHeader.SortOrder sortOrder) {
        final SAMRecordSetBuilder ret = new SAMRecordSetBuilder(sortForMe, sortOrder);
        ret.addPair("readB", 20, 200, 300);
        ret.addPair("readA", 20, 100, 150);
        ret.addFrag("readC", 20, 140, true);
        ret.addFrag("readD", 20, 140, false);
        return ret;
    }

    /**
     * Parse some SAM text into a SAM object, then write as BAM.  If SAM text was presorted, then the BAM file can
     * be read and compared with the SAM object.
     * @param samRecordSetBuilder source of input SAMFileReader to be written and compared with
     * @param sortOrder How the BAM should be written
     * @param presorted If true, samText is in the order specified by sortOrder
     */
    private void testHelper(final SAMRecordSetBuilder samRecordSetBuilder, final SAMFileHeader.SortOrder sortOrder, final boolean presorted) throws Exception {
        SAMFileReader samReader = samRecordSetBuilder.getSamReader();
        final File bamFile = File.createTempFile("test.", ".bam");
        bamFile.deleteOnExit();
        samReader.getFileHeader().setSortOrder(sortOrder);
        final SAMFileWriter bamWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samReader.getFileHeader(), presorted, bamFile);
        CloseableIterator<SAMRecord> it = samReader.iterator();
        while (it.hasNext()) {
            bamWriter.addAlignment(it.next());
        }
        bamWriter.close();
        it.close();
        samReader.close();

        if (presorted) {
            // If SAM text input was presorted, then we can compare SAM object to BAM object
            final SAMFileReader bamReader = new SAMFileReader(bamFile);
            samReader = samRecordSetBuilder.getSamReader();
            samReader.getFileHeader().setSortOrder(bamReader.getFileHeader().getSortOrder());
            Assert.assertEquals(bamReader.getFileHeader(), samReader.getFileHeader());
            it = samReader.iterator();
            final CloseableIterator<SAMRecord> bamIt = bamReader.iterator();
            while (it.hasNext()) {
                Assert.assertTrue(bamIt.hasNext());
                final SAMRecord samRecord = it.next();
                final SAMRecord bamRecord = bamIt.next();

                // SAMRecords don't have this set, so stuff it in there
                samRecord.setIndexingBin(bamRecord.getIndexingBin());

                // Force reference index attributes to be populated
                samRecord.getReferenceIndex();
                bamRecord.getReferenceIndex();
                samRecord.getMateReferenceIndex();
                bamRecord.getMateReferenceIndex();

                Assert.assertEquals(bamRecord, samRecord);
            }
            Assert.assertFalse(bamIt.hasNext());
        }
    }

    @DataProvider(name = "test1")
    public Object[][] createTestData() {
        return new Object[][] {
                {"coordinate sorted", getSAMReader(false, SAMFileHeader.SortOrder.unsorted), SAMFileHeader.SortOrder.coordinate, false},
                {"query sorted", getSAMReader(false, SAMFileHeader.SortOrder.unsorted), SAMFileHeader.SortOrder.queryname, false},
                {"unsorted", getSAMReader(false, SAMFileHeader.SortOrder.unsorted), SAMFileHeader.SortOrder.unsorted, false},
                {"coordinate presorted", getSAMReader(true, SAMFileHeader.SortOrder.coordinate), SAMFileHeader.SortOrder.coordinate, true},
                {"query presorted", getSAMReader(true, SAMFileHeader.SortOrder.queryname), SAMFileHeader.SortOrder.queryname, true},
        };
    }

    @Test(dataProvider = "test1")
    public void testPositive(final String testName, final SAMRecordSetBuilder samRecordSetBuilder, final SAMFileHeader.SortOrder order, final boolean presorted) throws Exception {

        testHelper(samRecordSetBuilder, order, presorted);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativePresorted() throws Exception {

        testHelper(getSAMReader(true, SAMFileHeader.SortOrder.coordinate), SAMFileHeader.SortOrder.queryname, true);
        Assert.fail("Exception should be thrown");
    }
}
