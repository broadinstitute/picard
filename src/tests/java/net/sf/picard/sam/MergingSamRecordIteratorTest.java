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
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.samtools.*;
import net.sf.samtools.util.SequenceUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

/**
 * Tests for MultiIterator
 *
 * @author Dave Tefft
 */
public class MergingSamRecordIteratorTest {

    @Test
    public void testVanillaCoordinateMultiIterator() throws Exception {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder();
        builder1.addFrag("read_28833_29006_6945", 20, 28833, false); // ok
        builder1.addFrag("read_28701_28881_323b", 22, 28834, false); // ok

        final SAMFileReader samReader = builder1.getSamReader();
        samReader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder();
        builder2.addFrag("read_28833_29006_6945", 20, 30000, false); // ok
        builder2.addFrag("read_28701_28881_323b", 22, 28835, false); // ok
        builder2.addFrag("read_28701_28881_323c", 22, 28835, false); // ok

        final SAMFileReader samReader2 = builder2.getSamReader();
        samReader2.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);


        final List<SAMFileReader> readerList = new ArrayList<SAMFileReader>();
        readerList.add(samReader);
        readerList.add(samReader2);

        final SamFileHeaderMerger fileHeaderMerger = new SamFileHeaderMerger(readerList, SAMFileHeader.SortOrder.coordinate);

        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(fileHeaderMerger, false);


        int i = 0;

        // This is the correct order for start bases.  The first two are on chr20, the next three on chr22
        final int[] startBasesInOrder = {28833, 30000, 28834, 28835, 28835};

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            System.out.println(rec.format());
            Assert.assertEquals(rec.getAlignmentStart(), startBasesInOrder[i]);
            i++;
        }
        samReader.close();
    }

    @Test
    public void testVanillaReadOrderMultiIterator() throws Exception {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        builder1.addFrag("a", 20, 28833, false); // ok
        builder1.addFrag("e", 19, 28834, false); // ok

        final SAMFileReader samReader = builder1.getSamReader();

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        builder2.addFrag("b", 20, 30000, false); // ok
        builder2.addFrag("c", 22, 28835, false); // ok
        builder2.addFrag("d", 20, 28835, false); // ok

        final SAMFileReader samReader2 = builder2.getSamReader();


        final List<SAMFileReader> readerList = new ArrayList<SAMFileReader>();
        readerList.add(samReader);
        readerList.add(samReader2);

        final SamFileHeaderMerger fileHeaderMerger = new SamFileHeaderMerger(readerList, SAMFileHeader.SortOrder.queryname);

        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(fileHeaderMerger, false);


        int i = 0;

        // This is the correct order for start bases.  The first two are on chr20, the next three on chr22
        final String[] orderedReadNames = {"a", "b", "c", "d", "e"};

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            System.out.println(rec.getReadName());
            Assert.assertEquals(rec.getReadName(), orderedReadNames[i]);
            i++;
        }
        samReader.close();
    }

    @Test
    public void testVanillaUnsortedMultiIterator() throws Exception {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.unsorted);
        builder1.addFrag("b", 20, 28833, false); // ok
        builder1.addFrag("a", 19, 28834, false); // ok

        final SAMFileReader samReader = builder1.getSamReader();

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.unsorted);
        builder2.addFrag("d", 20, 30000, false); // ok
        builder2.addFrag("e", 22, 28835, false); // ok
        builder2.addFrag("c", 20, 28835, false); // ok

        final SAMFileReader samReader2 = builder2.getSamReader();


        final List<SAMFileReader> readerList = new ArrayList<SAMFileReader>();
        readerList.add(samReader);
        readerList.add(samReader2);

        final SamFileHeaderMerger fileHeaderMerger = new SamFileHeaderMerger(readerList, SAMFileHeader.SortOrder.unsorted);

        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(fileHeaderMerger, false);


        int i = 0;

        // With unsorted option there is no garantee that order of the names to come back from the iterator
        final String[] readNames = {"b", "a", "d", "e", "c"};

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            System.out.println(rec.getReadName());
            i++;
        }
        Assert.assertEquals(i, readNames.length);
        samReader.close();
    }

    @Test(expectedExceptions = SequenceUtil.SequenceListsDifferException.class)
    public void testConflictingHeaders() throws Exception {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder();
        builder1.addFrag("read_28833_29006_6945", 20, 28833, false); // ok
        builder1.addFrag("read_28701_28881_323b", 22, 28834, false); // ok

        final SAMFileReader samReader = builder1.getSamReader();
        samReader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder();
        builder2.addFrag("read_28833_29006_6945", 20, 30000, false); // ok
        builder2.addFrag("read_28701_28881_323b", 22, 28835, false); // ok
        builder2.addFrag("read_28701_28881_323c", 22, 28835, false); // ok

        final SAMFileReader samReader2 = builder2.getSamReader();
        samReader2.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);

        //Change one of the header so they are no longer compatible
        final SAMSequenceRecord sRec = new SAMSequenceRecord("BADSEQ", 0);
        samReader2.getFileHeader().addSequence(sRec);


        final List<SAMFileReader> readerList = new ArrayList<SAMFileReader>();
        readerList.add(samReader);
        readerList.add(samReader2);

        final SamFileHeaderMerger samFileHeaderMerger = new SamFileHeaderMerger(readerList, SAMFileHeader.SortOrder.coordinate);

        new MergingSamRecordIterator(samFileHeaderMerger, false);
        Assert.fail("This method should throw exception before getting to this point");


    }

    @Test
    public void distinctReadGroups() throws Exception {
        final String[] groupIds = {"GROUP1", "GROUP2"};

        final List<SAMReadGroupRecord> samGroups = createGroupRecords(groupIds);
        int i = 0;
        for (final SAMReadGroupRecord g : samGroups) {
            System.out.println("g.getReadGroupId() = " + g.getReadGroupId());
            Assert.assertEquals(groupIds[i], g.getReadGroupId());
            i++;
        }
    }

    @Test
    public void mergeReadGroups() throws Exception {
        final String[] groupIds = {"GROUP1", "GROUP1"};
        final String[] expectedGroupIds = {"A", "B"};

        final List<SAMReadGroupRecord> samGroups = createGroupRecords(groupIds);
        int i = 0;
        for (final SAMReadGroupRecord g : samGroups) {
            System.out.println("g.getReadGroupId() = " + g.getReadGroupId());
            Assert.assertEquals(expectedGroupIds[i], g.getReadGroupId());
            i++;
        }
    }

    @Test(expectedExceptions = PicardException.class)
    public void filesNotSortedCorrectly() throws Exception {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.coordinate);
        builder1.addFrag("read_28833_29006_6945", 20, 28833, false); // ok
        builder1.addFrag("read_28701_28881_323b", 22, 28834, false); // ok

        final SAMFileReader samReader = builder1.getSamReader();

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.unsorted);

        final SAMFileReader samReader2 = builder2.getSamReader();

        builder2.addFrag("read_28701_28881_323b", 22, 28835, false); // ok
        builder2.addFrag("read_28833_29006_6945", 20, 30000, false); // ok
        builder2.addFrag("read_28701_28881_323c", 22, 28835, false); // ok

        final List<SAMFileReader> readerList = new ArrayList<SAMFileReader>();
        readerList.add(samReader);
        readerList.add(samReader2);

        final SamFileHeaderMerger fileHeaderMerger = new SamFileHeaderMerger(readerList, SAMFileHeader.SortOrder.coordinate);

        new MergingSamRecordIterator(fileHeaderMerger, false);
        Assert.fail("This method should throw exception before getting to this point");

    }

    private List<SAMReadGroupRecord> createGroupRecords(final String[] groupIds) {
        final List<SAMFileReader> readerList = new ArrayList<SAMFileReader>();

        for (final String groupId : groupIds) {
            final SAMFileReader samReader = getAFileReader();
            final List<SAMReadGroupRecord> records1 = new ArrayList<SAMReadGroupRecord>();
            final SAMReadGroupRecord r1 = new SAMReadGroupRecord(groupId);
            records1.add(r1);
            samReader.getFileHeader().setReadGroups(records1);
            readerList.add(samReader);
        }

        final SamFileHeaderMerger fileHeaderMerger = new SamFileHeaderMerger(readerList, SAMFileHeader.SortOrder.coordinate);
        final SAMFileHeader header = fileHeaderMerger.getMergedHeader();
        return header.getReadGroups();
    }

    private SAMFileReader getAFileReader() {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder();
        builder1.addFrag("read_28833_29006_6945", 20, 28833, false); // ok
        return builder1.getSamReader();
    }

    /**
     * List of program groups from the input files are merged, and renumbered, and SAMRecords
     * with PG tags get assigned the updated PG ID
     */
    @Test
    public void testMergingProgramGroups() {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        final SAMProgramRecord program1 = new SAMProgramRecord("0");
        program1.setCommandLine("Hi, Mom!");
        builder1.setProgramRecord(program1);
        builder1.addFrag("read1", 20, 28833, false);

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        final SAMProgramRecord program2 = new SAMProgramRecord("0");
        program2.setCommandLine("Hi, Dad!");
        program2.setProgramVersion("123");
        builder2.setProgramRecord(program2);
        builder2.addFrag("read2", 19, 28833, false);
        // No PG tag on this record
        builder2.setProgramRecord(null);
        builder2.addFrag("read3", 19, 28833, false);

        final List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        readers.add(builder1.getSamReader());
        readers.add(builder2.getSamReader());

        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SAMFileHeader.SortOrder.queryname);
        final List<SAMProgramRecord> outputProgramGroups = headerMerger.getMergedHeader().getProgramRecords();
        Assert.assertEquals(outputProgramGroups.size(), 2);
        Assert.assertTrue(outputProgramGroups.get(0).equivalent(program1));
        Assert.assertTrue(outputProgramGroups.get(1).equivalent(program2));

        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, false);
        SAMRecord samRecord = iterator.next();
        Assert.assertEquals(samRecord.getAttribute(SAMTag.PG.name()), "0");
        samRecord = iterator.next();
        Assert.assertEquals(samRecord.getAttribute(SAMTag.PG.name()), "1");
        samRecord = iterator.next();
        Assert.assertEquals(samRecord.getAttribute(SAMTag.PG.name()), null);
        Assert.assertFalse(iterator.hasNext());
    }

    /**
     * Program groups in input files that are equivalent (identical except for ID) are folded into
     * one program group in the output, and renumbered.
     */
    @Test
    public void testMergingEquivalentProgramGroups() {
        final SAMRecordSetBuilder builder1 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        final SAMProgramRecord program1 = new SAMProgramRecord("0");
        program1.setCommandLine("Hi, Mom!");
        builder1.setProgramRecord(program1);
        builder1.addFrag("read1", 20, 28833, false);

        final SAMRecordSetBuilder builder2 = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.queryname);
        final SAMProgramRecord program2 = new SAMProgramRecord("55");
        program2.setCommandLine("Hi, Mom!");
        builder2.setProgramRecord(program2);
        builder2.addFrag("read2", 19, 28833, false);

        final List<SAMFileReader> readers = new ArrayList<SAMFileReader>();
        readers.add(builder1.getSamReader());
        readers.add(builder2.getSamReader());

        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(readers, SAMFileHeader.SortOrder.queryname);
        final List<SAMProgramRecord> outputProgramGroups = headerMerger.getMergedHeader().getProgramRecords();
        Assert.assertEquals(outputProgramGroups.size(), 1);
        Assert.assertTrue(outputProgramGroups.get(0).equivalent(program1));
        Assert.assertTrue(outputProgramGroups.get(0).equivalent(program2));

        final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, false);
        SAMRecord samRecord = iterator.next();
        Assert.assertEquals(samRecord.getAttribute(SAMTag.PG.name()), "0");
        samRecord = iterator.next();
        Assert.assertEquals(samRecord.getAttribute(SAMTag.PG.name()), "0");
        Assert.assertFalse(iterator.hasNext());
    }

}
