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

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.picard.util.CigarUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import sun.tools.asm.TryData;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.*;
import java.io.*;

/**
 *  Test for the MergeBamAlignment class
 *
 * @author ktibbett@broadinstitute.org
 */
public class MergeBamAlignmentTest {

    private static final File unmappedBam = new File("testdata/net/sf/picard/sam/unmapped.sam");
    private static final File alignedBam = new File("testdata/net/sf/picard/sam/aligned.sam");
    private static final File oneHalfAlignedBam = new File("testdata/net/sf/picard/sam/onehalfaligned.sam");
    private static final File otherHalfAlignedBam = new File("testdata/net/sf/picard/sam/otherhalfaligned.sam");
    private static final File mergingUnmappedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/unmapped.sam");
    private static final File firstReadAlignedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/allread1.trimmed.aligned.sam");
    private static final File secondReadAlignedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/allread2.trimmed.aligned.sam");
    private static final File firstReadAlignedBam_firstHalf = new File("testdata/net/sf/picard/sam/MergeBamAlignment/firsthalf.read1.trimmed.aligned.sam");
    private static final File firstReadAlignedBam_secondHalf = new File("testdata/net/sf/picard/sam/MergeBamAlignment/secondhalf.read1.trimmed.aligned.sam");
    private static final File secondReadAlignedBam_firstHalf = new File("testdata/net/sf/picard/sam/MergeBamAlignment/firsthalf.read2.trimmed.aligned.sam");
    private static final File secondReadAlignedBam_secondHalf = new File("testdata/net/sf/picard/sam/MergeBamAlignment/secondhalf.read2.trimmed.aligned.sam");
    private static final File alignedQuerynameSortedBam =
            new File("testdata/net/sf/picard/sam/aligned_queryname_sorted.sam");
    private static final File fasta = new File("testdata/net/sf/picard/sam/merger.fasta");
    private static final File badorderUnmappedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/unmapped.badorder.sam");
    private static final File badorderAlignedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/aligned.badorder.sam");

    @Test
    public void testMerger() throws Exception {
        File output = File.createTempFile("mergeTest", ".sam");
        output.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(new File[]{alignedBam});
        merger.ALIGNED_READS_ONLY = false;
        merger.CLIP_ADAPTERS = true;
        merger.IS_BISULFITE_SEQUENCE = false;
        merger.MAX_INSERTIONS_OR_DELETIONS = 1;
        merger.PROGRAM_RECORD_ID = "0";
        merger.PROGRAM_GROUP_VERSION = "1.0";
        merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
        merger.PROGRAM_GROUP_NAME = "myAligner";
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = output;

        Assert.assertEquals(merger.doWork(), 0, "Merge did not succeed");
        SAMFileReader result = new SAMFileReader(output);
        Assert.assertEquals(result.getFileHeader().getSequenceDictionary().getSequences().size(), 8,
                "Number of sequences did not match");
        SAMProgramRecord pg = result.getFileHeader().getProgramRecords().get(0);
        Assert.assertEquals(pg.getProgramGroupId(), merger.PROGRAM_RECORD_ID);
        Assert.assertEquals(pg.getProgramVersion(), merger.PROGRAM_GROUP_VERSION);
        Assert.assertEquals(pg.getCommandLine(), merger.PROGRAM_GROUP_COMMAND_LINE);
        Assert.assertEquals(pg.getProgramName(), merger.PROGRAM_GROUP_NAME);

        SAMReadGroupRecord rg = result.getFileHeader().getReadGroups().get(0);
        Assert.assertEquals(rg.getReadGroupId(), "0");
        Assert.assertEquals(rg.getSample(), "Hi,Mom!");

        Assert.assertEquals(result.getFileHeader().getSortOrder(), SAMFileHeader.SortOrder.coordinate);

        for (CloseableIterator<SAMRecord> it = result.iterator(); it.hasNext();) {
            SAMRecord sam = it.next();
            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_adapter") ||
                sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " +
                        sam.getReadName());
                }
                else {
                    Assert.assertEquals(sam.getCigarString(), "96M5S", "Incorrect CIGAR string for " +
                        sam.getReadName());
                }
            }
            // This tests that we DON'T clip when we run off the end if there are equal to or more than
            // MIN_ADAPTER_BASES hanging off the end
            else if (sam.getReadName().equals("both_reads_align_min_adapter_bases_exceeded")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                Assert.assertTrue(sam.getCigarString().indexOf("S") == -1,
                        "Read was clipped when it should not be.");
            }
            else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                     sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                }
                else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            }
            else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }

        }

        // Quick test to make sure the program record gets picked up from the file if not specified
        // on the command line.
        merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(new File[]{alignedBam});
        merger.ALIGNED_READS_ONLY = false;
        merger.CLIP_ADAPTERS = true;
        merger.IS_BISULFITE_SEQUENCE = false;
        merger.MAX_INSERTIONS_OR_DELETIONS = 1;
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = output;

        Assert.assertEquals(merger.doWork(), 0);

        result = new SAMFileReader(output);
        pg = result.getFileHeader().getProgramRecords().get(0);
        Assert.assertEquals(pg.getProgramGroupId(), "1",
                "Program group ID not picked up correctly from aligned BAM");
        Assert.assertEquals(pg.getProgramVersion(), "2.0",
                "Program version not picked up correctly from aligned BAM");
        Assert.assertNull(pg.getCommandLine(),
                "Program command line not picked up correctly from aligned BAM");
        Assert.assertEquals(pg.getProgramName(), "Hey!",
                "Program name not picked up correctly from aligned BAM");

    }


    @Test
    public void testMergerFromMultipleFiles() throws Exception {
        File output = File.createTempFile("mergeTest", ".sam");
        output.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(new File[]{oneHalfAlignedBam, otherHalfAlignedBam});
        merger.ALIGNED_READS_ONLY = false;
        merger.CLIP_ADAPTERS = true;
        merger.IS_BISULFITE_SEQUENCE = false;
        merger.MAX_INSERTIONS_OR_DELETIONS = 1;
        merger.PROGRAM_RECORD_ID = "0";
        merger.PROGRAM_GROUP_VERSION = "1.0";
        merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
        merger.PROGRAM_GROUP_NAME = "myAligner";
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = output;

        Assert.assertEquals(merger.doWork(), 0, "Merge did not succeed");
        SAMFileReader result = new SAMFileReader(output);

        for (CloseableIterator<SAMRecord> it = result.iterator(); it.hasNext();) {
            SAMRecord sam = it.next();
            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_adapter") ||
                sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " +
                        sam.getReadName());
                }
                else {
                    Assert.assertEquals(sam.getCigarString(), "96M5S", "Incorrect CIGAR string for " +
                        sam.getReadName());
                }
            }
            // This tests that we DON'T clip when we run off the end if there are equal to or more than
            // MIN_ADAPTER_BASES hanging off the end
            else if (sam.getReadName().equals("both_reads_align_min_adapter_bases_exceeded")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                Assert.assertTrue(sam.getCigarString().indexOf("S") == -1,
                        "Read was clipped when it should not be.");
            }
            else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                     sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                }
                else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            }
            else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }

        }

    }

    @Test(dataProvider="data")
    public void testSortingOnSamAlignmentMerger(File unmapped, File aligned, boolean sorted, String testName)
        throws IOException {

        File target = File.createTempFile("target", "bam");
        SamAlignmentMerger merger = new SamAlignmentMerger(unmapped,  target, fasta, null, true, false,
                 true, false, false, aligned, 1);
        merger.mergeAlignment();
        Assert.assertEquals(sorted, !merger.getForceSort());
        SAMRecordIterator it = new SAMFileReader(target).iterator();
        int aln = 0;
        int ct = 0;
        while (it.hasNext()) {
            SAMRecord rec = it.next();
            if (!rec.getReadUnmappedFlag()) {
                aln++;
            }
        }
        Assert.assertEquals(aln, 6, "Incorrect number of aligned reads in merged BAM file");
        target.delete();
    }

    @DataProvider(name="data")
    public Object[][] getDataForSortingTest() {
        return new Object[][] {
                {unmappedBam, alignedQuerynameSortedBam, true, "Basic test with pre-sorted alignment"},
                {unmappedBam, alignedBam, false, "Basic test with unsorted alignment"}

        };
    }


    /**
     * Minimal test of merging data from separate read 1 and read 2 alignments
     */
    @Test(dataProvider="separateTrimmed")
    public void testMergingFromSeparatedReadTrimmedAlignments(File unmapped, List<File> r1Align, List<File> r2Align, int r1Trim, int r2Trim, String testName) throws Exception {
         File output = File.createTempFile("mergeMultileAlignmentsTest", ".sam");
         output.deleteOnExit();

System.out.println(testName);
         MergeBamAlignment merger = new MergeBamAlignment();
         merger.UNMAPPED_BAM = unmapped;
         merger.ALIGNED_READS_ONLY = false;
         merger.CLIP_ADAPTERS = true;
         merger.IS_BISULFITE_SEQUENCE = false;
         merger.MAX_INSERTIONS_OR_DELETIONS = 1;
         merger.PROGRAM_RECORD_ID = "0";
         merger.PROGRAM_GROUP_VERSION = "1.0";
         merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
         merger.PROGRAM_GROUP_NAME = "myAligner";
         merger.PAIRED_RUN = true;
         merger.REFERENCE_SEQUENCE = fasta;
         merger.READ1_ALIGNED_BAM = r1Align;
         merger.READ2_ALIGNED_BAM = r2Align;
         merger.READ1_TRIM = r1Trim;
         merger.READ2_TRIM = r2Trim;
         merger.OUTPUT = output;

         Assert.assertEquals(merger.doWork(), 0, "Merge did not succeed: " + testName);
         SAMFileReader result = new SAMFileReader(output);
         SAMProgramRecord pg = result.getFileHeader().getProgramRecords().get(0);


         for (CloseableIterator<SAMRecord> it = result.iterator(); it.hasNext();) {
             SAMRecord sam = it.next();

             // Get the alignment record
             List<File> rFiles = sam.getFirstOfPairFlag() ? r1Align : r2Align;
             SAMRecord alignment = null;
             for (File f : rFiles) {
                 Iterator<SAMRecord> it2 = new SAMFileReader(f).iterator();
                 while (it2.hasNext()) {
                     SAMRecord tmp = it2.next();
                     if (tmp.getReadName().equals(sam.getReadName())) {
                         alignment = tmp;
                         break;
                     }
                 }
                 if (alignment != null) break;
             }

             int trim = sam.getFirstOfPairFlag() ? r1Trim : r2Trim;
             int notWrittenToFastq = sam.getReadLength() - (trim + alignment.getReadLength());
             int beginning = sam.getReadNegativeStrandFlag() ? notWrittenToFastq : trim;
             int end = sam.getReadNegativeStrandFlag() ? trim : notWrittenToFastq;

             if (!sam.getReadUnmappedFlag()) {
                 CigarElement firstMergedCigarElement = sam.getCigar().getCigarElement(0);
                 CigarElement lastMergedCigarElement = sam.getCigar().getCigarElement(sam.getCigar().getCigarElements().size()-1);
                 CigarElement firstAlignedCigarElement = alignment.getCigar().getCigarElement(0);
                 CigarElement lastAlignedCigarElement = alignment.getCigar().getCigarElement(alignment.getCigar().getCigarElements().size()-1);

                 if (beginning > 0) {
                     Assert.assertEquals(firstMergedCigarElement.getOperator(), CigarOperator.S, "First element is not a soft clip");
                     Assert.assertEquals(firstMergedCigarElement.getLength(), beginning + ((firstAlignedCigarElement.getOperator() == CigarOperator.S) ? firstAlignedCigarElement.getLength() : 0));
                 }
                 if (end > 0) {
                     Assert.assertEquals(lastMergedCigarElement.getOperator(), CigarOperator.S, "Last element is not a soft clip");
                     Assert.assertEquals(lastMergedCigarElement.getLength(), end + ((lastAlignedCigarElement.getOperator() == CigarOperator.S) ? lastAlignedCigarElement.getLength() : 0));
                 }
             }
         }

    }

    @DataProvider(name="separateTrimmed")
    public Object[][] getDataForMergingTest() {
        return new Object[][] {
            // Tests using one file for each read, queryname sorted, trimmed, not all bases written
            {mergingUnmappedBam, Arrays.asList(new File[]{firstReadAlignedBam}), Arrays.asList(new File[]{secondReadAlignedBam}), 17, 20, "one file per read"},
            // Tests using multiple files for each read, coordinate sorted, trimmed, not all bases written
            {mergingUnmappedBam, Arrays.asList(new File[]{firstReadAlignedBam_firstHalf, firstReadAlignedBam_secondHalf}),
                    Arrays.asList(new File[]{secondReadAlignedBam_firstHalf, secondReadAlignedBam_secondHalf}), 17, 20, "two files per read"}
        };
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testOldQuerynameSortFails() throws IOException {

        File merged = File.createTempFile("merged", ".bam");
        merged.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = badorderUnmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(new File[]{badorderAlignedBam});
        merger.ALIGNED_READS_ONLY = false;
        merger.CLIP_ADAPTERS = true;
        merger.IS_BISULFITE_SEQUENCE = false;
        merger.MAX_INSERTIONS_OR_DELETIONS = 1;
        merger.PROGRAM_RECORD_ID = "0";
        merger.PROGRAM_GROUP_VERSION = "1.0";
        merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
        merger.PROGRAM_GROUP_NAME = "myAligner";
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = merged;
        merger.PAIRED_RUN = true;
        merger.doWork();
        Assert.fail("Merger should have failed because unmapped reads are not in queryname order but didn't");
    }

}
