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
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 *  Test for the MergeBamAlignment class
 *
 * @author ktibbett@broadinstitute.org
 */
public class MergeBamAlignmentTest {

    private static final File unmappedBam = new File("testdata/net/sf/picard/sam/unmapped.sam");
    private static final File alignedBam = new File("testdata/net/sf/picard/sam/aligned.sam");
    private static final File alignedQuerynameSortedBam =
            new File("testdata/net/sf/picard/sam/aligned_queryname_sorted.sam");
    private static final File fasta = new File("testdata/net/sf/picard/sam/merger.fasta");

    @Test
    public void testMerger() throws Exception {
        File output = File.createTempFile("mergeTest", ".sam");
        output.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedBam;
        merger.ALIGNED_BAM = alignedBam;
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
        merger.ALIGNED_BAM = alignedBam;
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


}
