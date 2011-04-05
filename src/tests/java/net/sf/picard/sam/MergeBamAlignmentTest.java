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
import net.sf.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
    private static final String bigSequenceName = "chr7"; // The longest sequence in merger.fasta
    private static final File sequenceDict = new File("testdata/net/sf/picard/sam/merger.dict");
    private static final File badorderUnmappedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/unmapped.badorder.sam");
    private static final File badorderAlignedBam = new File("testdata/net/sf/picard/sam/MergeBamAlignment/aligned.badorder.sam");
    private static final File TESTDATA_DIR = new File("testdata/net/sf/picard/sam/MergeBamAlignment");

    @Test
    public void testMerger() throws Exception {
        File output = File.createTempFile("mergeTest", ".sam");
        output.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(alignedBam);
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
        merger.EXPECTED_ORIENTATIONS=Arrays.asList(SamPairUtil.PairOrientation.FR);

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

        for (SAMRecord sam : result) {
            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_adapter") ||
                    sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " +
                            sam.getReadName());
                } else {
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
            } else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                    sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                } else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            } else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }

        }

        // Quick test to make sure the program record gets picked up from the file if not specified
        // on the command line.
        merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(alignedBam);
        merger.ALIGNED_READS_ONLY = false;
        merger.CLIP_ADAPTERS = true;
        merger.IS_BISULFITE_SEQUENCE = false;
        merger.MAX_INSERTIONS_OR_DELETIONS = 1;
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = output;
        merger.EXPECTED_ORIENTATIONS=Arrays.asList(SamPairUtil.PairOrientation.FR);

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
        merger.ALIGNED_BAM = Arrays.asList(oneHalfAlignedBam, otherHalfAlignedBam);
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
        merger.EXPECTED_ORIENTATIONS=Arrays.asList(SamPairUtil.PairOrientation.FR);

        Assert.assertEquals(merger.doWork(), 0, "Merge did not succeed");
        SAMFileReader result = new SAMFileReader(output);

        for (SAMRecord sam : result) {
            // This tests that we clip both (a) when the adapter is marked in the unmapped BAM file and
            // (b) when the insert size is less than the read length
            if (sam.getReadName().equals("both_reads_align_clip_adapter") ||
                    sam.getReadName().equals("both_reads_align_clip_marked")) {
                Assert.assertEquals(sam.getReferenceName(), "chr7");
                if (sam.getReadNegativeStrandFlag()) {
                    Assert.assertEquals(sam.getCigarString(), "5S96M", "Incorrect CIGAR string for " +
                            sam.getReadName());
                } else {
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
            } else if (sam.getReadName().equals("neither_read_aligns_or_present")) {
                Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should be unmapped but isn't");
            }
            // Two pairs in which only the first read should align
            else if (sam.getReadName().equals("both_reads_present_only_first_aligns") ||
                    sam.getReadName().equals("read_2_too_many_gaps")) {
                if (sam.getFirstOfPairFlag()) {
                    Assert.assertEquals(sam.getReferenceName(), "chr7", "Read should be mapped but isn't");
                } else {
                    Assert.assertTrue(sam.getReadUnmappedFlag(), "Read should not be mapped but is");
                }
            } else {
                throw new Exception("Unexpected read name: " + sam.getReadName());
            }

        }

    }

    @Test(dataProvider="data")
    public void testSortingOnSamAlignmentMerger(File unmapped, File aligned, boolean sorted, String testName)
        throws IOException {

        File target = File.createTempFile("target", "bam");
        target.deleteOnExit();
        SamAlignmentMerger merger = new SamAlignmentMerger(unmapped,  target, fasta, null, true, false,
                 true, false, Arrays.asList(aligned), 1, null, null, null, null, null,
                Arrays.asList(SamPairUtil.PairOrientation.FR));

        merger.mergeAlignment();
        Assert.assertEquals(sorted, !merger.getForceSort());
        SAMRecordIterator it = new SAMFileReader(target).iterator();
        int aln = 0;
        while (it.hasNext()) {
            SAMRecord rec = it.next();
            if (!rec.getReadUnmappedFlag()) {
                aln++;
            }
        }
        Assert.assertEquals(aln, 6, "Incorrect number of aligned reads in merged BAM file");
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
         merger.EXPECTED_ORIENTATIONS=Arrays.asList(SamPairUtil.PairOrientation.FR);

         Assert.assertEquals(merger.doWork(), 0, "Merge did not succeed: " + testName);
         SAMFileReader result = new SAMFileReader(output);
         SAMProgramRecord pg = result.getFileHeader().getProgramRecords().get(0);


        for (SAMRecord sam : result) {
            // Get the alignment record
            List<File> rFiles = sam.getFirstOfPairFlag() ? r1Align : r2Align;
            SAMRecord alignment = null;
            for (File f : rFiles) {
                for (SAMRecord tmp : new SAMFileReader(f)) {
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
                CigarElement lastMergedCigarElement = sam.getCigar().getCigarElement(sam.getCigar().getCigarElements().size() - 1);
                CigarElement firstAlignedCigarElement = alignment.getCigar().getCigarElement(0);
                CigarElement lastAlignedCigarElement = alignment.getCigar().getCigarElement(alignment.getCigar().getCigarElements().size() - 1);

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
            {mergingUnmappedBam, Arrays.asList(firstReadAlignedBam), Arrays.asList(secondReadAlignedBam), 17, 20, "one file per read"},
            // Tests using multiple files for each read, coordinate sorted, trimmed, not all bases written
            {mergingUnmappedBam, Arrays.asList(firstReadAlignedBam_firstHalf, firstReadAlignedBam_secondHalf),
                    Arrays.asList(secondReadAlignedBam_firstHalf, secondReadAlignedBam_secondHalf), 17, 20, "two files per read"}
        };
    }

    @Test(expectedExceptions = {IllegalStateException.class, PicardException.class})
    public void testOldQuerynameSortFails() throws IOException {

        File merged = File.createTempFile("merged", ".bam");
        merged.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = badorderUnmappedBam;
        merger.ALIGNED_BAM = Arrays.asList(badorderAlignedBam);
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
        merger.EXPECTED_ORIENTATIONS=Arrays.asList(SamPairUtil.PairOrientation.FR);
        merger.doWork();
        Assert.fail("Merger should have failed because unmapped reads are not in queryname order but didn't");
    }

    @Test
    public void testMultiHit() throws IOException {
        final File unmappedSam = new File(TESTDATA_DIR, "multihit.unmapped.sam");
        final File alignedSam = new File(TESTDATA_DIR, "multihit.aligned.sam");
        File merged = File.createTempFile("merged", ".sam");
        merged.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedSam;
        merger.ALIGNED_BAM = Arrays.asList(alignedSam);
        merger.REFERENCE_SEQUENCE = fasta;
        merger.ALIGNED_READS_ONLY = false;
        merger.PROGRAM_RECORD_ID = "0";
        merger.PROGRAM_GROUP_VERSION = "1.0";
        merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
        merger.PROGRAM_GROUP_NAME = "myAligner";
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = merged;
        merger.PAIRED_RUN = true;
        merger.doWork();

        // Iterate over the merged output and gather some statistics
        final Map<String, AlignmentAccumulator> accumulatorMap = new HashMap<String, AlignmentAccumulator>();

        final SAMFileReader reader = new SAMFileReader(merged);
        for (final SAMRecord rec : reader) {
            final String readName;
            if (!rec.getReadPairedFlag()) readName = rec.getReadName();
            else if (rec.getFirstOfPairFlag()) readName = rec.getReadName() + "/1";
            else readName = rec.getReadName() + "/2";
            AlignmentAccumulator acc = accumulatorMap.get(readName);
            if (acc == null) {
                acc = new AlignmentAccumulator();
                accumulatorMap.put(readName, acc);
            }
            if (rec.getReadUnmappedFlag()) {
                Assert.assertFalse(acc.seenUnaligned, readName);
                acc.seenUnaligned = true;
            } else {
                ++acc.numAlignments;
                if (!rec.getNotPrimaryAlignmentFlag()) {
                    Assert.assertNull(acc.primaryAlignmentSequence, readName);
                    acc.primaryAlignmentSequence = rec.getReferenceName();
                }
            }
        }

        // Set up expected results
        final Map<String, AlignmentAccumulator> expectedMap = new HashMap<String, AlignmentAccumulator>();
        expectedMap.put("frag_hit", new AlignmentAccumulator(false, 1, "chr1"));
        expectedMap.put("frag_multihit", new AlignmentAccumulator(false, 3, "chr2"));
        expectedMap.put("frag_no_hit", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_both_hit/1", new AlignmentAccumulator(false, 1, "chr7"));
        expectedMap.put("pair_both_hit/2", new AlignmentAccumulator(false, 1, "chr7"));
        expectedMap.put("pair_both_multihit/1", new AlignmentAccumulator(false, 3, "chr8"));
        expectedMap.put("pair_both_multihit/2", new AlignmentAccumulator(false, 3, "chr8"));
        expectedMap.put("pair_first_hit/1", new AlignmentAccumulator(false, 1, "chr1"));
        expectedMap.put("pair_first_hit/2", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_first_multihit/1", new AlignmentAccumulator(false, 3, "chr1"));
        expectedMap.put("pair_first_multihit/2", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_no_hit/1", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_no_hit/2", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_second_hit/1", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_second_hit/2", new AlignmentAccumulator(false, 1, "chr1"));
        expectedMap.put("pair_second_multihit/1", new AlignmentAccumulator(true, 0, null));
        expectedMap.put("pair_second_multihit/2", new AlignmentAccumulator(false, 3, "chr3"));

        Assert.assertEquals(accumulatorMap.size(), expectedMap.size());

        for (Map.Entry<String, AlignmentAccumulator> entry : expectedMap.entrySet()) {
            final AlignmentAccumulator actual = accumulatorMap.get(entry.getKey());
            Assert.assertEquals(actual, entry.getValue(), entry.getKey());
        }

        assertSamValid(merged);
    }

    private void assertSamValid(final File sam) {
        final ValidateSamFile validator = new ValidateSamFile();
        validator.INPUT = sam;
        Assert.assertEquals(validator.doWork(), 0);
    }

    private static class AlignmentAccumulator {
        boolean seenUnaligned = false;
        int numAlignments = 0;
        String primaryAlignmentSequence = null;

        private AlignmentAccumulator() {}

        private AlignmentAccumulator(boolean seenUnaligned, int numAlignments, String primaryAlignmentSequence) {
            this.seenUnaligned = seenUnaligned;
            this.numAlignments = numAlignments;
            this.primaryAlignmentSequence = primaryAlignmentSequence;
        }

        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            AlignmentAccumulator that = (AlignmentAccumulator) o;

            if (numAlignments != that.numAlignments) return false;
            if (seenUnaligned != that.seenUnaligned) return false;
            if (primaryAlignmentSequence != null ? !primaryAlignmentSequence.equals(that.primaryAlignmentSequence) : that.primaryAlignmentSequence != null)
                return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = (seenUnaligned ? 1 : 0);
            result = 31 * result + numAlignments;
            result = 31 * result + (primaryAlignmentSequence != null ? primaryAlignmentSequence.hashCode() : 0);
            return result;
        }

        @Override
        public String toString() {
            return "AlignmentAccumulator{" +
                    "seenUnaligned=" + seenUnaligned +
                    ", numAlignments=" + numAlignments +
                    ", primaryAlignmentSequence='" + primaryAlignmentSequence + '\'' +
                    '}';
        }
    }

    /**
     * Read a single paired-end read from a file, and create one or more aligned records for the read pair based on
     * the lists, merge with the original paired-end read, and assert expected results.
     * @param description
     * @param firstOfPair List that describes the aligned SAMRecords to create for the first end.
     * @param secondOfPair List that describes the aligned SAMRecords to create for the second end.
     * @param expectedPrimaryHitIndex Expected value for the HI tag in the primary alignment in the merged output.
     * @param expectedNumFirst Expected number of first-of-pair SAMRecords in the merged output.
     * @param expectedNumSecond Expected number of second-of-pair SAMRecords in the merged output.
     * @param expectedPrimaryMapq Sum of mapqs of both ends of primary alignment in the merged output.
     * @throws Exception
     */
    @Test(dataProvider = "testPairedMultiHitWithFilteringTestCases")
    public void testPairedMultiHitWithFiltering(String description, List<HitSpec> firstOfPair, List<HitSpec> secondOfPair,
                                                Integer expectedPrimaryHitIndex, int expectedNumFirst,
                                                int expectedNumSecond, int expectedPrimaryMapq) throws Exception {

        // Create the aligned file by copying bases, quals, readname from the unmapped read, and conforming to each HitSpec.
        final File unmappedSam = new File(TESTDATA_DIR, "multihit.filter.unmapped.sam");
        SAMRecordIterator unmappedSamFileIterator = new SAMFileReader(unmappedSam).iterator();
        final SAMRecord firstUnmappedRec = unmappedSamFileIterator.next();
        final SAMRecord secondUnmappedRec = unmappedSamFileIterator.next();
        unmappedSamFileIterator.close();
        final File alignedSam = File.createTempFile("aligned.", ".sam");
        alignedSam.deleteOnExit();
        final SAMFileHeader alignedHeader = new SAMFileHeader();
        alignedHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        alignedHeader.setSequenceDictionary(new SAMFileReader(sequenceDict).getFileHeader().getSequenceDictionary());
        final SAMFileWriter alignedWriter = new SAMFileWriterFactory().makeSAMWriter(alignedHeader, true, alignedSam);
        for (int i = 0; i < Math.max(firstOfPair.size(), secondOfPair.size()); ++i) {
            final HitSpec firstHitSpec = firstOfPair.isEmpty()? null: firstOfPair.get(i);
            final HitSpec secondHitSpec = secondOfPair.isEmpty()? null: secondOfPair.get(i);
            final SAMRecord first = makeRead(alignedHeader, firstUnmappedRec, firstHitSpec, true, i);
            final SAMRecord second = makeRead(alignedHeader, secondUnmappedRec, secondHitSpec, false, i);
            if (first != null && second != null) SamPairUtil.setMateInfo(first, second, alignedHeader);
            if (first != null) {
                if (second == null) first.setMateUnmappedFlag(true);
                alignedWriter.addAlignment(first);
            }
            if (second != null) {
                if (first == null) second.setMateUnmappedFlag(true);
                alignedWriter.addAlignment(second);
            }
        }
        alignedWriter.close();

        // Merge aligned file with original unmapped file.
        final File mergedSam = File.createTempFile("merged.", ".sam");
        mergedSam.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedSam;
        merger.ALIGNED_BAM = Arrays.asList(alignedSam);
        merger.REFERENCE_SEQUENCE = fasta;
        merger.ALIGNED_READS_ONLY = false;
        merger.PROGRAM_RECORD_ID = "0";
        merger.PROGRAM_GROUP_VERSION = "1.0";
        merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
        merger.PROGRAM_GROUP_NAME = "myAligner";
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = mergedSam;
        merger.PAIRED_RUN = true;
        merger.doWork();

        assertSamValid(mergedSam);

        // Tally metrics and check for agreement with expected.
        final SAMFileReader mergedReader = new SAMFileReader(mergedSam);
        int numFirst = 0;
        int numSecond = 0;
        Integer primaryHitIndex = null;
        int primaryMapq = 0;
        for (final SAMRecord rec : mergedReader) {
            if (rec.getFirstOfPairFlag()) ++numFirst;
            if (rec.getSecondOfPairFlag()) ++numSecond;
            if (!rec.getNotPrimaryAlignmentFlag()  && !rec.getReadUnmappedFlag()) {
                final Integer hitIndex = rec.getIntegerAttribute(SAMTag.HI.name());
                int newHitIndex = (hitIndex == null? -1: hitIndex);
                if (primaryHitIndex == null) primaryHitIndex = newHitIndex;
                else Assert.assertEquals(newHitIndex, primaryHitIndex.intValue());
                primaryMapq += rec.getMappingQuality();
            }
        }
        Assert.assertEquals(numFirst, expectedNumFirst);
        Assert.assertEquals(numSecond, expectedNumSecond);
        Assert.assertEquals(primaryHitIndex, expectedPrimaryHitIndex);
        Assert.assertEquals(primaryMapq, expectedPrimaryMapq);
    }

    private SAMRecord makeRead(SAMFileHeader alignedHeader, SAMRecord unmappedRec, HitSpec hitSpec,
                               boolean firstOfPair, int hitIndex) {
        if (hitSpec == null) return null;
        final SAMRecord rec = makeRead(alignedHeader, unmappedRec, hitSpec, hitIndex);
        rec.setReadPairedFlag(true);
        if (firstOfPair) {
            rec.setFirstOfPairFlag(true);
            rec.setAlignmentStart(hitIndex + 1);
        } else {
            rec.setSecondOfPairFlag(true);
            rec.setAlignmentStart(hitIndex + 201);
        }
        return rec;
    }

    private SAMRecord makeRead(SAMFileHeader alignedHeader, SAMRecord unmappedRec, HitSpec hitSpec, int hitIndex) {
        final SAMRecord rec = new SAMRecord(alignedHeader);
        rec.setReadName(unmappedRec.getReadName());
        rec.setReadBases(unmappedRec.getReadBases());
        rec.setBaseQualities(unmappedRec.getBaseQualities());
        rec.setMappingQuality(hitSpec.mapq);
        if (!hitSpec.primary) rec.setNotPrimaryAlignmentFlag(true);
        final Cigar cigar = new Cigar();
        int readLength = rec.getReadLength();
        if (hitSpec.filtered) {
            // Add two insertions so alignment is filtered.
            cigar.add(new CigarElement(readLength/2, CigarOperator.M));
            cigar.add(new CigarElement(1, CigarOperator.I));
            cigar.add(new CigarElement(1, CigarOperator.M));
            cigar.add(new CigarElement(1, CigarOperator.I));
            cigar.add(new CigarElement(readLength - (readLength/2 + 3), CigarOperator.M));
        } else {
            cigar.add(new CigarElement(readLength, CigarOperator.M));
        }
        rec.setCigar(cigar);

        rec.setReferenceName(bigSequenceName);
        rec.setAttribute(SAMTag.HI.name(), hitIndex);
        rec.setAlignmentStart(hitIndex + 1);

        return rec;
    }

    class HitSpec {
        final boolean primary;
        final boolean filtered;
        final int mapq;

        HitSpec(boolean primary, boolean filtered, int mapq) {
            this.primary = primary;
            this.filtered = filtered;
            this.mapq = mapq;
        }
    }

    @DataProvider(name="testPairedMultiHitWithFilteringTestCases")
    public Object[][] testPairedMultiHitWithFilteringTestCases() {
        ArrayList<Object[]> ret = new ArrayList<Object[]>();
        List<HitSpec> firstOfPair;
        List<HitSpec> secondOfPair;

        // One end aligned test cases...

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = Collections.emptyList();
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Only 1st end has alignments, primary alignment is filtered, one 2ndary is filtered", firstOfPair, secondOfPair, -1, 1, 1, 9});

        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair = Collections.emptyList();
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, true, 10));
        secondOfPair.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"Only 2nd end has alignments, primary alignment is filtered, one 2ndary is filtered", firstOfPair, secondOfPair, -1, 1, 1, 11});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = Collections.emptyList();
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"Only 1st end has alignments, all are filtered", firstOfPair, secondOfPair, null, 1, 1, 0});

        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair = Collections.emptyList();
        secondOfPair.add(new HitSpec(false, true, 10));
        secondOfPair.add(new HitSpec(true, true, 10));
        secondOfPair.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"Only 2nd end has alignments, all are filtered", firstOfPair, secondOfPair, null, 1, 1, 0});

        // Both ends aligned test cases...

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 9));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        secondOfPair.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Both ends aligned, one end of primary filtered", firstOfPair, secondOfPair, 1, 3, 3, 30});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, true, 9));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        secondOfPair.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Both ends aligned, one end of primary filtered, one end of secondary filtered", firstOfPair, secondOfPair, 1, 2, 3, 30});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 10));
        secondOfPair.add(new HitSpec(false, false, 9));
        secondOfPair.add(new HitSpec(true, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Both ends aligned, both ends of primary filtered", firstOfPair, secondOfPair, 1, 2, 2, 21});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 10));
        secondOfPair.add(new HitSpec(false, false, 12));
        secondOfPair.add(new HitSpec(true, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Both ends aligned, both ends of primary filtered, two secondary alignments with identical mapq",
                firstOfPair, secondOfPair, 1, 2, 2, 21});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, true, 9));
        firstOfPair.add(new HitSpec(true, true, 10));
        firstOfPair.add(new HitSpec(false, false, 10));
        secondOfPair.add(new HitSpec(false, false, 30));
        secondOfPair.add(new HitSpec(true, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Both ends aligned, both ends of primary filtered, one end of secondary filtered, but with high mapq",
                firstOfPair, secondOfPair, 0, 2, 2, 30});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, true, 9));
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, false, 10));
        secondOfPair.add(new HitSpec(false, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        ret.add(new Object[]{"Both ends aligned, secondary filtered on both ends", firstOfPair, secondOfPair, 1, 2, 2, 40});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(false, false, 10));
        firstOfPair.add(new HitSpec(true, false, 10));
        secondOfPair.add(new HitSpec(false, true, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        ret.add(new Object[]{"Both ends aligned, secondary filtered on second end", firstOfPair, secondOfPair, 2, 3, 2, 40});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, false, 9));
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, false, 10));
        secondOfPair.add(new HitSpec(false, false, 30));
        secondOfPair.add(new HitSpec(false, false, 11));
        secondOfPair.add(new HitSpec(true, false, 30));
        ret.add(new Object[]{"Both ends aligned, secondary filtered on first end", firstOfPair, secondOfPair, 2, 2, 3, 40});

        firstOfPair = new ArrayList<HitSpec>();
        secondOfPair = new ArrayList<HitSpec>();
        firstOfPair.add(new HitSpec(false, true, 9));
        firstOfPair.add(new HitSpec(false, true, 10));
        firstOfPair.add(new HitSpec(true, true, 10));
        secondOfPair.add(new HitSpec(false, true, 30));
        secondOfPair.add(new HitSpec(false, true, 11));
        secondOfPair.add(new HitSpec(true, true, 30));
        ret.add(new Object[]{"Both ends aligned, all filtered", firstOfPair, secondOfPair, null, 1, 1, 0});

        return ret.toArray(new Object[0][]);
    }




    /**
     * Read a single fragment read from a file, and create one or more aligned records for the read pair based on
     * the lists, merge with the original read, and assert expected results.
     * @param description
     * @param hitSpecs List that describes the aligned SAMRecords to create.
     * @param expectedPrimaryHitIndex Expected value for the HI tag in the primary alignment in the merged output.
     * @param expectedNumReads Expected number of SAMRecords in the merged output.
     * @param expectedPrimaryMapq Mapq of both ends of primary alignment in the merged output.
     * @throws Exception
     */
    @Test(dataProvider = "testFragmentMultiHitWithFilteringTestCases")
    public void testFragmentMultiHitWithFiltering(String description, List<HitSpec> hitSpecs,
                                                Integer expectedPrimaryHitIndex, int expectedNumReads,
                                                int expectedPrimaryMapq) throws Exception {

        // Create the aligned file by copying bases, quals, readname from the unmapped read, and conforming to each HitSpec.
        final File unmappedSam = new File(TESTDATA_DIR, "multihit.filter.fragment.unmapped.sam");
        SAMRecordIterator unmappedSamFileIterator = new SAMFileReader(unmappedSam).iterator();
        final SAMRecord unmappedRec = unmappedSamFileIterator.next();
        unmappedSamFileIterator.close();
        final File alignedSam = File.createTempFile("aligned.", ".sam");
        alignedSam.deleteOnExit();
        final SAMFileHeader alignedHeader = new SAMFileHeader();
        alignedHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        alignedHeader.setSequenceDictionary(new SAMFileReader(sequenceDict).getFileHeader().getSequenceDictionary());
        final SAMFileWriter alignedWriter = new SAMFileWriterFactory().makeSAMWriter(alignedHeader, true, alignedSam);
        for (int i = 0; i < hitSpecs.size(); ++i) {
            final HitSpec hitSpec = hitSpecs.get(i);
            final SAMRecord mappedRec = makeRead(alignedHeader, unmappedRec, hitSpec, i);
            if (mappedRec != null) {
                alignedWriter.addAlignment(mappedRec);
            }
        }
        alignedWriter.close();

        // Merge aligned file with original unmapped file.
        final File mergedSam = File.createTempFile("merged.", ".sam");
        mergedSam.deleteOnExit();

        MergeBamAlignment merger = new MergeBamAlignment();
        merger.UNMAPPED_BAM = unmappedSam;
        merger.ALIGNED_BAM = Arrays.asList(alignedSam);
        merger.REFERENCE_SEQUENCE = fasta;
        merger.ALIGNED_READS_ONLY = false;
        merger.PROGRAM_RECORD_ID = "0";
        merger.PROGRAM_GROUP_VERSION = "1.0";
        merger.PROGRAM_GROUP_COMMAND_LINE = "align!";
        merger.PROGRAM_GROUP_NAME = "myAligner";
        merger.PAIRED_RUN = true;
        merger.REFERENCE_SEQUENCE = fasta;
        merger.OUTPUT = mergedSam;
        merger.PAIRED_RUN = false;
        merger.doWork();

        assertSamValid(mergedSam);

        // Tally metrics and check for agreement with expected.
        final SAMFileReader mergedReader = new SAMFileReader(mergedSam);
        int numReads = 0;
        Integer primaryHitIndex = null;
        int primaryMapq = 0;
        for (final SAMRecord rec : mergedReader) {
            ++numReads;
            if (!rec.getNotPrimaryAlignmentFlag()  && !rec.getReadUnmappedFlag()) {
                final Integer hitIndex = rec.getIntegerAttribute(SAMTag.HI.name());
                int newHitIndex = (hitIndex == null? -1: hitIndex);
                Assert.assertNull(primaryHitIndex);
                primaryHitIndex = newHitIndex;
                primaryMapq = rec.getMappingQuality();
            }
        }
        Assert.assertEquals(numReads, expectedNumReads);
        Assert.assertEquals(primaryHitIndex, expectedPrimaryHitIndex);
        Assert.assertEquals(primaryMapq, expectedPrimaryMapq);
    }

    @DataProvider(name="testFragmentMultiHitWithFilteringTestCases")
    public Object[][] testFragmentMultiHitWithFilteringTestCases() {
        ArrayList<Object[]> ret = new ArrayList<Object[]>();
        List<HitSpec> hitSpecs;

        hitSpecs = new ArrayList<HitSpec>();
        hitSpecs.add(new HitSpec(false, true, 10));
        hitSpecs.add(new HitSpec(true, false, 8));
        hitSpecs.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"One secondary filtered", hitSpecs, 0, 2, 8});

        hitSpecs = new ArrayList<HitSpec>();
        hitSpecs.add(new HitSpec(false, false, 10));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Primary filtered", hitSpecs, 1, 2, 11});

        hitSpecs = new ArrayList<HitSpec>();
        hitSpecs.add(new HitSpec(false, false, 11));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, false, 11));
        ret.add(new Object[]{"Primary filtered, two secondaries with identical mapq", hitSpecs, 1, 2, 11});

        hitSpecs = new ArrayList<HitSpec>();
        hitSpecs.add(new HitSpec(false, true, 10));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, false, 9));
        ret.add(new Object[]{"Primary and one secondary filtered", hitSpecs, -1, 1, 9});

        hitSpecs = new ArrayList<HitSpec>();
        hitSpecs.add(new HitSpec(false, true, 10));
        hitSpecs.add(new HitSpec(true, true, 8));
        hitSpecs.add(new HitSpec(false, true, 9));
        ret.add(new Object[]{"All filtered", hitSpecs, null, 1, 0});

        return ret.toArray(new Object[0][]);
    }
}
