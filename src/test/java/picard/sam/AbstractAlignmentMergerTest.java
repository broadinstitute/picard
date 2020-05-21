package picard.sam;

import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.StringUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.cmdline.argumentcollections.RequiredReferenceArgumentCollection;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * Tests related to code in AbstractAlignmentMerger
 */
public class AbstractAlignmentMergerTest extends CommandLineProgramTest {

    @DataProvider(name = "overlapReadData")
    public Object[][] overlapReadData() {
        // The spaces here are deliberate to illustrate the region the two default reads match
        final String default120LongR1Bases =                     "ATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCAGATTCTCCTGTCAGTTTGC";
        final String default120LongR2Bases = "CGTTGGCAATGCCGGGCACAATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACC";

        final String default110LongR1Bases =           "ATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACCAGATTCTCCT";
        final String default110LongR2Bases = "GCCGGGCACAATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACC";

        final String sharedBases = "ATCACACCAGTGTCTGCGTTCACAGCAGGCATCATCAGTAGCCTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAGGCTGCAGTGGCTGACC";

        final String default120LongR1ClippedBases = "AGATTCTCCTGTCAGTTTGC";
        final String default120LongR2ClippedBases = "TGTGCCCGGCATTGCCAACG";

        final String default110LongR1ClippedBases = "AGATTCTCCT";
        final String default110LongR2ClippedBases = "TGTGCCCGGC";

        final String default120LongR1BaseQualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFF.FFF.FFF";
        final String default120LongR2BaseQualities ="FFFFFF.FFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
        final String default110LongR1BaseQualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFF.FFF";
        final String default110LongR2BaseQualities = "FFFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";

        final String sharedQualities = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";

        // The use of the reverse method here is so that it is easy to match by eye the base qualities above.
        final String r1ClippedQualities10 = "FF.FFF.FFF";
        final String r2ClippedQualities10 = new StringBuilder("FFFFFF.FFF").reverse().toString();
        final String r1ClippedQualities20 = "FFFFFFFF.FFF.FFF.FFF";
        final String r2ClippedQualities20 = new StringBuilder("FFFFFF.FFFFF.FFFFFFF").reverse().toString();

        return new Object[][] {
                {110, 100, 200, "110M", "110M", false, true, 100, 200, "110M", "110M", false,
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Non overlapping reads

                {110, 100, 200, "110M", "110M", false, true, 100, 200, "110M", "110M", true,
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null},

                {110, 100, 90, "110M", "110M", false, true, 100, 100, "100M10S", "10S100M", false,
                        default110LongR1Bases, default110LongR2Bases, default110LongR1Bases, default110LongR2Bases, null, null,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, default110LongR1BaseQualities, default110LongR2BaseQualities, null, null}, // Basic overlapped read

                {110, 100, 90, "110M", "110M", false, true, 100, 100, "100M10H", "10H100M", true,
                        default110LongR1Bases, default110LongR2Bases, sharedBases, sharedBases, default110LongR1ClippedBases, default110LongR2ClippedBases,
                        default110LongR1BaseQualities, default110LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities10, r2ClippedQualities10},

                {120, 100, 95, "110M10S5H", "5H15S105M", false, true, 100, 100, "100M20S5H", "5H20S100M", false,
                        default120LongR1Bases, default120LongR2Bases, default120LongR1Bases, default120LongR2Bases, null, null,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, default120LongR1BaseQualities, default120LongR2BaseQualities, null, null}, // Already hard and soft clipped

                {120, 100, 95, "110M10S5H", "5H15S105M", false, true, 100, 100, "100M25H", "25H100M", true,
                        default120LongR1Bases, default120LongR2Bases, sharedBases, sharedBases, default120LongR1ClippedBases, default120LongR2ClippedBases,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities20, r2ClippedQualities20},

                {120, 100, 95, "110M10S", "15S105M", false, true, 100, 100, "100M20S", "20S100M", false,
                        default120LongR1Bases, default120LongR2Bases, default120LongR1Bases, default120LongR2Bases, null, null,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, default120LongR1BaseQualities, default120LongR2BaseQualities, null, null}, // Already soft clipped

                {120, 100, 95, "110M10S", "15S105M", false, true, 100, 100, "100M20H", "20H100M", true,
                        default120LongR1Bases, default120LongR2Bases, sharedBases, sharedBases, default120LongR1ClippedBases, default120LongR2ClippedBases,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities20, r2ClippedQualities20},

                {120, 100, 95, "95M25S", "15S105M", false, true, 100, 100, "95M5S20H", "20H100M", true,
                        default120LongR1Bases, default120LongR2Bases, sharedBases, sharedBases, default120LongR1ClippedBases, default120LongR2ClippedBases,
                        default120LongR1BaseQualities, default120LongR2BaseQualities, sharedQualities, sharedQualities, r1ClippedQualities20, r2ClippedQualities20}

        };
    }

    @Test(dataProvider = "overlapReadData")
    public void testOverlappedReadClipping(final int readLength, final int start1, final int start2, final String cigar1, final String cigar2,
                                           final boolean strand1, final boolean strand2,
                                           final int r1ExpectedAlignmentStart, final int r2ExpectedAlignmentStart,
                                           final String expectedR1Cigar, final String expectedR2Cigar, final boolean hardClipOverlappingReads,
                                           final String read1Bases, final String read2Bases, final String expectedR1Bases, final String expectedR2Bases,
                                           final String expectedR1ClippedBases, final String expectedR2ClippedBases, final String read1Qualities,
                                           final String read2Qualities, final String expectedR1Qualities, final String expectedR2Qualities,
                                           final String expectedR1ClippedQualities, final String expectedR2ClippedQualities) {

        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(readLength);
        final List<SAMRecord> recs = set.addPair("q1", 0, start1, start2, false, false, cigar1, cigar2, strand1, strand2, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);

        r1.setReadBases(StringUtil.stringToBytes(read1Bases));
        r2.setReadBases(StringUtil.stringToBytes(read2Bases));

        r1.setBaseQualities(SAMUtils.fastqToPhred(read1Qualities));
        r2.setBaseQualities(SAMUtils.fastqToPhred(read2Qualities));

        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2, hardClipOverlappingReads);
        Assert.assertEquals(r1.getAlignmentStart(), r1ExpectedAlignmentStart);
        Assert.assertEquals(r1.getCigarString(), expectedR1Cigar);
        Assert.assertEquals(r2.getAlignmentStart(), r2ExpectedAlignmentStart);
        Assert.assertEquals(r2.getCigarString(), expectedR2Cigar);
        Assert.assertEquals(r1.getReadString(), expectedR1Bases);
        Assert.assertEquals(r2.getReadString(), expectedR2Bases);
        Assert.assertEquals(SAMUtils.phredToFastq(r1.getBaseQualities()), expectedR1Qualities);
        Assert.assertEquals(SAMUtils.phredToFastq(r2.getBaseQualities()), expectedR2Qualities);

        Assert.assertEquals(r1.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASES_TAG), expectedR1ClippedBases);
        Assert.assertEquals(r2.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASES_TAG), expectedR2ClippedBases);

        Assert.assertEquals(r1.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASE_QUALITIES_TAG), expectedR1ClippedQualities);
        Assert.assertEquals(r2.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASE_QUALITIES_TAG), expectedR2ClippedQualities);


        // Swap first and second read to ensure logic is correct for both F1R2 and F2R1
        final SAMRecordSetBuilder setSwapped = new SAMRecordSetBuilder();
        setSwapped.setReadLength(readLength);
        final List<SAMRecord> recsSwapped = set.addPair("q1", 0, start2, start1, false, false, cigar2, cigar1, strand2, strand1, 30);
        final SAMRecord r1Swapped = recsSwapped.get(0);
        final SAMRecord r2Swapped = recsSwapped.get(1);

        r1Swapped.setReadBases(StringUtil.stringToBytes(read2Bases));
        r2Swapped.setReadBases(StringUtil.stringToBytes(read1Bases));

        r1Swapped.setBaseQualities(SAMUtils.fastqToPhred(read2Qualities));
        r2Swapped.setBaseQualities(SAMUtils.fastqToPhred(read1Qualities));

        AbstractAlignmentMerger.clipForOverlappingReads(r1Swapped, r2Swapped, hardClipOverlappingReads);
        Assert.assertEquals(r1Swapped.getAlignmentStart(), r2ExpectedAlignmentStart);
        Assert.assertEquals(r1Swapped.getCigarString(), expectedR2Cigar);
        Assert.assertEquals(r2Swapped.getAlignmentStart(), r1ExpectedAlignmentStart);
        Assert.assertEquals(r2Swapped.getCigarString(), expectedR1Cigar);
        Assert.assertEquals(r1Swapped.getReadString(), expectedR2Bases);
        Assert.assertEquals(r2Swapped.getReadString(), expectedR1Bases);
        Assert.assertEquals(SAMUtils.phredToFastq(r1Swapped.getBaseQualities()), expectedR2Qualities);
        Assert.assertEquals(SAMUtils.phredToFastq(r2Swapped.getBaseQualities()), expectedR1Qualities);

        Assert.assertEquals(r1Swapped.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASES_TAG), expectedR2ClippedBases);
        Assert.assertEquals(r2Swapped.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASES_TAG), expectedR1ClippedBases);

        Assert.assertEquals(r1Swapped.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASE_QUALITIES_TAG), expectedR2ClippedQualities);
        Assert.assertEquals(r2Swapped.getAttribute(AbstractAlignmentMerger.HARD_CLIPPED_BASE_QUALITIES_TAG), expectedR1ClippedQualities);

    }




    @Test
    public void testUnmapBacterialContamination() throws IOException {
        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder(true, SAMFileHeader.SortOrder.queryname);
        final SAMFileHeader header = builder.getHeader();
        final SAMFileHeader.SortOrder sortOrder = header.getSortOrder();
        final SAMFileHeader newHeader = SAMRecordSetBuilder.makeDefaultHeader(sortOrder, 100000,true);
        builder.setHeader(newHeader);

        final File reference = File.createTempFile("reference",".fasta");
        reference.deleteOnExit();

        builder.writeRandomReference(reference.toPath());

        builder.addPair("overlappingpair", 0,500,500, false,false,"20S20M60S","20S20M60M",true,false,45);
        builder.addPair("overlappingpairFirstAligned", 0,500,500, false,true,"20S20M60S",null,true,false,45);
        builder.addPair("overlappingpairSecondAligned", 0,500,500, true,false,null,"20S20M60S",true,false,45);
        builder.addPair("overlappingpairFirstAlignedB", 0,500,500, false,true,"20S20M60S",null,false,true,45);
        builder.addPair("overlappingpairSecondAlignedB", 0,500,500, true,false,null,"20S20M60S",false,true,45);

//        builder.addFrag("frag",1,500,false,false,"20S20M60S",null, 45);
//        builder.addFrag("frag2",1,500,true,false,"20S20M60S",null, 45);
//        builder.addFrag("frag3",1,500,false,false,"20S20M60S",null, 45);
//        builder.addFrag("frag4",1,500,true,false,"20S20M60S",null, 45);

        final File file = newTempSamFile("aligned");

        try (SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(builder.getHeader(), true, file, null)) {
            builder.getRecords().forEach(writer::addAlignment);
        }

        final RevertSam revertSam = new RevertSam();

        revertSam.INPUT = file;
        final File fileUnaligned = newTempSamFile("unaligned");
        revertSam.OUTPUT = fileUnaligned;

        revertSam.SANITIZE = false;
        revertSam.REMOVE_ALIGNMENT_INFORMATION=true;
        revertSam.REMOVE_DUPLICATE_INFORMATION=true;

        revertSam.SORT_ORDER = SAMFileHeader.SortOrder.queryname;

        Assert.assertEquals(revertSam.doWork(),0);

        MergeBamAlignment mergeBamAlignment = new MergeBamAlignment();

        mergeBamAlignment.ALIGNED_BAM = Collections.singletonList(file);
        mergeBamAlignment.UNMAPPED_BAM = fileUnaligned;
        mergeBamAlignment.UNMAP_CONTAMINANT_READS = true;

        //yuck!
        final RequiredReferenceArgumentCollection requiredReferenceArgumentCollection = new RequiredReferenceArgumentCollection();
        requiredReferenceArgumentCollection.REFERENCE_SEQUENCE = reference;
        mergeBamAlignment.referenceSequence = requiredReferenceArgumentCollection;

        final File fileMerged = newTempSamFile("merged");

        mergeBamAlignment.OUTPUT = fileMerged;

        // merge file with itself.
        Assert.assertEquals(mergeBamAlignment.doWork(),0);

        // check that all reads have been unmapped due to bacterial contamination as needed.
        try (SamReader mergedReader = SamReaderFactory.makeDefault().open(fileMerged)) {
            for (SAMRecord mergedRecord : mergedReader) {
                Assert.assertTrue(mergedRecord.getReadUnmappedFlag(), mergedRecord.toString());
                Assert.assertTrue(!mergedRecord.getReadPairedFlag() || mergedRecord.getMateUnmappedFlag(), mergedRecord.toString());
            }
        }
    }

    @Override
    public String getCommandLineProgramName() {
        return this.getClass().getSimpleName();
    }


    private static File newTempSamFile(final String filename) throws IOException {
        final File file = File.createTempFile(filename,".sam");
        file.deleteOnExit();
        return file;
    }
    @DataProvider(name = "readPositionIgnoringSoftClips")
    public Object[][] readPositionIgnoringSoftClips() {


        return new Object[][] {
                {"26S58M62S", 55533688, 55533827, 0}, // This is from the read that made us aware of a bug
                {"26S58M62S", 55533688, 55533665, 4},
                {"26S58M62S", 55533688, 55533660, 0}, // Before soft clip
                {"10S100M2S", 5, 10, 16},
                {"10S100M2S", 5, 3, 9},
                {"10S100M2S", 10, 12, 13},
                {"10S100M2S", 5, 107, 0},
        };
    }
    @Test(dataProvider = "readPositionIgnoringSoftClips")
    public void testGetReadPositionIgnoringSoftClips(final String cigarString, final int startPosition, final int queryPosition, final int expectedReadPosititon) {
        final SAMFileHeader newHeader = SAMRecordSetBuilder.makeDefaultHeader(SAMFileHeader.SortOrder.queryname, 100000,false);
        final SAMRecord rec = new SAMRecord(newHeader);

        rec.setCigarString(cigarString);
        rec.setAlignmentStart(startPosition);

        final int readPosition = AbstractAlignmentMerger.getReadPositionAtReferencePositionIgnoreSoftClips(rec, queryPosition);

        Assert.assertEquals(readPosition, expectedReadPosititon);
    }
}

