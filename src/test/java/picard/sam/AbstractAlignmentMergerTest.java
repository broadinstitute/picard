package picard.sam;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang.ArrayUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.cmdline.argumentcollections.RequiredReferenceArgumentCollection;
import picard.illumina.CustomAdapterPair;
import picard.illumina.parser.ReadDescriptor;
import picard.illumina.parser.ReadStructure;
import picard.illumina.parser.ReadType;
import picard.util.AdapterPair;
import picard.util.IlluminaUtil;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Tests related to code in AbstractAlignmentMerger
 */
public class AbstractAlignmentMergerTest extends CommandLineProgramTest {
    @Test public void tesOverlappedReadClippingWithNonOverlappedReads() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(110);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 200, false, false, "110M", "110M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "110M");
        Assert.assertEquals(r2.getAlignmentStart(), 200);
        Assert.assertEquals(r2.getCigarString(), "110M");
    }

    @Test
    public void testOverlappedReadHardClippingWithNonOverlappedReads() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(110);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 200, false, false, "110M", "110M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2, CigarOperator.HARD_CLIP, null, null, Collections.emptyList());
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "110M");
        Assert.assertEquals(r2.getAlignmentStart(), 200);
        Assert.assertEquals(r2.getCigarString(), "110M");
        Assert.assertEquals(r1.getReadLength(), 110);
        Assert.assertEquals(r2.getReadLength(), 110);
    }

    @Test public void testBasicOverlappedReadClipping() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(110);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 90, false, false, "110M", "110M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "100M10S");
        Assert.assertEquals(r2.getAlignmentStart(), 100);
        Assert.assertEquals(r2.getCigarString(), "10S100M");
    }

    @DataProvider(name = "hardClippingDataProvider")
    public Object [][] getHardClippingData() {


        final List<Object[]> ret = new ArrayList<>();
        final List<AdapterPair> illuminaAdapters = Arrays.asList(IlluminaUtil.IlluminaAdapterPair.values());
        final byte[] templateBases = StringUtil.stringToBytes("ACTGCATGCTAGCTTAGGACAGATACGATAGCTAGACAGACATAATTTAGCGGATGACATTCGGACAGATCGGACGAGCTAGACAGACTGAGACAGCTAGCAGATCGAGG");
        final byte[] templateBasesRC = Arrays.copyOf(templateBases, templateBases.length);
        SequenceUtil.reverseComplement(templateBasesRC);
        //All Illumina barcode combinations should work
        
        for (int nAdapterBases = 1; nAdapterBases <= 10; nAdapterBases+=3) {
            final int readLength = templateBases.length + nAdapterBases;
            final String cigarF = readLength + "M";
            final String cigarR = readLength + "M";
            final String expectedCigarF = templateBases.length + "M" + nAdapterBases + "H";
            final String expectedCigarR = nAdapterBases + "H" + templateBases.length + "M";

            for (final AdapterPair adapterPair : illuminaAdapters) {
                SequenceUtil.reverseComplement(templateBasesRC);
                ret.add(new Object[]{cigarF, cigarR, buildReadBases(templateBases, adapterPair.get3PrimeAdapterBytesInReadOrder(), readLength, false),
                        buildReadBases(templateBasesRC, adapterPair.get5PrimeAdapterBytesInReadOrder(), readLength, true),
                        false, true, 100, 100 - nAdapterBases, illuminaAdapters, expectedCigarF, expectedCigarR, 100, 100, null, null}); //F1R2
                ret.add(new Object[]{cigarR, cigarF, buildReadBases(templateBases, adapterPair.get3PrimeAdapterBytesInReadOrder(), readLength, true),
                        buildReadBases(templateBasesRC, adapterPair.get5PrimeAdapterBytesInReadOrder(), readLength, false),
                        true, false, 100 - nAdapterBases, 100, illuminaAdapters, expectedCigarR, expectedCigarF, 100, 100, null, null}); //F2R1
            }

            final String expectedCigarSoftF = templateBases.length + "M" + nAdapterBases + "S";
            final String expectedCigarSoftR = nAdapterBases +"S" + templateBases.length +"M";

            //3' of one barcode, 5' of another, should only softclip
            ret.add(new Object[] {cigarF, cigarR, buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), readLength, false),
                    buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.NEXTERA_V1.get5PrimeAdapterBytesInReadOrder(), readLength, true),
                    false, true, 100, 100 - nAdapterBases, illuminaAdapters, expectedCigarSoftF, expectedCigarSoftR, 100, 100, null, null}); //F1R2
            ret.add(new Object[] {cigarR, cigarF, buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), readLength, true),
                    buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.NEXTERA_V1.get5PrimeAdapterBytesInReadOrder(), readLength, false),
                    true, false, 100 - nAdapterBases, 100, illuminaAdapters, expectedCigarSoftR, expectedCigarSoftF, 100, 100, null, null}); //F2R1

        }

        //already soft-clipped
        ret.add(new Object[] {"118M2S", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "110M10H", "10H110M", 100, 100, null, null});
        //already soft-clipped more than adapters
        ret.add(new Object[] {"108M12S", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "108M2S10H", "10H110M", 100, 100, null, null});

        //already hard-clipped
        ret.add(new Object[] {"118M2H", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 118, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "110M10H", "10H110M", 100, 100, null, null});

        //soft-clipped beginning
        ret.add(new Object[] {"2S118M", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 102, 90, illuminaAdapters, "2S108M10H", "10H110M", 102, 100, null, null});


        //hard-clipped beginning
        ret.add(new Object[]{"2H118M", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 118, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 102, 90, illuminaAdapters, "2H110M8H", "10H110M", 102, 100, null, null});

        //insertion in reads
        ret.add(new Object[]{"50M2I68M", "50M2I68M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "50M2I58M10H", "10H40M2I68M", 100, 100, null, null});

        //insertion in adapter portion of read
        ret.add(new Object[]{"114M2I4M", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "110M10H", "10H110M", 100, 100, null, null});

        //deletion in reads
        ret.add(new Object[]{"50M2D70M", "50M2D70M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "50M2D60M10H", "10H40M2D70M", 100, 100, null, null});

        //deletion in adapter portion of read
        ret.add(new Object[]{"114M2D6M", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                false, true, 100, 90, illuminaAdapters, "110M10H", "10H110M", 100, 100, null, null});

        //Ns in adapter
        final CustomAdapterPair customAdapterPair = new CustomAdapterPair("GTGCTTGCANNN", "NNNNAGTCGATTGC");
        ret.add(new Object[] {"120M", "120M", buildReadBases(templateBases, StringUtil.stringToBytes("ACGTAGTCGATTGC"), 120, false),
                buildReadBases(templateBasesRC, StringUtil.stringToBytes("ACGTGCAAGCAC"), 120, true),
                false, true, 100, 90, Collections.singletonList(customAdapterPair), "110M10H", "10H110M", 100, 100, null, null});

        //read structures and UMIs
        for(int umiLength = 0; umiLength<=6; umiLength+=3) {
            for (int gapBeforeUMI = 0; gapBeforeUMI<3; gapBeforeUMI++) {
                for (int gapAfterUMI = 0; gapAfterUMI<3 && gapAfterUMI<umiLength; gapAfterUMI++) {
                    final int totalBasesToRemove = umiLength + gapAfterUMI + gapBeforeUMI;
                    final byte[] templateBasesClipped = Arrays.copyOfRange(templateBases, totalBasesToRemove, templateBases.length);
                    final byte[] templateBasesRCClipped = Arrays.copyOfRange(templateBasesRC, totalBasesToRemove, templateBasesRC.length);


                    final byte[] umi1 = Arrays.copyOfRange(templateBases, gapBeforeUMI, gapBeforeUMI + umiLength);
                    final byte[] umi2 = Arrays.copyOfRange(templateBasesRC, gapBeforeUMI, gapBeforeUMI + umiLength);

                    final String umi1String = StringUtil.bytesToString(umi1);
                    final String umi2String = StringUtil.bytesToString(umi2);

                    SequenceUtil.reverseComplement(umi1);
                    SequenceUtil.reverseComplement(umi2);

                    final byte[] templateBasesWithUMI = ArrayUtils.addAll(templateBasesClipped, umi1);
                    final byte[] templateBasesRCWithUMI = ArrayUtils.addAll(templateBasesRCClipped, umi2);

                    final List<ReadDescriptor> readDescriptors = new ArrayList<>();
                    if (gapBeforeUMI > 0) {
                        readDescriptors.add(new ReadDescriptor(gapBeforeUMI, ReadType.Skip));
                    }
                    if (umiLength > 0) {
                        readDescriptors.add(new ReadDescriptor(umiLength, ReadType.MolecularIndex));
                    }
                    if (gapAfterUMI > 0) {
                        readDescriptors.add(new ReadDescriptor(gapAfterUMI, ReadType.Skip));
                    }
                    readDescriptors.add(new ReadDescriptor(110 - totalBasesToRemove, ReadType.Template));

                    final ReadStructure readStructure = new ReadStructure(readDescriptors);

                    ret.add(new Object[] {"120M", "120M", buildReadBases(templateBases, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapterBytesInReadOrder(), 120, false),
                            buildReadBases(templateBasesRC, IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapterBytesInReadOrder(), 120, true),
                            false, true, 100, 90 - totalBasesToRemove, illuminaAdapters, (110 - totalBasesToRemove) + "M" + (10 + totalBasesToRemove) + "H", (10 + totalBasesToRemove) + "H" + (110 - totalBasesToRemove) + "M", 100, 100,
                            umi1String + "-" + umi2String, readStructure});

                }
            }
        }

        return ret.toArray(new Object[][]{});
    }

    private byte[] buildReadBases(final byte[] templateBasesReadOrder, final byte[] adapterBasesReadOrder, final int readLength, final boolean negativeStrand) {
        final byte[] bases = ArrayUtils.addAll(templateBasesReadOrder, adapterBasesReadOrder);
        final byte[] readBases = Arrays.copyOf(bases, readLength);
        if (negativeStrand) {
            SequenceUtil.reverseComplement(readBases);
        }
        return readBases;
    }

    @Test (dataProvider = "hardClippingDataProvider")
    public void testOverlappedReadHardClipping(final String originalCigar1, final String originalCigar2, final byte[] read1Bases, final byte[] read2Bases, final boolean strand1, final boolean strand2,
                                               final int start1, final int start2, final List<AdapterPair> adapters, final String expectedCigar1, final String expectedCigar2, final int expectedStart1, final int expectedStart2,
                                               final String umiTag, final ReadStructure readStructure) {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        final List<SAMRecord> recs = set.addPair("q1", 0, start1, start2, false, false, originalCigar1, originalCigar2, strand1, strand2, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);

        r1.setReadBases(read1Bases);
        r2.setReadBases(read2Bases);

        if (umiTag != null) {
            r1.setAttribute(SAMTag.RX.toString(), umiTag);
            r2.setAttribute(SAMTag.RX.toString(), umiTag);
        }

        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2, CigarOperator.HARD_CLIP, readStructure, readStructure, adapters);
        Assert.assertEquals(r1.getAlignmentStart(), expectedStart1);
        Assert.assertEquals(r1.getCigarString(), expectedCigar1);

        Assert.assertEquals(r2.getAlignmentStart(), expectedStart2);
        Assert.assertEquals(r2.getCigarString(), expectedCigar2);
    }

    @DataProvider(name = "getReadPosToClipFromDataProvider")
    public Object[][] getReadPosToClipFromData() {
        return new Object[][] {
                {"120M", false, 450, 502, 53},
                {"120M", true, 450, 502, 68},
                {"120M", false, 450, 580, -1},
                {"120M", false, 450, 440, -1},
                {"120M", true, 450, 580, -1},
                {"120M", true, 450, 440, -1},

                {"100M3I100M", false, 300, 425, 129},
                {"100M3I100M", true, 300, 425, 75},

                {"100M3D100M", false, 300, 425, 123},
                {"100M3D100M", true, 300, 425, 78},

                {"100M16S", false, 300, 410, 111},
                {"100M16S", true, 300, 400, 16},
                {"16S100M", false, 316, 410, 111},
                {"16S100M", true, 316, 310, 106},

                {"100M16H", false, 300, 410, -1},
                {"100M16H", true, 300, 410, -1},
                {"100M16H", true, 300, 390, 26},
                {"16H100M", false, 316, 310, -1},
                {"16H100M", true, 316, 310, -1},
                {"16H100M", true, 316, 350, 66}
        };
    }

    @Test(dataProvider = "getReadPosToClipFromDataProvider")
    public void testGetReadPositionToClipFrom(final String cigarString, final boolean negativeStrand, final int start, final int refPosToClipFrom, final int expectedReadPosToCLipFrom) {
        final SAMFileHeader header = new SAMFileHeader();
        final SAMRecord rec = new SAMRecord(header);
        rec.setCigarString(cigarString);
        rec.setReadNegativeStrandFlag(negativeStrand);
        rec.setAlignmentStart(start);

        final int readPositionToClipFrom = AbstractAlignmentMerger.getReadPositionToClipFrom(rec, refPosToClipFrom);
        Assert.assertEquals(readPositionToClipFrom, expectedReadPosToCLipFrom);
    }

    @Test public void testOverlappedReadClippingWithExistingSoftClipping() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(120);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 95, false, false, "110M10S", "15S105M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "100M20S");
        Assert.assertEquals(r2.getAlignmentStart(), 100);
        Assert.assertEquals(r2.getCigarString(), "20S100M");
    }

    @Test public void testOverlappedReadClippingWithExistingSoftClippingAndHardClipping() {
        final SAMRecordSetBuilder set = new SAMRecordSetBuilder();
        set.setReadLength(120);
        final List<SAMRecord> recs = set.addPair("q1", 0, 100, 95, false, false, "110M10S5H", "5H15S105M", false, true, 30);
        final SAMRecord r1 = recs.get(0);
        final SAMRecord r2 = recs.get(1);
        AbstractAlignmentMerger.clipForOverlappingReads(r1, r2);
        Assert.assertEquals(r1.getAlignmentStart(), 100);
        Assert.assertEquals(r1.getCigarString(), "100M20S5H");
        Assert.assertEquals(r2.getAlignmentStart(), 100);
        Assert.assertEquals(r2.getCigarString(), "5H20S100M");
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
}

