package picard.sam.BamErrorMetric;

import htsjdk.samtools.*;
import htsjdk.samtools.util.SamLocusIterator;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * A slew of unit tests for the various stratifiers
 */
public class ReadBaseStratificationTest {

    @BeforeClass
    public void setup() {
        ReadBaseStratification.setLongHomopolymer(6);
    }

    @DataProvider
    public Object[][] baseStratifierData() {
        return new Object[][]{

                {-1, true, null, ReadBaseStratification.currentReadBaseStratifier},
                {0, true, 'C', ReadBaseStratification.currentReadBaseStratifier},
                {1, true, 'A', ReadBaseStratification.currentReadBaseStratifier},
                {2, true, 'T', ReadBaseStratification.currentReadBaseStratifier},
                {3, true, 'G', ReadBaseStratification.currentReadBaseStratifier},
                {4, true, 'G', ReadBaseStratification.currentReadBaseStratifier},
                {5, true, 'G', ReadBaseStratification.currentReadBaseStratifier},
                {6, true, 'G', ReadBaseStratification.currentReadBaseStratifier},
                {7, true, 'A', ReadBaseStratification.currentReadBaseStratifier},
                {8, true, null, ReadBaseStratification.currentReadBaseStratifier},

                {-1, false, null, ReadBaseStratification.currentReadBaseStratifier},
                {0, false, 'G', ReadBaseStratification.currentReadBaseStratifier},
                {1, false, 'T', ReadBaseStratification.currentReadBaseStratifier},
                {2, false, 'A', ReadBaseStratification.currentReadBaseStratifier},
                {3, false, 'C', ReadBaseStratification.currentReadBaseStratifier},
                {4, false, 'C', ReadBaseStratification.currentReadBaseStratifier},
                {5, false, 'C', ReadBaseStratification.currentReadBaseStratifier},
                {6, false, 'C', ReadBaseStratification.currentReadBaseStratifier},
                {7, false, 'T', ReadBaseStratification.currentReadBaseStratifier},
                {8, false, null, ReadBaseStratification.currentReadBaseStratifier},

                {0, true, null, ReadBaseStratification.previousReadBaseStratifier},
                {1, true, 'C', ReadBaseStratification.previousReadBaseStratifier},
                {2, true, 'A', ReadBaseStratification.previousReadBaseStratifier},
                {3, true, 'T', ReadBaseStratification.previousReadBaseStratifier},
                {4, true, 'G', ReadBaseStratification.previousReadBaseStratifier},
                {5, true, 'G', ReadBaseStratification.previousReadBaseStratifier},
                {6, true, 'G', ReadBaseStratification.previousReadBaseStratifier},
                {7, true, 'G', ReadBaseStratification.previousReadBaseStratifier},
                {8, true, 'A', ReadBaseStratification.previousReadBaseStratifier},
                {9, true, null, ReadBaseStratification.previousReadBaseStratifier},

                {-2, false, null, ReadBaseStratification.previousReadBaseStratifier},
                {-1, false, 'G', ReadBaseStratification.previousReadBaseStratifier},
                {0, false, 'T', ReadBaseStratification.previousReadBaseStratifier},
                {1, false, 'A', ReadBaseStratification.previousReadBaseStratifier},
                {2, false, 'C', ReadBaseStratification.previousReadBaseStratifier},
                {3, false, 'C', ReadBaseStratification.previousReadBaseStratifier},
                {4, false, 'C', ReadBaseStratification.previousReadBaseStratifier},
                {5, false, 'C', ReadBaseStratification.previousReadBaseStratifier},
                {6, false, 'T', ReadBaseStratification.previousReadBaseStratifier},
                {7, false, null, ReadBaseStratification.previousReadBaseStratifier},

                {-2, true, null, ReadBaseStratification.nextReadBaseStratifier},
                {-1, true, 'C', ReadBaseStratification.nextReadBaseStratifier},
                {0, true, 'A', ReadBaseStratification.nextReadBaseStratifier},
                {1, true, 'T', ReadBaseStratification.nextReadBaseStratifier},
                {2, true, 'G', ReadBaseStratification.nextReadBaseStratifier},
                {3, true, 'G', ReadBaseStratification.nextReadBaseStratifier},
                {4, true, 'G', ReadBaseStratification.nextReadBaseStratifier},
                {5, true, 'G', ReadBaseStratification.nextReadBaseStratifier},
                {6, true, 'A', ReadBaseStratification.nextReadBaseStratifier},
                {7, true, null, ReadBaseStratification.nextReadBaseStratifier},

                {0, false, null, ReadBaseStratification.nextReadBaseStratifier},
                {1, false, 'G', ReadBaseStratification.nextReadBaseStratifier},
                {2, false, 'T', ReadBaseStratification.nextReadBaseStratifier},
                {3, false, 'A', ReadBaseStratification.nextReadBaseStratifier},
                {4, false, 'C', ReadBaseStratification.nextReadBaseStratifier},
                {5, false, 'C', ReadBaseStratification.nextReadBaseStratifier},
                {6, false, 'C', ReadBaseStratification.nextReadBaseStratifier},
                {7, false, 'C', ReadBaseStratification.nextReadBaseStratifier},
                {8, false, 'T', ReadBaseStratification.nextReadBaseStratifier},
                {9, false, null, ReadBaseStratification.nextReadBaseStratifier},

                {-1, true, null, ReadBaseStratification.referenceBaseStratifier},
                {0, true, 'C', ReadBaseStratification.referenceBaseStratifier},
                {1, true, 'A', ReadBaseStratification.referenceBaseStratifier},
                {2, true, 'T', ReadBaseStratification.referenceBaseStratifier},
                {3, true, 'G', ReadBaseStratification.referenceBaseStratifier},
                {4, true, 'G', ReadBaseStratification.referenceBaseStratifier},
                {5, true, 'G', ReadBaseStratification.referenceBaseStratifier},
                {6, true, 'G', ReadBaseStratification.referenceBaseStratifier},
                {7, true, 'A', ReadBaseStratification.referenceBaseStratifier},
                {8, true, null, ReadBaseStratification.referenceBaseStratifier},

                {-1, false, null, ReadBaseStratification.referenceBaseStratifier},
                {0, false, 'G', ReadBaseStratification.referenceBaseStratifier},
                {1, false, 'T', ReadBaseStratification.referenceBaseStratifier},
                {2, false, 'A', ReadBaseStratification.referenceBaseStratifier},
                {3, false, 'C', ReadBaseStratification.referenceBaseStratifier},
                {4, false, 'C', ReadBaseStratification.referenceBaseStratifier},
                {5, false, 'C', ReadBaseStratification.referenceBaseStratifier},
                {6, false, 'C', ReadBaseStratification.referenceBaseStratifier},
                {7, false, 'T', ReadBaseStratification.referenceBaseStratifier},
                {8, false, null, ReadBaseStratification.referenceBaseStratifier},

                {-1, true, null, ReadBaseStratification.postDiNucleotideStratifier},
                {0, true, new Pair<>('C', 'A'), ReadBaseStratification.postDiNucleotideStratifier},
                {1, true, new Pair<>('A', 'T'), ReadBaseStratification.postDiNucleotideStratifier},
                {2, true, new Pair<>('T', 'G'), ReadBaseStratification.postDiNucleotideStratifier},
                {3, true, new Pair<>('G', 'G'), ReadBaseStratification.postDiNucleotideStratifier},
                {4, true, new Pair<>('G', 'G'), ReadBaseStratification.postDiNucleotideStratifier},
                {5, true, new Pair<>('G', 'G'), ReadBaseStratification.postDiNucleotideStratifier},
                {6, true, new Pair<>('G', 'A'), ReadBaseStratification.postDiNucleotideStratifier},
                {7, true, null, ReadBaseStratification.postDiNucleotideStratifier},

                {0, false, null, ReadBaseStratification.postDiNucleotideStratifier},
                {1, false, new Pair<>('T', 'G'), ReadBaseStratification.postDiNucleotideStratifier},
                {2, false, new Pair<>('A', 'T'), ReadBaseStratification.postDiNucleotideStratifier},
                {3, false, new Pair<>('C', 'A'), ReadBaseStratification.postDiNucleotideStratifier},
                {4, false, new Pair<>('C', 'C'), ReadBaseStratification.postDiNucleotideStratifier},
                {5, false, new Pair<>('C', 'C'), ReadBaseStratification.postDiNucleotideStratifier},
                {6, false, new Pair<>('C', 'C'), ReadBaseStratification.postDiNucleotideStratifier},
                {7, false, new Pair<>('T', 'C'), ReadBaseStratification.postDiNucleotideStratifier},
                {8, false, null, ReadBaseStratification.postDiNucleotideStratifier},

                {0, true, null, ReadBaseStratification.preDiNucleotideStratifier},
                {1, true, new Pair<>('C', 'A'), ReadBaseStratification.preDiNucleotideStratifier},
                {2, true, new Pair<>('A', 'T'), ReadBaseStratification.preDiNucleotideStratifier},
                {3, true, new Pair<>('T', 'G'), ReadBaseStratification.preDiNucleotideStratifier},
                {4, true, new Pair<>('G', 'G'), ReadBaseStratification.preDiNucleotideStratifier},
                {5, true, new Pair<>('G', 'G'), ReadBaseStratification.preDiNucleotideStratifier},
                {6, true, new Pair<>('G', 'G'), ReadBaseStratification.preDiNucleotideStratifier},
                {7, true, new Pair<>('G', 'A'), ReadBaseStratification.preDiNucleotideStratifier},
                {8, true, null, ReadBaseStratification.preDiNucleotideStratifier},

                {-1, false, null, ReadBaseStratification.preDiNucleotideStratifier},
                {0, false, new Pair<>('T', 'G'), ReadBaseStratification.preDiNucleotideStratifier},
                {1, false, new Pair<>('A', 'T'), ReadBaseStratification.preDiNucleotideStratifier},
                {2, false, new Pair<>('C', 'A'), ReadBaseStratification.preDiNucleotideStratifier},
                {3, false, new Pair<>('C', 'C'), ReadBaseStratification.preDiNucleotideStratifier},
                {4, false, new Pair<>('C', 'C'), ReadBaseStratification.preDiNucleotideStratifier},
                {5, false, new Pair<>('C', 'C'), ReadBaseStratification.preDiNucleotideStratifier},
                {6, false, new Pair<>('T', 'C'), ReadBaseStratification.preDiNucleotideStratifier},
                {7, false, null, ReadBaseStratification.preDiNucleotideStratifier},

                {-1, true, null, ReadBaseStratification.homoPolymerLengthStratifier},
                {0, true, 0, ReadBaseStratification.homoPolymerLengthStratifier},
                {1, true, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {2, true, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {3, true, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {4, true, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {5, true, 2, ReadBaseStratification.homoPolymerLengthStratifier},
                {6, true, 3, ReadBaseStratification.homoPolymerLengthStratifier},
                {7, true, 4, ReadBaseStratification.homoPolymerLengthStratifier},
                {8, true, null, ReadBaseStratification.homoPolymerLengthStratifier},

                {-1, false, null, ReadBaseStratification.homoPolymerLengthStratifier},
                {0, false, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {1, false, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {2, false, 4, ReadBaseStratification.homoPolymerLengthStratifier},
                {3, false, 3, ReadBaseStratification.homoPolymerLengthStratifier},
                {4, false, 2, ReadBaseStratification.homoPolymerLengthStratifier},
                {5, false, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {6, false, 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {7, false, 0, ReadBaseStratification.homoPolymerLengthStratifier},
                {8, false, null, ReadBaseStratification.homoPolymerLengthStratifier},

                {0, true, null, ReadBaseStratification.homopolymerStratifier},
                {1, true, new Pair<>( 1, new Pair<>('C', 'A')), ReadBaseStratification.homopolymerStratifier},
                {2, true, new Pair<>( 1, new Pair<>('A', 'T')), ReadBaseStratification.homopolymerStratifier},
                {3, true, new Pair<>( 1, new Pair<>('T', 'G')), ReadBaseStratification.homopolymerStratifier},
                {4, true, new Pair<>( 1, new Pair<>('G', 'G')), ReadBaseStratification.homopolymerStratifier},
                {5, true, new Pair<>( 2, new Pair<>('G', 'G')), ReadBaseStratification.homopolymerStratifier},
                {6, true, new Pair<>( 3, new Pair<>('G', 'G')), ReadBaseStratification.homopolymerStratifier},
                {7, true, new Pair<>( 4, new Pair<>('G', 'A')), ReadBaseStratification.homopolymerStratifier},
                {8, true, null, ReadBaseStratification.homopolymerStratifier},

                {-1, false, null, ReadBaseStratification.homopolymerStratifier},
                {0, false, new Pair<>(1, new Pair<>('T', 'G')), ReadBaseStratification.homopolymerStratifier},
                {1, false, new Pair<>(1, new Pair<>('A', 'T')), ReadBaseStratification.homopolymerStratifier},
                {2, false, new Pair<>(4, new Pair<>('C', 'A')), ReadBaseStratification.homopolymerStratifier},
                {3, false, new Pair<>(3, new Pair<>('C', 'C')), ReadBaseStratification.homopolymerStratifier},
                {4, false, new Pair<>(2, new Pair<>('C', 'C')), ReadBaseStratification.homopolymerStratifier},
                {5, false, new Pair<>(1, new Pair<>('C', 'C')), ReadBaseStratification.homopolymerStratifier},
                {6, false, new Pair<>(1, new Pair<>('T', 'C')), ReadBaseStratification.homopolymerStratifier},
                {7, false, null, ReadBaseStratification.homopolymerStratifier},

                {0, true, null, ReadBaseStratification.oneBasePaddedContextStratifier},
                {1, true, "CAT", ReadBaseStratification.oneBasePaddedContextStratifier},
                {2, true, "ATG", ReadBaseStratification.oneBasePaddedContextStratifier},
                {3, true, "TGG", ReadBaseStratification.oneBasePaddedContextStratifier},
                {4, true, "GGG", ReadBaseStratification.oneBasePaddedContextStratifier},
                {5, true, "GGG", ReadBaseStratification.oneBasePaddedContextStratifier},
                {6, true, "GGA", ReadBaseStratification.oneBasePaddedContextStratifier},
                {7, true, null, ReadBaseStratification.oneBasePaddedContextStratifier},

                {0, false, null, ReadBaseStratification.oneBasePaddedContextStratifier},
                {1, false, "ATG", ReadBaseStratification.oneBasePaddedContextStratifier},
                {2, false, "CAT", ReadBaseStratification.oneBasePaddedContextStratifier},
                {3, false, "CCA", ReadBaseStratification.oneBasePaddedContextStratifier},
                {4, false, "CCC", ReadBaseStratification.oneBasePaddedContextStratifier},
                {5, false, "CCC", ReadBaseStratification.oneBasePaddedContextStratifier},
                {6, false, "TCC", ReadBaseStratification.oneBasePaddedContextStratifier},
                {7, false, null, ReadBaseStratification.oneBasePaddedContextStratifier},

                {0, true, null, ReadBaseStratification.twoBasePaddedContextStratifier},
                {1, true, null, ReadBaseStratification.twoBasePaddedContextStratifier},
                {2, true, "CATGG", ReadBaseStratification.twoBasePaddedContextStratifier},
                {3, true, "ATGGG", ReadBaseStratification.twoBasePaddedContextStratifier},
                {4, true, "TGGGG", ReadBaseStratification.twoBasePaddedContextStratifier},
                {5, true, "GGGGA", ReadBaseStratification.twoBasePaddedContextStratifier},
                {6, true, null, ReadBaseStratification.twoBasePaddedContextStratifier},
                {7, true, null, ReadBaseStratification.twoBasePaddedContextStratifier},

                {0, false, null, ReadBaseStratification.twoBasePaddedContextStratifier},
                {1, false, null, ReadBaseStratification.twoBasePaddedContextStratifier},
                {2, false, "CCATG", ReadBaseStratification.twoBasePaddedContextStratifier},
                {3, false, "CCCAT", ReadBaseStratification.twoBasePaddedContextStratifier},
                {4, false, "CCCCA", ReadBaseStratification.twoBasePaddedContextStratifier},
                {5, false, "TCCCC", ReadBaseStratification.twoBasePaddedContextStratifier},
                {6, false, null, ReadBaseStratification.twoBasePaddedContextStratifier},
                {7, false, null, ReadBaseStratification.twoBasePaddedContextStratifier},

                {0, true, "all", ReadBaseStratification.nonStratifier},
                {1, true, "all", ReadBaseStratification.nonStratifier},
                {2, true, "all", ReadBaseStratification.nonStratifier},
                {3, true, "all", ReadBaseStratification.nonStratifier},
                {4, true, "all", ReadBaseStratification.nonStratifier},
                {5, true, "all", ReadBaseStratification.nonStratifier},
                {6, true, "all", ReadBaseStratification.nonStratifier},
                {7, true, "all", ReadBaseStratification.nonStratifier},
                {0, false, "all", ReadBaseStratification.nonStratifier},
                {1, false, "all", ReadBaseStratification.nonStratifier},
                {2, false, "all", ReadBaseStratification.nonStratifier},
                {3, false, "all", ReadBaseStratification.nonStratifier},
                {4, false, "all", ReadBaseStratification.nonStratifier},
                {5, false, "all", ReadBaseStratification.nonStratifier},
                {6, false, "all", ReadBaseStratification.nonStratifier},
                {7, false, "all", ReadBaseStratification.nonStratifier},

                {7, false, ReadBaseStratification.ReadDirection.NEGATIVE, ReadBaseStratification.readDirectionStratifier},
                {7, true, ReadBaseStratification.ReadDirection.POSITIVE, ReadBaseStratification.readDirectionStratifier},

                {1, false, (Math.round(100.0 * 5 / 8)) / 100.0, ReadBaseStratification.gcContentStratifier},

                {1, false, 93, ReadBaseStratification.flowCellTileStratifier},

                {1, false, "rgID", ReadBaseStratification.readgroupStratifier},

                {1, false, ReadBaseStratification.ReadOrdinality.FIRST, ReadBaseStratification.readOrdinalityStratifier},

                {1, true, ReadBaseStratification.ReadDirection.POSITIVE, ReadBaseStratification.readDirectionStratifier},
                {1, false, ReadBaseStratification.ReadDirection.NEGATIVE, ReadBaseStratification.readDirectionStratifier},

                {1, true, ReadBaseStratification.PairOrientation.F1R2, ReadBaseStratification.readOrientationStratifier},
                {1, false, ReadBaseStratification.PairOrientation.F2R1, ReadBaseStratification.readOrientationStratifier},

                {0, true, ReadBaseStratification.CycleBin.QUINTILE_1, ReadBaseStratification.binnedReadCycleStratifier},
                {1, true, ReadBaseStratification.CycleBin.QUINTILE_2, ReadBaseStratification.binnedReadCycleStratifier},
                {2, true, ReadBaseStratification.CycleBin.QUINTILE_2, ReadBaseStratification.binnedReadCycleStratifier},
                {3, true, ReadBaseStratification.CycleBin.QUINTILE_3, ReadBaseStratification.binnedReadCycleStratifier},
                {4, true, ReadBaseStratification.CycleBin.QUINTILE_4, ReadBaseStratification.binnedReadCycleStratifier},
                {5, true, ReadBaseStratification.CycleBin.QUINTILE_4, ReadBaseStratification.binnedReadCycleStratifier},
                {6, true, ReadBaseStratification.CycleBin.QUINTILE_5, ReadBaseStratification.binnedReadCycleStratifier},
                {7, true, ReadBaseStratification.CycleBin.QUINTILE_5, ReadBaseStratification.binnedReadCycleStratifier},

                {7, false, ReadBaseStratification.CycleBin.QUINTILE_1, ReadBaseStratification.binnedReadCycleStratifier},
                {6, false, ReadBaseStratification.CycleBin.QUINTILE_2, ReadBaseStratification.binnedReadCycleStratifier},
                {5, false, ReadBaseStratification.CycleBin.QUINTILE_2, ReadBaseStratification.binnedReadCycleStratifier},
                {4, false, ReadBaseStratification.CycleBin.QUINTILE_3, ReadBaseStratification.binnedReadCycleStratifier},
                {3, false, ReadBaseStratification.CycleBin.QUINTILE_4, ReadBaseStratification.binnedReadCycleStratifier},
                {2, false, ReadBaseStratification.CycleBin.QUINTILE_4, ReadBaseStratification.binnedReadCycleStratifier},
                {1, false, ReadBaseStratification.CycleBin.QUINTILE_5, ReadBaseStratification.binnedReadCycleStratifier},
                {0, false, ReadBaseStratification.CycleBin.QUINTILE_5, ReadBaseStratification.binnedReadCycleStratifier},

                {0, true, 1, ReadBaseStratification.baseCycleStratifier},
                {1, true, 2, ReadBaseStratification.baseCycleStratifier},
                {2, true, 3, ReadBaseStratification.baseCycleStratifier},
                {3, true, 4, ReadBaseStratification.baseCycleStratifier},
                {4, true, 5, ReadBaseStratification.baseCycleStratifier},
                {5, true, 6, ReadBaseStratification.baseCycleStratifier},
                {6, true, 7, ReadBaseStratification.baseCycleStratifier},
                {7, true, 8, ReadBaseStratification.baseCycleStratifier},

                {7, false, 1, ReadBaseStratification.baseCycleStratifier},
                {6, false, 2, ReadBaseStratification.baseCycleStratifier},
                {5, false, 3, ReadBaseStratification.baseCycleStratifier},
                {4, false, 4, ReadBaseStratification.baseCycleStratifier},
                {3, false, 5, ReadBaseStratification.baseCycleStratifier},
                {2, false, 6, ReadBaseStratification.baseCycleStratifier},
                {1, false, 7, ReadBaseStratification.baseCycleStratifier},
                {0, false, 8, ReadBaseStratification.baseCycleStratifier},

                //in order to curb huge insert lengths that arise from chimeric reads,
                //the insert-length stratifier caps the insert length by 10 * read-length
                {0, true, Math.min(100,10*8), ReadBaseStratification.insertLengthStratifier},

                {0, true, (byte)'A' , ReadBaseStratification.baseQualityStratifier},
                {1, true, (byte)'B' , ReadBaseStratification.baseQualityStratifier},
                {2, true, (byte)'C' , ReadBaseStratification.baseQualityStratifier},
                {3, true, (byte)'D' , ReadBaseStratification.baseQualityStratifier},
                {4, true, (byte)'E' , ReadBaseStratification.baseQualityStratifier},
                {5, true, (byte)'F' , ReadBaseStratification.baseQualityStratifier},
                {6, true, (byte)'G' , ReadBaseStratification.baseQualityStratifier},
                {7, true, (byte)'H' , ReadBaseStratification.baseQualityStratifier},

                {0, true, 40 , ReadBaseStratification.mappingQualityStratifier},
        };
    }

    @Test(dataProvider = "baseStratifierData")
    public void testReadBaseStratifier(final int offset, final boolean readStrandPositive, final Object expectedStratum, final ReadBaseStratification.RecordAndOffsetStratifier<?> recordAndOffsetStratifier) {
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr1", 2_000_000);
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord("rgID");
        samFileHeader.addReadGroup(readGroupRecord);

        final SAMRecord samRecord = new SAMRecord(samFileHeader);

        samRecord.setReadName("62A40AAXX101028:2:93:3981:7576");
        samRecord.setAttribute(SAMTag.RG.name(), readGroupRecord.getId());
        samRecord.setAttribute(SAMTag.NM.name(), 1);

        samRecord.setReadBases("CATGGGGA".getBytes());
        samRecord.setBaseQualities("ABCDEFGH".getBytes());
        samRecord.setInferredInsertSize(100);
        samRecord.setFirstOfPairFlag(true);
        samRecord.setSecondOfPairFlag(false);
        samRecord.setReadUnmappedFlag(false);
        samRecord.setMateUnmappedFlag(false);
        samRecord.setReadPairedFlag(true);
        samRecord.setMateNegativeStrandFlag(readStrandPositive);
        samRecord.setReadNegativeStrandFlag(!readStrandPositive);
        samRecord.setMappingQuality(40);

        samRecord.setAlignmentStart(1);
        SamLocusIterator.RecordAndOffset recordAndOffset = new SamLocusIterator.RecordAndOffset(samRecord, offset);
        SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, 1);

        final boolean offSetOOB = offset < 0 || offset >= samRecord.getReadBases().length;
        final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, offSetOOB ? (byte) 'N' : samRecord.getReadBases()[offset]);

        Assert.assertEquals(recordAndOffsetStratifier.stratify(recordAndOffset, locusAndReference), expectedStratum, recordAndOffsetStratifier.getSuffix());
    }

    @DataProvider
    public Object[][] homopolymerStratifierData() {
        return new Object[][]{
                {-1, true, (byte) 'G', null, ReadBaseStratification.homoPolymerLengthStratifier},
                {0, true, (byte) 'G', 0, ReadBaseStratification.homoPolymerLengthStratifier},
                {1, true, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {2, true, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {3, true, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {4, true, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {5, true, (byte) 'G', 2, ReadBaseStratification.homoPolymerLengthStratifier},
                {6, true, (byte) 'G', 3, ReadBaseStratification.homoPolymerLengthStratifier},
                {7, true, (byte) 'G', 4, ReadBaseStratification.homoPolymerLengthStratifier},
                {8, true, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},

                {-1, false, (byte) 'G', null, ReadBaseStratification.homoPolymerLengthStratifier},
                {0, false, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {1, false, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {2, false, (byte) 'G', 4, ReadBaseStratification.homoPolymerLengthStratifier},
                {3, false, (byte) 'G', 3, ReadBaseStratification.homoPolymerLengthStratifier},
                {4, false, (byte) 'G', 2, ReadBaseStratification.homoPolymerLengthStratifier},
                {5, false, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {6, false, (byte) 'G', 9, ReadBaseStratification.homoPolymerLengthStratifier},
                {7, false, (byte) 'G', 8, ReadBaseStratification.homoPolymerLengthStratifier},
                {8, false, (byte) 'G', 7, ReadBaseStratification.homoPolymerLengthStratifier},
                {9, false, (byte) 'G', 6, ReadBaseStratification.homoPolymerLengthStratifier},
                {10, false, (byte) 'G', 5, ReadBaseStratification.homoPolymerLengthStratifier},
                {11, false, (byte) 'G', 4, ReadBaseStratification.homoPolymerLengthStratifier},
                {12, false, (byte) 'G', 3, ReadBaseStratification.homoPolymerLengthStratifier},
                {13, false, (byte) 'G', 2, ReadBaseStratification.homoPolymerLengthStratifier},
                {14, false, (byte) 'G', 1, ReadBaseStratification.homoPolymerLengthStratifier},
                {15, false, (byte) 'G', 0, ReadBaseStratification.homoPolymerLengthStratifier},
                {16, false, (byte) 'G', null, ReadBaseStratification.homoPolymerLengthStratifier},

                {7, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.SHORT_HOMOPOLYMER, new Pair<>('G', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {8, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.SHORT_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {9, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.SHORT_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {10, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.SHORT_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {11, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.SHORT_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {12, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.SHORT_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {13, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.LONG_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {14, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.LONG_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {15, true, (byte) 'A', new Pair<>(ReadBaseStratification.LongShortHomopolymer.LONG_HOMOPOLYMER, new Pair<>('A', 'A')), ReadBaseStratification.binnedHomopolymerStratifier.get()},
                {16, true, (byte) 'A', null, ReadBaseStratification.binnedHomopolymerStratifier.get()},
        };
    }

    @Test(dataProvider = "homopolymerStratifierData")
    public void testHomopolymerLength(final int offset, final boolean readStrandPositive, final byte referenceBase, final Object expectedStratum, final ReadBaseStratification.RecordAndOffsetStratifier<?> recordAndOffsetStratifier) {
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr1", 200);
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        final SAMRecord samRecord = new SAMRecord(samFileHeader);

        samRecord.setReadBases("CATGGGGAAAAAAAAA".getBytes());
        samRecord.setReadUnmappedFlag(false);
        samRecord.setReadNegativeStrandFlag(!readStrandPositive);
        samRecord.setAlignmentStart(1);
        SamLocusIterator.RecordAndOffset recordAndOffset = new SamLocusIterator.RecordAndOffset(samRecord, offset);
        SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, 1);

        final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, referenceBase);

        Assert.assertEquals(recordAndOffsetStratifier.stratify(recordAndOffset, locusAndReference), expectedStratum);
    }

    @Test()
    public void testHomopolymerLength2() {
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr1", 2_000_000);
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        final SAMRecord samRecord = new SAMRecord(samFileHeader);

        final String refString="GGT";
        final String recString="GGa";
        samRecord.setReadBases(recString.getBytes());
        samRecord.setReadUnmappedFlag(false);
        samRecord.setReadNegativeStrandFlag(false);
        samRecord.setAlignmentStart(1);
        SamLocusIterator.RecordAndOffset recordAndOffset = new SamLocusIterator.RecordAndOffset(samRecord, 2);
        SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, 1);

        final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, refString.getBytes()[2]);

        Assert.assertEquals(ReadBaseStratification.homoPolymerLengthStratifier.stratify(recordAndOffset, locusAndReference), (Integer) 2);
    }
}
