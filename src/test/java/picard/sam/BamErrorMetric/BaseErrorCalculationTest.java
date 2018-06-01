package picard.sam.BamErrorMetric;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SamLocusIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;

/**
 * A slew of unit tests for the various Calculators
 */

public class BaseErrorCalculationTest {

    @Test
    public void testSimpleErrorCalculator() {

        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr1", 200);
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        samFileHeader.addSequence(samSequenceRecord);
        final SAMRecord samRecord = new SAMRecord(samFileHeader);

        samRecord.setReadBases("CgTGtGGAcAAAgAAA".getBytes());
        final byte[] refBases = "CATGGGGAAAAAAAAA".getBytes();
        final int n = refBases.length;

        samRecord.setReadUnmappedFlag(false);
        samRecord.setReadNegativeStrandFlag(true);
        samRecord.setAlignmentStart(1);
        samRecord.setReferenceIndex(0);

        final BaseErrorCalculation.SimpleErrorCalculator baseErrorCalculator = new BaseErrorCalculation.SimpleErrorCalculator();

        for (int i = 0; i < n; i++) {

            SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, i + 1);
            final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, refBases[i]);

            SamLocusIterator.RecordAndOffset recordAndOffset = new SamLocusIterator.RecordAndOffset(samRecord, i);
            baseErrorCalculator.addBase(recordAndOffset, locusAndReference);
        }
        final ErrorMetrics.SimpleErrorMetric metric = baseErrorCalculator.getMetric();
        metric.calculateDerivedFields();
        Assert.assertEquals(metric.TOTAL_BASES, n);
        Assert.assertEquals(metric.ERROR_BASES, 4L);
    }

    @Test
    public void testOverlappingErrorCalculator() {

        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord("chr1", 200);
        final SAMFileHeader samFileHeader = new SAMFileHeader();
        samFileHeader.addSequence(samSequenceRecord);

        final SAMRecordSetBuilder builder = new SAMRecordSetBuilder();

        final List<SAMRecord> samRecords = builder.addPair("Read1234", 0, 1, 1,
                false, false, "16M", "16M", true, false, 20);
        final SAMRecord samRecord1 = samRecords.get(0);
        final SAMRecord samRecord2 = samRecords.get(1);

        samRecord1.setReadBases("CgTGtGGAcAAAgAAA".getBytes());
        samRecord2.setReadBases("CcTGGtGAcAAAgAAA".getBytes());
        final byte[] refBases = "CATGGGGAAAAAAAAA".getBytes();

        final BaseErrorCalculation.OverlappingReadsErrorCalculator overlappingErrorCalculator1 = new BaseErrorCalculation.OverlappingReadsErrorCalculator();
        final BaseErrorCalculation.OverlappingReadsErrorCalculator overlappingErrorCalculator2 = new BaseErrorCalculation.OverlappingReadsErrorCalculator();
        final int length = refBases.length;

        for (int i = 0; i < length; i++) {

            SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, i + 1);
            SamLocusIterator.RecordAndOffset recordAndOffset1 = new SamLocusIterator.RecordAndOffset(samRecord1, i);
            SamLocusIterator.RecordAndOffset recordAndOffset2 = new SamLocusIterator.RecordAndOffset(samRecord2, i);
            locusInfo.add(recordAndOffset1);
            locusInfo.add(recordAndOffset2);

            final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, refBases[i]);

            overlappingErrorCalculator1.addBase(recordAndOffset1, locusAndReference);
            overlappingErrorCalculator2.addBase(recordAndOffset2, locusAndReference);
        }

        final ErrorMetrics.OverlappingErrorMetric metric1 = overlappingErrorCalculator1.getMetric();
        metric1.calculateDerivedFields();
        Assert.assertEquals(metric1.TOTAL_BASES, length);
        Assert.assertEquals(metric1.NUM_BASES_WITH_OVERLAPPING_READS, length);
        Assert.assertEquals(metric1.NUM_DISAGREES_WITH_REFERENCE_ONLY, 2L);
        Assert.assertEquals(metric1.NUM_DISAGREES_WITH_REF_AND_MATE, 1L);
        Assert.assertEquals(metric1.NUM_THREE_WAYS_DISAGREEMENT, 1L);

        final ErrorMetrics.OverlappingErrorMetric metric2 = overlappingErrorCalculator1.getMetric();
        metric2.calculateDerivedFields();
        Assert.assertEquals(metric2.TOTAL_BASES, length);
        Assert.assertEquals(metric2.NUM_BASES_WITH_OVERLAPPING_READS,  length);
        Assert.assertEquals(metric2.NUM_DISAGREES_WITH_REFERENCE_ONLY, 2L);
        Assert.assertEquals(metric2.NUM_DISAGREES_WITH_REF_AND_MATE, 1L);
        Assert.assertEquals(metric2.NUM_THREE_WAYS_DISAGREEMENT, 1L);
    }
}
