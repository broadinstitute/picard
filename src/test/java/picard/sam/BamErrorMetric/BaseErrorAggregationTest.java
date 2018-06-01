package picard.sam.BamErrorMetric;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.SamLocusIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Created by farjoun on 5/31/18.
 */
public class BaseErrorAggregationTest {

    @Test
    public void testBaseErrorAggregation() {
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

        BaseErrorAggregation<?> baseErrorAggregation = new BaseErrorAggregation<>(BaseErrorCalculation.OverlappingReadsErrorCalculator::new, ReadBaseStratification.readOrdinalityStratifier);
        final int length = refBases.length;

        for (int i = 0; i < length; i++) {

            SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, i + 1);
            SamLocusIterator.RecordAndOffset recordAndOffset1 = new SamLocusIterator.RecordAndOffset(samRecord1, i);
            SamLocusIterator.RecordAndOffset recordAndOffset2 = new SamLocusIterator.RecordAndOffset(samRecord2, i);
            locusInfo.add(recordAndOffset1);
            locusInfo.add(recordAndOffset2);

            final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, refBases[i]);

            baseErrorAggregation.addBase(recordAndOffset1, locusAndReference);
            baseErrorAggregation.addBase(recordAndOffset2, locusAndReference);
        }

        final ErrorMetrics.ErrorMetricBase[] metrics = baseErrorAggregation.getMetrics();
        final ErrorMetrics.OverlappingErrorMetric metric1 = (ErrorMetrics.OverlappingErrorMetric) metrics[0];
        metric1.calculateDerivedFields();
        Assert.assertEquals(metric1.COVARIATE, ReadBaseStratification.ReadOrdinality.FIRST.name());
        Assert.assertEquals(metric1.TOTAL_BASES, length);
        Assert.assertEquals(metric1.NUM_BASES_WITH_OVERLAPPING_READS, length);
        Assert.assertEquals(metric1.NUM_DISAGREES_WITH_REFERENCE_ONLY, 2L);
        Assert.assertEquals(metric1.NUM_DISAGREES_WITH_REF_AND_MATE, 1L);
        Assert.assertEquals(metric1.NUM_THREE_WAYS_DISAGREEMENT, 1L);

        final ErrorMetrics.OverlappingErrorMetric metric2 = (ErrorMetrics.OverlappingErrorMetric) metrics[1];
        metric2.calculateDerivedFields();
        Assert.assertEquals(metric2.COVARIATE, ReadBaseStratification.ReadOrdinality.SECOND.name());
        Assert.assertEquals(metric2.TOTAL_BASES, length);
        Assert.assertEquals(metric2.NUM_BASES_WITH_OVERLAPPING_READS, length);
        Assert.assertEquals(metric2.NUM_DISAGREES_WITH_REFERENCE_ONLY, 2L);
        Assert.assertEquals(metric2.NUM_DISAGREES_WITH_REF_AND_MATE, 1L);
        Assert.assertEquals(metric2.NUM_THREE_WAYS_DISAGREEMENT, 1L);
    }

    @Test
    public void testBaseErrorAggregation2() {
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

        BaseErrorAggregation<?> baseErrorAggregation = new BaseErrorAggregation<>(BaseErrorCalculation.OverlappingReadsErrorCalculator::new, ReadBaseStratification.readDirectionStratifier);
        final int length = refBases.length;

        for (int i = 0; i < length; i++) {

            SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, i + 1);
            SamLocusIterator.RecordAndOffset recordAndOffset1 = new SamLocusIterator.RecordAndOffset(samRecord1, i);
            SamLocusIterator.RecordAndOffset recordAndOffset2 = new SamLocusIterator.RecordAndOffset(samRecord2, i);
            locusInfo.add(recordAndOffset1);
            locusInfo.add(recordAndOffset2);

            final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, refBases[i]);

            baseErrorAggregation.addBase(recordAndOffset1, locusAndReference);
            baseErrorAggregation.addBase(recordAndOffset2, locusAndReference);
        }

        final ErrorMetrics.ErrorMetricBase[] metrics = baseErrorAggregation.getMetrics();
        final ErrorMetrics.OverlappingErrorMetric metric1 = (ErrorMetrics.OverlappingErrorMetric) metrics[0];
        metric1.calculateDerivedFields();
        Assert.assertEquals(metric1.COVARIATE, ReadBaseStratification.ReadDirection.POSITIVE.toString());
        Assert.assertEquals(metric1.TOTAL_BASES, length);
        Assert.assertEquals(metric1.NUM_BASES_WITH_OVERLAPPING_READS, length);
        Assert.assertEquals(metric1.NUM_DISAGREES_WITH_REFERENCE_ONLY, 2L);
        Assert.assertEquals(metric1.NUM_DISAGREES_WITH_REF_AND_MATE, 1L);
        Assert.assertEquals(metric1.NUM_THREE_WAYS_DISAGREEMENT, 1L);

        final ErrorMetrics.OverlappingErrorMetric metric2 = (ErrorMetrics.OverlappingErrorMetric) metrics[1];
        metric2.calculateDerivedFields();
        Assert.assertEquals(metric2.COVARIATE, ReadBaseStratification.ReadDirection.NEGATIVE.toString());
        Assert.assertEquals(metric2.TOTAL_BASES, length);
        Assert.assertEquals(metric2.NUM_BASES_WITH_OVERLAPPING_READS, length);
        Assert.assertEquals(metric2.NUM_DISAGREES_WITH_REFERENCE_ONLY, 2L);
        Assert.assertEquals(metric2.NUM_DISAGREES_WITH_REF_AND_MATE, 1L);
        Assert.assertEquals(metric2.NUM_THREE_WAYS_DISAGREEMENT, 1L);
    }

    @Test
    public void testBaseErrorAggregation3() {
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

        BaseErrorAggregation<?> baseErrorAggregation = new BaseErrorAggregation<>(BaseErrorCalculation.SimpleErrorCalculator::new, ReadBaseStratification.referenceBaseStratifier);
        final int length = refBases.length;

        for (int i = 0; i < length; i++) {

            SamLocusIterator.LocusInfo locusInfo = new SamLocusIterator.LocusInfo(samSequenceRecord, i + 1);
            SamLocusIterator.RecordAndOffset recordAndOffset1 = new SamLocusIterator.RecordAndOffset(samRecord1, i);
            SamLocusIterator.RecordAndOffset recordAndOffset2 = new SamLocusIterator.RecordAndOffset(samRecord2, i);
            locusInfo.add(recordAndOffset1);
            locusInfo.add(recordAndOffset2);

            final SAMLocusAndReferenceIterator.SAMLocusAndReference locusAndReference = new SAMLocusAndReferenceIterator.SAMLocusAndReference(locusInfo, refBases[i]);

            baseErrorAggregation.addBase(recordAndOffset1, locusAndReference);
            baseErrorAggregation.addBase(recordAndOffset2, locusAndReference);
        }

        final ErrorMetrics.ErrorMetricBase[] metrics = baseErrorAggregation.getMetrics();
        final ErrorMetrics.SimpleErrorMetric metricA = Arrays.stream(metrics).map(a -> (ErrorMetrics.SimpleErrorMetric) a)
                .filter(a -> a.COVARIATE.equals("A")).findFirst().orElse(null);
        metricA.calculateDerivedFields();
        Assert.assertEquals(metricA.COVARIATE, "A");
        Assert.assertEquals(metricA.TOTAL_BASES, 11);
        Assert.assertEquals(metricA.ERROR_BASES, 3);

        final ErrorMetrics.SimpleErrorMetric metricC = Arrays.stream(metrics).map(a -> (ErrorMetrics.SimpleErrorMetric) a)
                .filter(a -> a.COVARIATE.equals("C")).findFirst().orElse(null);
        metricA.calculateDerivedFields();
        Assert.assertEquals(metricC.COVARIATE, "C");
        Assert.assertEquals(metricC.TOTAL_BASES, 5);
        Assert.assertEquals(metricC.ERROR_BASES, 1);

        final ErrorMetrics.SimpleErrorMetric metricG = Arrays.stream(metrics).map(a -> (ErrorMetrics.SimpleErrorMetric) a)
                .filter(a -> a.COVARIATE.equals("G")).findFirst().orElse(null);
        metricA.calculateDerivedFields();
        Assert.assertEquals(metricG.COVARIATE, "G");
        Assert.assertEquals(metricG.TOTAL_BASES, 5);
        Assert.assertEquals(metricG.ERROR_BASES, 1);

        final ErrorMetrics.SimpleErrorMetric metricT = Arrays.stream(metrics).map(a -> (ErrorMetrics.SimpleErrorMetric) a)
                .filter(a -> a.COVARIATE.equals("T")).findFirst().orElse(null);
        metricA.calculateDerivedFields();
        Assert.assertEquals(metricT.COVARIATE, "T");
        Assert.assertEquals(metricT.TOTAL_BASES, 11);
        Assert.assertEquals(metricT.ERROR_BASES, 3);

    }
}