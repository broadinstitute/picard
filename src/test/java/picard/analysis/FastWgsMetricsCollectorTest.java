package picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.AbstractLocusInfo;
import htsjdk.samtools.util.AbstractLocusIterator;
import htsjdk.samtools.util.EdgingRecordAndOffset;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import static org.testng.Assert.assertEquals;
import static picard.analysis.CollectWgsMetricsTestUtils.createIntervalList;
import static picard.analysis.CollectWgsMetricsTestUtils.createReadEndsIterator;
import static picard.analysis.CollectWgsMetricsTestUtils.exampleSamTwoReads;

public class FastWgsMetricsCollectorTest {
    final private byte[] highQualities = {30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 50, 50, 50, 50, 50, 50, 50, 50, 50, 60, 60, 60, 70, 70, 70, 80, 80, 90, 90, 90};
    //unclear why there are so few qualities here....I guess they are not needed, but this is not really valid.
    final private byte[] qualities = {2, 2, 3, 3, 3, 4, 4, 4, 4, 1};
    final private byte[] readBases = "ACCTACGTTCAATATTCTTCGAGTCDGTCDAGTCTTCGAGTCTTCGCTTCGAGTCDGTCDAGTCTTCGAGTCTTCGCTTCGAGTCDGTCDGAGTCDGTCD".getBytes();
    private SAMRecord record;
    private SAMSequenceRecord sequence;
    private SAMRecord secondRecord;
    private SAMRecord thirdRecord;
    private ReferenceSequence ref;

    @BeforeTest
    public void setUp() {
        final String referenceString = ">chrM\nACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTC";
        ref = new ReferenceSequence("chrM", 0, referenceString.getBytes());
        sequence = new SAMSequenceRecord("chrM", 100);
        record = new SAMRecord(new SAMFileHeader());
        record.setReadName("test");
        record.setBaseQualities(qualities);
        record.setReadBases(readBases);
        secondRecord = generateRecord("test1");
        thirdRecord = generateRecord("test2");
    }

    private SAMRecord generateRecord(String name) {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(name);
        record.setBaseQualities(highQualities);
        record.setReadBases(readBases);
        return record;
    }

    @Test
    public void testAddInfoForQuality() throws Exception {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 10);
        final EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record1);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        final EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record2);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        assertEquals(20, collector.basesExcludedByBaseq);
    }

    @Test
    public void testAddInfoForOverlap() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 5);
        final AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 10);
        final EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 10, 0);
        firstInfo.add(record1);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        final EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 5, 5, 0);
        secondInfo.add(record2);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        collector.addInfo(thirdInfo, ref, false);
        assertEquals(5, collector.basesExcludedByOverlap, "Excluded by overlap:");
    }

    @Test
    public void testAddInfoForCapping() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 1, createIntervalList());
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 10);
        final EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 10, 0);
        firstInfo.add(record1);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        final EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 10, 0);
        firstInfo.add(record2);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        assertEquals(1, collector.basesExcludedByCapping, "Excluded by capping:");
    }

    @Test
    public void testForExcludedForQualityHistogramArray() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        collectWgsMetrics.INCLUDE_BQ_HISTOGRAM = true;
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 30);
        final EdgingRecordAndOffset record = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 30, 0);
        firstInfo.add(record);
        firstInfo.add(EdgingRecordAndOffset.createEndRecord(record));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        final long[] expectedResult = new long[127];
        expectedResult[30] = 10;
        expectedResult[50] = 9;
        expectedResult[60] = 3;
        expectedResult[70] = 3;
        expectedResult[80] = 2;
        expectedResult[90] = 3;
        assertEquals(expectedResult, collector.unfilteredBaseQHistogramArray);
    }

    @Test
    public void testForBaseQualityHetSens() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        collectWgsMetrics.INCLUDE_BQ_HISTOGRAM = true;
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 10, createIntervalList());
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 1);
        final EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 30, 0);
        firstInfo.add(record1);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        final EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record2);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        collector.addInfo(thirdInfo, ref, false);

        final long[] expectedResult = new long[127];
        expectedResult[30] = 10;
        expectedResult[50] = 9;
        expectedResult[60] = 3;
        expectedResult[70] = 3;
        expectedResult[80] = 2;
        expectedResult[90] = 3;
        expectedResult[1] = 1;
        expectedResult[2] = 2;
        expectedResult[3] = 3;
        expectedResult[4] = 4;
    }

    @Test
    public void testForHistogramArray() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 10, createIntervalList());
        final long[] templateHistogramArray = {0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0};
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 10);
        final AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 20);
        final AbstractLocusInfo<EdgingRecordAndOffset> fourthInfo = new AbstractLocusInfo<>(sequence, 30);
        final EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record1);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 30, 0);
        firstInfo.add(record2);
        fourthInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        EdgingRecordAndOffset record3 = EdgingRecordAndOffset.createBeginRecord(thirdRecord, 0, 20, 0);
        firstInfo.add(record3);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record3));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        collector.addInfo(thirdInfo, ref, false);
        collector.addInfo(fourthInfo, ref, false);
        assertEquals(templateHistogramArray, collector.unfilteredDepthHistogramArray);
    }

    @Test
    public void testForCollectorWithoutData() {
        final long[] templateQualHistogram = new long[127];
        final long[] templateHistogramArray = new long[11];
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 10, createIntervalList());
        assertEquals(templateHistogramArray, collector.unfilteredDepthHistogramArray);
        assertEquals(templateQualHistogram, collector.unfilteredBaseQHistogramArray);
        assertEquals(0, collector.basesExcludedByCapping);
        assertEquals(0, collector.basesExcludedByOverlap);
        assertEquals(0, collector.basesExcludedByBaseq);
    }

    @Test
    public void testAddInfoWithoutOverlap() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        final AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        final AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 5);
        final AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 6);
        final AbstractLocusInfo<EdgingRecordAndOffset> fourthInfo = new AbstractLocusInfo<>(sequence, 10);
        final EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 5, 0);
        firstInfo.add(record1);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 6, 5, 6);
        thirdInfo.add(record2);
        fourthInfo.add(EdgingRecordAndOffset.createEndRecord(record2));

        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        assertEquals(0, collector.basesExcludedByOverlap, "Excluded by overlap:");
    }

    @Test
    public void testForComplicatedCigar() {
        final CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        final FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        final AbstractLocusIterator sli = createReadEndsIterator(exampleSamTwoReads);
        while (sli.hasNext()) {
            AbstractLocusInfo<EdgingRecordAndOffset> info = sli.next();
            collector.addInfo(info, ref, false);
        }
        assertEquals(collector.basesExcludedByOverlap, 13,"Excluded by overlap:");
    }
}
