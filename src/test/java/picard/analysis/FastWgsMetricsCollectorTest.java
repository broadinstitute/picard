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
    private byte[] highQualities = {30, 30, 30 ,30, 30, 30, 30, 30, 30, 30, 50, 50, 50 ,50, 50, 50, 50, 50, 50, 60, 60, 60, 70 ,70, 70, 80, 80, 90, 90, 90};
    private byte[] qualities = {2, 2, 3 ,3, 3, 4, 4, 4, 4, 1};
    private byte[] refBases = {'A', 'C', 'C', 'T', 'A', 'C', 'G', 'T', 'T', 'C', 'A', 'A', 'T', 'A', 'T', 'T', 'C', 'T', 'T', 'C', 'G', 'A', 'G', 'T', 'C', 'D' , 'G', 'T', 'C', 'D','A', 'G', 'T', 'C', 'T', 'T', 'C', 'G', 'A', 'G', 'T', 'C', 'T', 'T', 'C', 'G', 'C', 'T', 'T', 'C', 'G', 'A', 'G', 'T', 'C', 'D' , 'G', 'T', 'C', 'D','A', 'G', 'T', 'C', 'T', 'T', 'C', 'G', 'A', 'G', 'T', 'C', 'T', 'T', 'C', 'G', 'C', 'T', 'T', 'C', 'G', 'A', 'G', 'T', 'C', 'D' , 'G', 'T', 'C', 'D', 'G', 'A', 'G', 'T', 'C', 'D' , 'G', 'T', 'C', 'D'};
    private SAMRecord record;
    private SAMSequenceRecord sequence;
    private SAMRecord secondRecord;
    private SAMRecord thirdRecord;
    private ReferenceSequence ref;

    @BeforeTest
    public void setUp(){
        String referenceString = ">chrM\nACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTCACCTACGTTCAATATTCTTC";
        ref = new ReferenceSequence("chrM", 0, referenceString.getBytes());
        sequence = new SAMSequenceRecord("chrM", 100);
        record = new SAMRecord(new SAMFileHeader());
        record.setReadName("test");
        record.setBaseQualities(qualities);
        record.setReadBases(refBases);
        secondRecord = generateRecord("test1");
        thirdRecord = generateRecord("test2");
    }

    private SAMRecord generateRecord(String name) {
        SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName(name);
        record.setBaseQualities(highQualities);
        record.setReadBases(refBases);
        return record;
    }

    @Test
    public void testAddInfoForQuality() throws Exception {
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 10);
        EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record1);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record2);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        assertEquals(20, collector.basesExcludedByBaseq);
    }

    @Test
    public void testAddInfoForOverlap(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 5);
        AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 10);
        EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 10, 0);
        firstInfo.add(record1);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 5, 5, 0);
        secondInfo.add(record2);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        collector.addInfo(thirdInfo, ref, false);
        assertEquals(5, collector.basesExcludedByOverlap, "Excluded by overlap:");
    }

    @Test
    public void testAddInfoForCapping(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 1, createIntervalList());
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 10);
        EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 10, 0);
        firstInfo.add(record1);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 10, 0);
        firstInfo.add(record2);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        assertEquals(1, collector.basesExcludedByCapping, "Excluded by capping:");
    }

    @Test
    public void testForExcludedForQualityHistogramArray(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        collectWgsMetrics.INCLUDE_BQ_HISTOGRAM = true;
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 30);
        EdgingRecordAndOffset record = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 30, 0);
        firstInfo.add(record);
        firstInfo.add(EdgingRecordAndOffset.createEndRecord(record));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        long[] expectedResult =  new long[127];
        expectedResult[30] = 10; expectedResult[50] = 9; expectedResult[60] = 3;
        expectedResult[70] = 3;  expectedResult[80] = 2;  expectedResult[90] = 3;
        assertEquals(expectedResult, collector.unfilteredBaseQHistogramArray);
    }

    @Test
    public void testForBaseQualityHetSens(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        collectWgsMetrics.INCLUDE_BQ_HISTOGRAM = true;
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 10, createIntervalList());
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 1);
        EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 30, 0);
        firstInfo.add(record1);
        thirdInfo.add(EdgingRecordAndOffset.createEndRecord(record1));
        EdgingRecordAndOffset record2 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
        firstInfo.add(record2);
        secondInfo.add(EdgingRecordAndOffset.createEndRecord(record2));
        collector.addInfo(firstInfo, ref, false);
        collector.addInfo(secondInfo, ref, false);
        collector.addInfo(thirdInfo, ref, false);

        long[] expectedResult =  new long[127];
        expectedResult[30] = 10; expectedResult[50] = 9; expectedResult[60] = 3;
        expectedResult[70] = 3;  expectedResult[80] = 2;  expectedResult[90] = 3;
        expectedResult[1] = 1;  expectedResult[2] = 2;  expectedResult[3] = 3;
        expectedResult[4] = 4;
    }

    @Test
    public void testForHistogramArray(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 10, createIntervalList());
        long[] templateHistogramArray = {0,1,3,0,0,0,0,0,0,0,0};
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 10);
        AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 20);
        AbstractLocusInfo<EdgingRecordAndOffset> fourthInfo = new AbstractLocusInfo<>(sequence, 30);
        EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(record, 0, 10, 0);
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
    public void testForCollectorWithoutData(){
        long[] templateQualHistogram =  new long[127];
        long[] templateHistogramArray = new long[11];
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 10, createIntervalList());
        assertEquals(templateHistogramArray, collector.unfilteredDepthHistogramArray);
        assertEquals(templateQualHistogram, collector.unfilteredBaseQHistogramArray);
        assertEquals(0, collector.basesExcludedByCapping);
        assertEquals(0, collector.basesExcludedByOverlap);
        assertEquals(0, collector.basesExcludedByBaseq);
    }

    @Test
    public void testAddInfoWithoutOverlap(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        AbstractLocusInfo<EdgingRecordAndOffset> firstInfo = new AbstractLocusInfo<>(sequence, 1);
        AbstractLocusInfo<EdgingRecordAndOffset> secondInfo = new AbstractLocusInfo<>(sequence, 5);
        AbstractLocusInfo<EdgingRecordAndOffset> thirdInfo = new AbstractLocusInfo<>(sequence, 6);
        AbstractLocusInfo<EdgingRecordAndOffset> fourthInfo = new AbstractLocusInfo<>(sequence, 10);
        EdgingRecordAndOffset record1 = EdgingRecordAndOffset.createBeginRecord(secondRecord, 0, 5, 0);
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
    public void testForComplicatedCigar(){
        CollectWgsMetrics collectWgsMetrics = new CollectWgsMetrics();
        FastWgsMetricsCollector collector = new FastWgsMetricsCollector(collectWgsMetrics, 100, createIntervalList());
        AbstractLocusIterator sli = createReadEndsIterator(exampleSamTwoReads);
        while(sli.hasNext()) {
            AbstractLocusInfo<EdgingRecordAndOffset> info = sli.next();
            collector.addInfo(info, ref, false);
        }
        assertEquals(11, collector.basesExcludedByOverlap, "Excluded by overlap:");
    }
}
