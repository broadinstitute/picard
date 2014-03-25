package net.sf.picard.util;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CollectionUtil;
import org.testng.Assert;
import org.testng.annotations.*;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

/**
 * Tests the IntervalList class
 */
public class IntervalListTest {


    final SAMFileHeader fileHeader;

    final IntervalList list1, list2, list3;

    public IntervalListTest() {
        fileHeader = IntervalList.fromFile(new File("testdata/net/sf/picard/intervallist/IntervalListchr123_empty.interval_list")).getHeader();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        list1 = new IntervalList(fileHeader);
        list2 = new IntervalList(fileHeader);
        list3 = new IntervalList(fileHeader);


        list1.add(new Interval("1", 1, 100));     //de-facto: 1:1-200 1:202-300     2:100-150 2:200-300
        list1.add(new Interval("1", 101, 200));
        list1.add(new Interval("1", 202, 300));
        list1.add(new Interval("2", 200, 300));
        list1.add(new Interval("2", 100, 150));

        list2.add(new Interval("1", 50, 150));   //de-facto 1:50-150 1:301-500      2:1-150 2:250-270 2:290-400
        list2.add(new Interval("1", 301, 500));
        list2.add(new Interval("2", 1, 150));
        list2.add(new Interval("2", 250, 270));
        list2.add(new Interval("2", 290, 400));

        list3.add(new Interval("1", 25, 400));    //de-facto 1:25-400                2:200-600                            3:50-470
        list3.add(new Interval("2", 200, 600));
        list3.add(new Interval("3", 50, 470));
    }


    @DataProvider(name = "intersectData")
    public Object[][] intersectData() {
        final IntervalList intersect123 = new IntervalList(fileHeader);
        final IntervalList intersect12 = new IntervalList(fileHeader);
        final IntervalList intersect13 = new IntervalList(fileHeader);
        final IntervalList intersect23 = new IntervalList(fileHeader);

        intersect123.add(new Interval("1", 50, 150));
        intersect123.add(new Interval("2", 250, 270));
        intersect123.add(new Interval("2", 290, 300));

        intersect12.add(new Interval("1", 50, 150));
        intersect12.add(new Interval("2", 100, 150));
        intersect12.add(new Interval("2", 250, 270));
        intersect12.add(new Interval("2", 290, 300));

        intersect13.add(new Interval("1", 25, 200));
        intersect13.add(new Interval("1", 202, 300));
        intersect13.add(new Interval("2", 200, 300));

        intersect23.add(new Interval("1", 50, 150));
        intersect23.add(new Interval("1", 301, 400));
        intersect23.add(new Interval("2", 250, 270));
        intersect23.add(new Interval("2", 290, 400));


        return new Object[][]{
                new Object[]{Arrays.asList(list1, list2, list3), intersect123},
                new Object[]{Arrays.asList(list1, list2), intersect12},
                new Object[]{Arrays.asList(list2, list1), intersect12},
                new Object[]{Arrays.asList(list2, list3), intersect23},
                new Object[]{Arrays.asList(list3, list2), intersect23},
                new Object[]{Arrays.asList(list1, list3), intersect13},
                new Object[]{Arrays.asList(list3, list1), intersect13}
        };
    }

    @Test(dataProvider = "intersectData")
    public void testIntersectIntervalLists(final List<IntervalList> lists, final IntervalList list) {
        Assert.assertEquals(
                CollectionUtil.makeCollection(IntervalList.intersection(lists).iterator()),
                CollectionUtil.makeCollection(list.iterator()));
    }

    @DataProvider(name = "mergeData")
    public Object[][] mergeData() {
        final IntervalList merge123 = new IntervalList(fileHeader);
        final IntervalList merge12 = new IntervalList(fileHeader);
        final IntervalList merge23 = new IntervalList(fileHeader);
        final IntervalList merge13 = new IntervalList(fileHeader);


        merge123.add(new Interval("1", 1, 100));     //de-facto: 1:1-200 1:202-300     2:100-150 2:200-300
        merge123.add(new Interval("1", 101, 200));
        merge123.add(new Interval("1", 202, 300));
        merge123.add(new Interval("2", 200, 300));
        merge123.add(new Interval("2", 100, 150));

        merge123.add(new Interval("1", 50, 150));   //de-facto 1:50-150 1:301-500      2:1-150 2:250-270 2:290-400
        merge123.add(new Interval("1", 301, 500));
        merge123.add(new Interval("2", 1, 150));
        merge123.add(new Interval("2", 250, 270));
        merge123.add(new Interval("2", 290, 400));

        merge123.add(new Interval("1", 25, 400));    //de-facto 1:25-400                2:200-600                            3:50-470
        merge123.add(new Interval("2", 200, 600));
        merge123.add(new Interval("3", 50, 470));


        merge12.add(new Interval("1", 1, 100));     //de-facto: 1:1-200 1:202-300     2:100-150 2:200-300
        merge12.add(new Interval("1", 101, 200));
        merge12.add(new Interval("1", 202, 300));
        merge12.add(new Interval("2", 200, 300));
        merge12.add(new Interval("2", 100, 150));

        merge12.add(new Interval("1", 50, 150));   //de-facto 1:50-150 1:301-500      2:1-150 2:250-270 2:290-400
        merge12.add(new Interval("1", 301, 500));
        merge12.add(new Interval("2", 1, 150));
        merge12.add(new Interval("2", 250, 270));
        merge12.add(new Interval("2", 290, 400));

        merge23.add(new Interval("1", 50, 150));   //de-facto 1:50-150 1:301-500      2:1-150 2:250-270 2:290-400
        merge23.add(new Interval("1", 301, 500));
        merge23.add(new Interval("2", 1, 150));
        merge23.add(new Interval("2", 250, 270));
        merge23.add(new Interval("2", 290, 400));

        merge23.add(new Interval("1", 25, 400));    //de-facto 1:25-400                2:200-600                            3:50-470
        merge23.add(new Interval("2", 200, 600));
        merge23.add(new Interval("3", 50, 470));


        merge13.add(new Interval("1", 1, 100));     //de-facto: 1:1-200 1:202-300     2:100-150 2:200-300
        merge13.add(new Interval("1", 101, 200));
        merge13.add(new Interval("1", 202, 300));
        merge13.add(new Interval("2", 200, 300));
        merge13.add(new Interval("2", 100, 150));

        merge13.add(new Interval("1", 25, 400));    //de-facto 1:25-400                2:200-600                            3:50-470
        merge13.add(new Interval("2", 200, 600));
        merge13.add(new Interval("3", 50, 470));


        return new Object[][]{
                new Object[]{Arrays.asList(list1, list2, list3), merge123},
                new Object[]{Arrays.asList(list1, list2), merge12},
                new Object[]{Arrays.asList(list2, list3), merge23},
                new Object[]{Arrays.asList(list1, list3), merge13}

        };
    }


    @Test(dataProvider = "mergeData")
    public void testMergeIntervalLists(final List<IntervalList> lists, final IntervalList list) {
        Assert.assertEquals(
                CollectionUtil.makeCollection(IntervalList.concatenate(lists).iterator()),
                CollectionUtil.makeCollection(list.iterator()));
    }


    @DataProvider(name = "unionData")
    public Object[][] unionData() {
        final IntervalList union123 = new IntervalList(fileHeader);
        final IntervalList union12 = new IntervalList(fileHeader);
        final IntervalList union13 = new IntervalList(fileHeader);
        final IntervalList union23 = new IntervalList(fileHeader);

        union123.add(new Interval("1", 1, 500));
        union123.add(new Interval("2", 1, 150));
        union123.add(new Interval("2", 200, 600));
        union123.add(new Interval("3", 50, 470));

        union12.add(new Interval("1", 1, 200));
        union12.add(new Interval("1", 202, 500));
        union12.add(new Interval("2", 1, 150));
        union12.add(new Interval("2", 200, 400));


        union23.add(new Interval("1", 25, 500));
        union23.add(new Interval("2", 1, 150));
        union23.add(new Interval("2", 200, 600));
        union23.add(new Interval("3", 50, 470));

        union13.add(new Interval("1", 1, 400));
        union13.add(new Interval("2", 100, 150));
        union13.add(new Interval("2", 200, 600));
        union13.add(new Interval("3", 50, 470));


        return new Object[][]{
                new Object[]{Arrays.asList(list1, list2, list3), union123},
                new Object[]{Arrays.asList(list1, list2), union12},
                new Object[]{Arrays.asList(list1, list2), union12},
                new Object[]{Arrays.asList(list2, list3), union23},
                new Object[]{Arrays.asList(list2, list3), union23},
                new Object[]{Arrays.asList(list1, list3), union13},
                new Object[]{Arrays.asList(list1, list3), union13}
        };
    }

    @Test(dataProvider = "unionData", enabled = true)
    public void testUnionIntervalLists(final List<IntervalList> lists, final IntervalList list) {
        Assert.assertEquals(
                CollectionUtil.makeCollection(IntervalList.union(lists).iterator()),
                CollectionUtil.makeCollection(list.iterator()));
    }

    @DataProvider(name = "invertData")
    public Object[][] invertData() {
        final IntervalList invert1 = new IntervalList(fileHeader);
        final IntervalList invert2 = new IntervalList(fileHeader);
        final IntervalList invert3 = new IntervalList(fileHeader);

        final IntervalList full = new IntervalList(fileHeader);
        final IntervalList fullChopped = new IntervalList(fileHeader);
        final IntervalList empty = new IntervalList(fileHeader);


        invert1.add(new Interval("1", 201, 201));
        invert1.add(new Interval("1", 301, fileHeader.getSequence("1").getSequenceLength()));
        invert1.add(new Interval("2", 1, 99));
        invert1.add(new Interval("2", 151, 199));
        invert1.add(new Interval("2", 301, fileHeader.getSequence("2").getSequenceLength()));
        invert1.add(new Interval("3", 1, fileHeader.getSequence("3").getSequenceLength()));

        invert2.add(new Interval("1", 1, 49));
        invert2.add(new Interval("1", 151, 300));
        invert2.add(new Interval("1", 501, fileHeader.getSequence("1").getSequenceLength()));
        invert2.add(new Interval("2", 151, 249));
        invert2.add(new Interval("2", 271, 289));
        invert2.add(new Interval("2", 401, fileHeader.getSequence("2").getSequenceLength()));
        invert2.add(new Interval("3", 1, fileHeader.getSequence("3").getSequenceLength()));

        invert3.add(new Interval("1", 1, 24));
        invert3.add(new Interval("1", 401, fileHeader.getSequence("1").getSequenceLength()));
        invert3.add(new Interval("2", 1, 199));
        invert3.add(new Interval("2", 601, fileHeader.getSequence("2").getSequenceLength()));
        invert3.add(new Interval("3", 1, 49));
        invert3.add(new Interval("3", 471, fileHeader.getSequence("3").getSequenceLength()));

        for (final SAMSequenceRecord samSequenceRecord : fileHeader.getSequenceDictionary().getSequences()) {
            full.add(new Interval(samSequenceRecord.getSequenceName(), 1, samSequenceRecord.getSequenceLength()));

            fullChopped.add(new Interval(samSequenceRecord.getSequenceName(), 1, samSequenceRecord.getSequenceLength() / 2));
            fullChopped.add(new Interval(samSequenceRecord.getSequenceName(), samSequenceRecord.getSequenceLength() / 2 + 1, samSequenceRecord.getSequenceLength()));
        }


        return new Object[][]{
                new Object[]{list1, invert1},
                new Object[]{list2, invert2},
                new Object[]{list3, invert3},
                new Object[]{full, empty},
                new Object[]{empty, full},
                new Object[]{fullChopped, empty}
        };
    }


    @Test(dataProvider = "invertData")
    public void testInvertSquared(final IntervalList list, @SuppressWarnings("UnusedParameters") final IntervalList ignored) throws Exception {
        final IntervalList inverseSquared = IntervalList.invert(IntervalList.invert(list));
        final IntervalList originalClone = new IntervalList(list.getHeader());

        for (final Interval interval : list) {
            originalClone.add(interval);
        }

        Assert.assertEquals(
                CollectionUtil.makeCollection(inverseSquared.iterator()),
                CollectionUtil.makeCollection(originalClone.uniqued().iterator()));
    }

    @Test(dataProvider = "invertData")
    public void testInvert(final IntervalList list, final IntervalList inverse) throws Exception {
        Assert.assertEquals(
                CollectionUtil.makeCollection(IntervalList.invert(list).iterator()),
                CollectionUtil.makeCollection(inverse.iterator()));
    }


    @DataProvider(name = "subtractData")
    public Object[][] subtractData() {
        final IntervalList subtract12_from_3 = new IntervalList(fileHeader);
        final IntervalList subtract1_from_2 = new IntervalList(fileHeader);
        final IntervalList subtract2_from_3 = new IntervalList(fileHeader);
        final IntervalList subtract1_from_3 = new IntervalList(fileHeader);
        final IntervalList subtract3_from_1 = new IntervalList(fileHeader);


        subtract12_from_3.add(new Interval("1", 201, 201));
        subtract12_from_3.add(new Interval("2", 401, 600));
        subtract12_from_3.add(new Interval("3", 50, 470));

        subtract1_from_2.add(new Interval("1", 301, 500));
        subtract1_from_2.add(new Interval("2", 1, 99));
        subtract1_from_2.add(new Interval("2", 301, 400));


        subtract2_from_3.add(new Interval("1", 25, 49));
        subtract2_from_3.add(new Interval("1", 151, 300));
        subtract2_from_3.add(new Interval("2", 200, 249));
        subtract2_from_3.add(new Interval("2", 271, 289));
        subtract2_from_3.add(new Interval("2", 401, 600));
        subtract2_from_3.add(new Interval("3", 50, 470));

        subtract1_from_3.add(new Interval("1", 201, 201));
        subtract1_from_3.add(new Interval("1", 301, 400));
        subtract1_from_3.add(new Interval("2", 301, 600));
        subtract1_from_3.add(new Interval("3", 50, 470));

        subtract3_from_1.add(new Interval("1", 1, 49));    //de-facto 1:25-400                2:200-600                            3:50-470
        subtract3_from_1.add(new Interval("2", 100, 150));


        return new Object[][]{
                new Object[]{CollectionUtil.makeList(list3), CollectionUtil.makeList(list1, list2), subtract12_from_3},
                new Object[]{CollectionUtil.makeList(list2), CollectionUtil.makeList(list1), subtract1_from_2},
                new Object[]{CollectionUtil.makeList(list3), CollectionUtil.makeList(list2), subtract2_from_3},
                new Object[]{CollectionUtil.makeList(list3), CollectionUtil.makeList(list1), subtract1_from_3},
        };
    }


    @Test(dataProvider = "subtractData")
    public void testSubtractIntervalLists(final List<IntervalList> fromLists, final List<IntervalList> whatLists, final IntervalList list) {
        Assert.assertEquals(
                CollectionUtil.makeCollection(IntervalList.subtract(fromLists, whatLists).iterator()),
                CollectionUtil.makeCollection(list.iterator()));
    }


    @DataProvider(name = "VCFCompData")
    public Object[][] VCFCompData() {
        return new Object[][]{
                new Object[]{"testdata/net/sf/picard/intervallist/IntervalListFromVCFTest.vcf", "testdata/net/sf/picard/intervallist/IntervalListFromVCFTestComp.interval_list", false},
                new Object[]{"testdata/net/sf/picard/intervallist/IntervalListFromVCFTest.vcf", "testdata/net/sf/picard/intervallist/IntervalListFromVCFTestCompInverse.interval_list", true},
                new Object[]{"testdata/net/sf/picard/intervallist/IntervalListFromVCFTestManual.vcf", "testdata/net/sf/picard/intervallist/IntervalListFromVCFTestManualComp.interval_list", false},
                new Object[]{"testdata/net/sf/picard/intervallist/IntervalListFromVCFTestManual.vcf", "testdata/net/sf/picard/intervallist/IntervalListFromVCFTestCompInverseManual.interval_list", true}
        };
    }


    @Test(dataProvider = "VCFCompData")
    public void testFromVCF(final String vcf, final String compInterval, final boolean invertVCF) {

        final File vcfFile = new File(vcf);
        final File compIntervalFile = new File(compInterval);

        final IntervalList compList = IntervalList.fromFile(compIntervalFile);
        final IntervalList list = invertVCF ? IntervalList.invert(IntervalList.fromVcf(vcfFile)) : IntervalList.fromVcf(vcfFile);

        compList.getHeader().getSequenceDictionary().assertSameDictionary(list.getHeader().getSequenceDictionary());

        final Collection<Interval> intervals = CollectionUtil.makeCollection(list.iterator());
        final Collection<Interval> compIntervals = CollectionUtil.makeCollection(compList.iterator());

        //assert that the intervals correspond
        Assert.assertEquals(intervals, compIntervals);

        final List<String> intervalNames = new LinkedList<String>();
        final List<String> compIntervalNames = new LinkedList<String>();

        for (final Interval interval : intervals) {
            intervalNames.add(interval.getName());
        }
        for (final Interval interval : compIntervals) {
            compIntervalNames.add(interval.getName());
        }
        //assert that the names match
        Assert.assertEquals(intervalNames, compIntervalNames);

    }

}