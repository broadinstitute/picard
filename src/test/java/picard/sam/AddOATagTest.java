package picard.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.OverlapDetector;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class AddOATagTest {
    @DataProvider(name="testOAData")
    Object[][] testOAData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned_withoutOA.sam"), new File("testdata/picard/sam/aligned_withOA.sam")},
        };
    }

    @DataProvider(name="testIntervalListData")
    Object[][] testIntervalListData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned.sam"), null},
                new Object[]{new File("testdata/picard/sam/aligned.sam"), new File("testdata/picard/analysis/directed/CollectHsMetrics/chrM.interval_list")},
                new Object[]{new File("testdata/picard/sam/aligned.sam"), new File("testdata/picard/sam/contiguous.interval_list")}
        };
    }

    @DataProvider(name="testOverwriteData")
    Object[][] overwriteTestData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned.sam"), null, true},
                new Object[]{new File("testdata/picard/sam/aligned.sam"), null, false}
        };
    }

    @Test(dataProvider = "testOAData")
    public void testWritingOATag(final File testSam, final File truthSam) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(testSam, clpOutput, null, null);

        validateOATAG(clpOutput, truthSam);
    }

    @Test(dataProvider = "testIntervalListData")
    public void testIntervalList(final File inputSam, final File interval_list) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(inputSam, clpOutput, interval_list, null);

        validateBamHasOA(clpOutput, interval_list, null);
    }

    @Test(dataProvider = "testOverwriteData")
    public void testOverWrite(final File inputSam, final File interval_list, Boolean overWriteTag) throws IOException {
        final File firstPassOutput = File.createTempFile("FirstPassAddOATag", ".bam");
        firstPassOutput.deleteOnExit();
        final File secondPassOutput = File.createTempFile("SecondPassAddOATag", ".bam");
        secondPassOutput.deleteOnExit();

        // make first pass bam that only has one value per OA tag
        runAddOATag(inputSam, firstPassOutput, interval_list, true);

        // make bam we want to test the overwrite option with
        runAddOATag(firstPassOutput, secondPassOutput, interval_list, overWriteTag);

        validateBamHasOA(secondPassOutput, interval_list, overWriteTag);
    }

    private void runAddOATag(final File inputSam, final File output, final File intervalList, Boolean overwriteTag) throws IOException {
        final ArrayList<String> args = new ArrayList<String>(){
            {
                add("INPUT=" + inputSam);
                add("OUTPUT=" + output);
            }
        };
        if (intervalList != null) {
            args.add("INTERVAL_LIST=" + intervalList);
        }
        if (overwriteTag != null) {
            args.add("OVERWRITE_TAG=" + overwriteTag);
        }
        AddOATag addOATag = new AddOATag();
        Assert.assertEquals(addOATag.instanceMain(args.toArray(new String[args.size()])), 0, "Running addOATag did not succeed");
    }

    private void validateOATAG(final File testSam, final File truthSam) {
        final ArrayList<String> truthOAValues = new ArrayList<>();
        final ArrayList<String> testOAValues = new ArrayList<>();

        SAMRecordIterator iterator = SamReaderFactory.makeDefault().open(truthSam).iterator();
        while (iterator.hasNext()){
            SAMRecord rec = iterator.next();
            truthOAValues.add(rec.getStringAttribute("OA"));
        }

        iterator = SamReaderFactory.makeDefault().open(testSam).iterator();
        while (iterator.hasNext()){
            SAMRecord rec = iterator.next();
            testOAValues.add(rec.getStringAttribute("OA"));
        }

        Assert.assertTrue(truthOAValues.equals(testOAValues));
    }

    private void validateBamHasOA(final File input, final File intervalList, Boolean overwriteTag) {
        final SAMRecordIterator iterator = SamReaderFactory.makeDefault().open(input).iterator();

        OverlapDetector<Interval> detector = AddOATag.getOverlapDetectorFromIntervalListFile(intervalList, 0, 0);
        while (iterator.hasNext()){
            SAMRecord rec = iterator.next();
            if (!rec.getReadUnmappedFlag() && (detector == null || detector.overlapsAny(rec))){
                Assert.assertNotNull(rec.getAttribute("OA"));
                if (overwriteTag != null){
                    validateOverwriteTag(rec, overwriteTag);
                }
            } else {
                Assert.assertNull(rec.getAttribute("OA"));
            }
        }
    }

    private void validateOverwriteTag(SAMRecord rec, Boolean overwriteTag) {
        String OAValue = rec.getStringAttribute("OA");
        String[] OASplitOnColon = OAValue.split(";");
        if (overwriteTag) {
            Assert.assertTrue(OASplitOnColon.length == 1);
        } else {
            Assert.assertTrue(OASplitOnColon.length == 2);
        }
        for (String s : OASplitOnColon) {
            Assert.assertTrue(s.matches("^[\\w]+,[\\w]+,[-+],[\\w]+,[\\w]+,[\\w]*"));
        }
    }
}
