package picard.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
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

    @DataProvider(name="testData")
    Object[][] testData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned.sam"), null},
                new Object[]{new File("testdata/picard/sam/aligned.sam"), new File("testdata/picard/analysis/directed/CollectHsMetrics/chrM.interval_list")},
                new Object[]{new File("testdata/picard/sam/aligned.sam"), new File("testdata/picard/sam/contiguous.interval_list")}
        };
    }

    @DataProvider(name="overWriteTestData")
    Object[][] overwriteTestData() {
        return new Object[][]{
                new Object[]{new File("testdata/picard/sam/aligned.sam"), null, true},
                new Object[]{new File("testdata/picard/sam/aligned.sam"), null, false}
        };
    }

    @Test(dataProvider = "testData")
    public void testTagWriting(final File input, final File interval_list) throws IOException {
        final File clpOutput = File.createTempFile("AddOATag", ".bam");
        clpOutput.deleteOnExit();

        runAddOATag(input, clpOutput, interval_list, null);

        validateBamForOA(clpOutput, interval_list, null);
    }

    @Test(dataProvider = "overWriteTestData")
    public void testOverWrite(final File input, final File interval_list, Boolean overWriteTag) throws IOException {
        final File firstPassOutput = File.createTempFile("FirstPassAddOATag", ".bam");
        firstPassOutput.deleteOnExit();

        final File secondPassOutput = File.createTempFile("SecondPassAddOATag", ".bam");
        secondPassOutput.deleteOnExit();

        runAddOATag(input, firstPassOutput, interval_list, true);

        runAddOATag(firstPassOutput, secondPassOutput, interval_list, overWriteTag);

        validateBamForOA(secondPassOutput, interval_list, overWriteTag);
    }

    private void runAddOATag(final File input, final File output, final File intervalList, Boolean overwriteTag) throws IOException {
        final ArrayList<String> args = new ArrayList<String>(){
            {
                add("INPUT=" + input);
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

    private void validateBamForOA(final File input, final File intervalList, Boolean overwriteTag) {
        final SamReader reader = SamReaderFactory.makeDefault().open(input);
        final SAMRecordIterator iterator = reader.iterator();

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
            Assert.assertTrue(OASplitOnColon.length > 1);
        }
        for (String s : OASplitOnColon) {
            Assert.assertTrue(s.matches("^[\\w]+,[\\w]+,[-+],[\\w]+,[\\w]+,[\\w]*"));
        }

    }


}
