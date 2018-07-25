package picard.sam;

import org.apache.commons.math3.util.Pair;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import picard.sam.StabilizeQualityScores.*;

public class StabilizeQualityScoresTest extends CommandLineProgramTest {

    @DataProvider(name = "StabilizerTestProvider")
    public Object[][] StabilizerTestProvider() {
        final List<Integer> testOQuals = Arrays.asList(0, 2, 3, 19, 20, 5, 201, 2); //silly outlier drags up the average
        final List<RLEElem> testRLE = Arrays.asList(
                new RLEElem(2, 3),
                new RLEElem(20, 2),
                new RLEElem(2, 1),
                new RLEElem(20, 1),
                new RLEElem(2, 1));
        final RLEElem testPrevRLE = new RLEElem(20, 2);

        return new Object[][] {{testOQuals, testPrevRLE, testRLE}};
    }

    @Override
    public String getCommandLineProgramName() { return StabilizeQualityScores.class.getSimpleName(); }

    List<Integer> box(int[] in) {
        return Arrays.stream(in).boxed().collect(Collectors.toList());
    }

    //Arrays.stream() has no override for byte arrays what the hell
    List<Byte> box(byte[] in) {
        List<Byte> outList = new ArrayList<>();
        for (int i = 0; i < in.length; ++i) {
            outList.add(in[i]);
        }
        return outList;
    }

    @Test
    public void testThresholdBinner() {
        ThresholdBinner tbin = new ThresholdBinner(20, null);

        int[] inquals = {0,1,2,15,19,20,21,50};
        Assert.assertEquals(box(tbin.bin(inquals)), Arrays.asList(2,2,2,2,2,20,20,20));
    }

    @Test
    public void testThresholdUpBinner() {
        ThresholdBinner tbin = new ThresholdBinner(20, 30);

        int[] inquals = {0,1,2,15,19,20,21,50};
        Assert.assertEquals(box(tbin.bin(inquals)), Arrays.asList(2,2,2,2,2,30,30,30));
    }

    @Test
    public void testNearestBinner() {
        NearestBinner nbin = new NearestBinner(Arrays.asList(0, 10, 20, 30));

        int[] inquals = {0, 4, 5, 9, 10, 11, 14, 15, 20, 31, 900};
        Assert.assertEquals(box(nbin.bin(inquals)), Arrays.asList(0, 0, 0, 10, 10, 10, 10, 10, 20, 30, 30));
    }

    @Test
    public void testQualsToRLE() {
        int[] inquals = {2,4,4,4,3,3,3,3};
        List<RLEElem> expectedRLE = Arrays.asList(new RLEElem(2, 1), new RLEElem(4, 3), new RLEElem(3, 4));
        Assert.assertEquals(StabilizeQualityScores.toRLE(inquals), expectedRLE);
    }

    @Test
    public void testRLEtoQuals() {
        List<RLEElem> inRLE = Arrays.asList(new RLEElem(2, 1), new RLEElem(4, 3), new RLEElem(3, 4));
        byte[] expectedQuals = {2,4,4,4,3,3,3,3};

        Assert.assertEquals(box(StabilizeQualityScores.toQuals(inRLE)), box(expectedQuals));
    }

    @Test(dataProvider = "StabilizerTestProvider")
    public void testMergeStabilizer(final List<Integer> testOQuals, final RLEElem testPrevRLE, final List<RLEElem> testRLE) {
        Stabilizer merge = new MergeStabilizer(4);

        //should merge correctly with previous RLE
        List<RLEElem> mergePrevious = merge.stabilize(testOQuals, testPrevRLE, testRLE);
        Assert.assertEquals(mergePrevious, Arrays.asList(new RLEElem(testPrevRLE.qual, testPrevRLE.count + testOQuals.size())));

        //should merge correctly with no RLE
        List<RLEElem> mergeBeginning = merge.stabilize(testOQuals, null, testRLE);
        Assert.assertEquals(mergeBeginning, Arrays.asList(new RLEElem(testRLE.get(0).qual, testOQuals.size())));
    }

    @Test(dataProvider = "StabilizerTestProvider")
    public void testDropStabilizer(final List<Integer> testOQuals, final RLEElem testPrevRLE, final List<RLEElem> testRLE) {
        Stabilizer drop = new DropStabilizer(4);

        //should drop correctly with previous RLE
        List<RLEElem> dropPrevious = drop.stabilize(testOQuals, testPrevRLE, testRLE);
        Assert.assertEquals(dropPrevious, Arrays.asList(testPrevRLE, new RLEElem(2, testOQuals.size())));

        //should drop correctly with no RLE
        List<RLEElem> dropBeginning = drop.stabilize(testOQuals, null, testRLE);
        Assert.assertEquals(dropBeginning, Arrays.asList(new RLEElem(2, testOQuals.size())));
    }

    @Test(dataProvider = "StabilizerTestProvider")
    public void testAverageStabilizer(final List<Integer> testOQuals, final RLEElem testPrevRLE, final List<RLEElem> testRLE) {
        ThresholdBinner tbin = new ThresholdBinner(20, null);
        Stabilizer average = new AverageStabilizer(4, tbin);

        //Note that the original quals are set up in such a way that the average of the binned quals is 0,
        //but the binned average of the original quals is 20. This is deliberate :)

        //should average correctly with previous RLE
        List<RLEElem> averagePrevious = average.stabilize(testOQuals, testPrevRLE, testRLE);
        Assert.assertEquals(averagePrevious, Arrays.asList(testPrevRLE, new RLEElem(20, testOQuals.size())));

        //should average correctly with no RLE
        List<RLEElem> averageBeginning = average.stabilize(testOQuals, null, testRLE);
        Assert.assertEquals(averageBeginning, Arrays.asList(new RLEElem(20, testOQuals.size())));
    }

    @Test
    public void testGroupify() {
        final List<RLEElem> testRLE = Arrays.asList(
                new RLEElem(2, 3),
                new RLEElem(20, 4), //group of 2
                new RLEElem(2, 10), //group of 1
                new RLEElem(20, 1),
                new RLEElem(2, 1)); //group of 2

        List<Pair<Boolean, List<RLEElem>>> groups = StabilizeQualityScores.groupify(testRLE, 4);

        Assert.assertEquals(groups,
                Arrays.asList(
                        new Pair<>(true, Arrays.asList(
                                new RLEElem(2, 3),
                                new RLEElem(20, 4))),
                        new Pair<>(false, Arrays.asList(
                                new RLEElem(2, 10))),
                        new Pair<>(true, Arrays.asList(
                                new RLEElem(20, 1),
                                new RLEElem(2, 1)))
                        ));
    }

    //TODO: test actual file writing
}
