package picard.sam;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import picard.sam.StabilizeQualityScores.*;

public class StabilizeQualityScoresTest extends CommandLineProgramTest {
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
        ThresholdBinner tbin = new ThresholdBinner(20);

        int[] inquals = {0,1,2,15,19,20,21,50};
        Assert.assertEquals(box(tbin.bin(inquals)), Arrays.asList(2,2,2,2,2,20,20,20));
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

    @Test
    public void testMergeStabilizer() {

    }
}
