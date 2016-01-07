package picard.fingerprint;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

import static org.testng.Assert.*;

/**
 * Created by farjoun on 8/27/15.
 */
public class FingerprintCheckerTest {

    @Test
    public void testRandomSublist() throws Exception {

        List<Integer> list = new ArrayList<>();
        list.add(1);
        list.add(2);
        list.add(3);

        Assert.assertEquals(list, FingerprintChecker.randomSublist(list, 3));
        Assert.assertEquals(list, FingerprintChecker.randomSublist(list, 4));

        Assert.assertEquals(FingerprintChecker.randomSublist(list, 2).size(), 2);
    }
}