package picard.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by farjoun on 2/15/18.
 */
public class StringDistanceUtilsTest {

    private Object[][] levensteinDistanceData() {
        return new Object[][]{
                {"one", "ont", 1},
                {"o-e", "oe-", 1},
                {"on-e", "one-", 1},
                {"o-ne", "one-", 1},
                {"o-nce", "once-", 1},
                {"on-ce", "once-", 1},
                {"on-ce", "onc-e", 2},
                {"n-ce", "nc-e", 2},
                {"", "", 0},
                {"ABCA-CABC", "ABCACABC-", 1}, // gaps in the end should be ignored
                {"BC AB- ABC", "-BC AB ABC", 2}, // but not in the beginning
                {"HELLO!", "HELLO!", 0},
                {"very large distances are not measured well",
                        "so this should return the threshold plus 1", 6},
                {"-ELLO", "ELLO-", 1},

        };
    }
    @DataProvider
    private Object[][] levensteinDistanceWithNoCallsData() {
        return new Object[][] {
                //see that no-calls are not "charged" for snps:
                {"one", ".ne", 0},
                {"any", "...", 0},
                {"hello", "hell.", 0},
                // but they still count for indels
                {"goodbye--", "good.bye-", 1},
                {"goodddbye", "good.bye-", 1},
                {"--goodbye", "good.bye-", 3},
                {"CGCCTTCC","C.CCCTCC",1},
                {"AATCGCTG","A.CCGCTG",1},
                {"TCTGCAAGGCCAGAAG","T.TGCAAGGCCAGAAG",0}
        };
    }
    @DataProvider
    public Object[][] levensteinDistanceSymmetricData() {
        final List<Object[]> tests = Arrays.asList(levensteinDistanceData());

        final List<Object[]> symmetricTests = new ArrayList<>(2 * tests.size());
        symmetricTests.addAll(tests);
        tests.forEach(t -> symmetricTests.add(new Object[]{t[1], t[0], t[2]}));
        return symmetricTests.toArray(new Object[symmetricTests.size()][]);
    }

    @Test(dataProvider = "levensteinDistanceSymmetricData")
    public void levenshteinDistanceTest(final String string1, final String string2, final int expectedDistance) {
        distanceHelper(string1, string2, expectedDistance);
    }

    @Test(dataProvider = "levensteinDistanceWithNoCallsData")
    public void levenshteinDistanceWithNocallsTest(final String string1, final String string2, final int expectedDistance) {
        distanceHelper(string1, string2, expectedDistance);
    }

    private void distanceHelper(String string1, String string2, int expectedDistance) {
        final byte[] string1Bytes = string1.getBytes();
        final byte[] string2Bytes = string2.getBytes();

        Assert.assertEquals(StringDistanceUtils.levenshteinDistance(string1Bytes, string2Bytes, 5), expectedDistance);

        //make sure that the code doesn't modify the input strings...
        TestNGUtil.assertEquals(string1Bytes, string1.getBytes());
        TestNGUtil.assertEquals(string2Bytes, string2.getBytes());
    }
}