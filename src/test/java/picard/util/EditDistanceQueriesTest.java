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
public class EditDistanceQueriesTest {

    private Object[][] levensteinDistanceData() {
        return new Object[][]{
                {"one", "ont", 3, 1},
                {"o-e", "oe-", 3, 1},
                {"on-e", "one-", 3, 1},
                {"o-ne", "one-", 3, 1},
                {"o-nce", "once-", 3, 1},
                {"on-ce", "once-", 3, 1},
                {"on-ce", "onc-e", 3, 2},
                {"n-ce", "nc-e", 3, 2},
                {"", "", 2, 0},
                {"ABCA-CABC", "ABCACABC-", 2, 1}, // gaps in the end should be ignored
                {"BC AB- ABC", "-BC AB ABC", 2, 2}, // but not in the beginning
                {"HELLO!", "HELLO!", 2, 0},
                {"very large distances are not measured well",
                        "so this should return the threshold plus 1", 5, 6},
                {"-ELLO", "ELLO-", 3, 1}
        };
    }

    @DataProvider
    private Object[][] levensteinDistanceWithNoCallsData() {
        return new Object[][]{
                {"one", ".ne", 0},
                {"any", "...", 0},
                {"hello", "hell.", 0},
                {"TCTGCAAGGCCAGAAG", "T.TGCAAGGCCAGAAG", 0},
                {       "CGCCTTCC",
                        "C.CCcTCC",    1},
                {       "AATCGCTG",
                        "A.cCGCTG", 1},
                {"goodbye--", "good.bye-", 1},
                {"goodddbye", "good.bye-", 1},
                {"--goodbye", "good.bye-", 3},
        };
    }

    @DataProvider
    public Object[][] levensteinDistanceSymmetricData() {
        final List<Object[]> tests = Arrays.asList(levensteinDistanceData());

        final List<Object[]> symmetricTests = new ArrayList<>(2 * tests.size());
        symmetricTests.addAll(tests);

        tests.forEach(t -> {
            final Object[] newTest = Arrays.copyOf(t, t.length);
            newTest[0] = t[1];
            newTest[1] = t[0];
            symmetricTests.add(newTest);
        });
        return symmetricTests.toArray(new Object[symmetricTests.size()][]);
    }

    @Test(dataProvider = "levensteinDistanceSymmetricData")
    public void levenshteinDistanceTest(final String string1, final String string2, final int threshold, final int expectedDistance) {
        distanceHelper(string1, string2, threshold, expectedDistance);
    }

    @Test(dataProvider = "levensteinDistanceWithNoCallsData")
    public void levenshteinDistanceWithNocallsTest(final String string1, final String string2, final int expectedDistance) {
        distanceHelper(string1, string2, 5, expectedDistance);
    }

    private void distanceHelper(final String string1, final String string2, final int threshold, final int expectedDistance) {
        final byte[] string1Bytes = string1.getBytes();
        final byte[] string2Bytes = string2.getBytes();

        Assert.assertEquals(new SingleBarcodeDistanceMetric(string1Bytes, string2Bytes, null, 0, threshold).freeDistance(), expectedDistance);

        // make sure that the code doesn't modify the input strings...
        Assert.assertEquals(string1Bytes, string1.getBytes());
        Assert.assertEquals(string2Bytes, string2.getBytes());
    }
}