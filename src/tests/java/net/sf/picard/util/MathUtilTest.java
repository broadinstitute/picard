package net.sf.picard.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * @author mccowan
 */
public class MathUtilTest {
    @Test
    public void logMathTest() {
        Assert.assertEquals(MathUtil.LOG_10_MATH.getLogValue(10), 1d, 0.00001d);
        Assert.assertEquals(MathUtil.LOG_10_MATH.getNonLogValue(1), 10d, 0.00001d);
        Assert.assertEquals(MathUtil.LOG_10_MATH.mean(5, 5, 5), 5d, 0.00001d);
        // http://www.wolframalpha.com/input/?i=log10%2810%5E1.23%2B10%5E4.56%2B10%5E99999%29
        Assert.assertEquals(MathUtil.LOG_10_MATH.sum(1.23, 4.56, 2), 4.5613970317323586660874152202433434022756298235604568d, 0.00001d);
        // http://www.wolframalpha.com/input/?i=log10%2810%5E1.23+*+10%5E4.56+*+10%5E2%29
        Assert.assertEquals(MathUtil.LOG_10_MATH.product(1.23, 4.56, 2), 7.7899999999999999999999999999999999999999999999999999d, 0.00001d);
    }

    @DataProvider
    public Object[][] seqMethodTestCases() {
        return new Object[][] {
                new Object[] {0d, 5d, 1d,     new double[] {0,1,2,3,4,5}},
                new Object[] {0d, 0.5d, 0.1d, new double[] {0, 0.1, 0.2, 0.3, 0.4, 0.5}},
                new Object[] {0d, 0.5d, 0.11d,new double[] {0, 0.11, 0.22, 0.33, 0.44}},
                new Object[] {50d, 55d, 1.25d,new double[] {50, 51.25, 52.5, 53.75, 55}},
                new Object[] {10d, 0d, 02d,    new double[] {}},
                new Object[] {10d, 0d, -2d,    new double[] {10, 8, 6, 4, 2, 0}},
        };
    }

    @Test(dataProvider="seqMethodTestCases")
    public void testSeqGeneration(final double from, final double to, final double by, final double[] expected) {
        final double[] actual = MathUtil.seq(from, to, by);
        Assert.assertEquals(actual.length, expected.length);

        for (int i=0; i<expected.length; ++i) {
            Assert.assertTrue(Math.abs(actual[i] - expected[i]) < 0.0000001);
        }
    }
}
