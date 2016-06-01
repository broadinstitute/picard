package picard.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static picard.util.MathUtil.divide;

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

        TestNGUtil.assertEqualDoubleArrays(MathUtil.LOG_10_MATH.getLogValue(new double[]{0.1, 1.0, 10.0}), new double[]{-1.0, 0.0, 1.0}, 0.00001d);
        TestNGUtil.assertEqualDoubleArrays(MathUtil.LOG_2_MATH.getLogValue(new double[]{0.5, 1.0, 2.0}), new double[]{-1.0, 0.0, 1.0}, 0.00001d);
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

    @DataProvider
    public Object[][] divideMethodTestCases() {
        return new Object[][] {
                new Object[] {new double [] {1, 2, 3, 4}, new double [] {2, 3, 4, -5}, new double[] {.5, 2.0/3, 3.0/4, -4.0/5}},
                new Object[] {new double [] {100}, new double [] {200}, new double[] {.5}},
                new Object[] {new double [] {0, 4, -3, 2}, new double [] {200, 30, 32, 12}, new double[] {0, 4.0/30, -3.0/32, 2.0/12}},
                new Object[] {new double [] {}, new double [] {}, new double[] {}},

        };
    }

    @Test(dataProvider = "divideMethodTestCases")
    public void testDivide(final double [] numerators, final double [] denominators, final double [] expected) {
        assertEquals(divide(numerators, denominators), expected);
    }

    @DataProvider
    public Object[][] divideMethodFailTestCases() {
        return new Object[][] {
                new Object[] {new double [] { 1, 2, 3, 4}, new double [] {2, 3, 4}},
                new Object[] {new double [] {100}, new double [] {}},
        };
    }

    @Test(dataProvider = "divideMethodFailTestCases",expectedExceptions = RuntimeException.class)
    public void testDivideFail(final double [] lhs, final double [] rhs){
        divide(lhs, rhs);
    }

    //TestNG doesn't have a utility function with this signature....
    private void assertEquals(final double [] actual, final double [] expected) {
        Assert.assertEquals(actual.length,expected.length,"Arrays do not have equal lengths");

        for(int i=0;i<actual.length;++i){
            Assert.assertEquals(actual[i], expected[i],"Array differ at position " +i);
        }
    }
}
