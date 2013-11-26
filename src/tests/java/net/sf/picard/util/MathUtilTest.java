package net.sf.picard.util;

import org.testng.Assert;
import org.testng.annotations.Test;

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
}
