package net.sf.picard.util;

import static java.lang.Math.abs;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 *
 */
public class HistogramTest {

    @Test(dataProvider = "histogramData")
    public void testHistogramFunctions(final int[] values, final double mean, final double stdev) {
        final Histogram<Integer> histo = new Histogram<Integer>();
        for (int value : values) {
            histo.increment(value);
        }

        final double m = histo.getMean();
        final double sd = histo.getStandardDeviation();

        Assert.assertEquals(round(mean), round(m), "Means are not equal");
        Assert.assertEquals(round(stdev), round(sd), "Stdevs are not equal");
    }

    @DataProvider(name = "histogramData")
    public Object[][] histogramData() {
        return new Object[][] {
            new Object[] {new int[] {1,2,3,4,5,6,7,8,9,10} , 5.5d, 3.027650d },
            new Object[] {new int[] {1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9}, 6.333333d, 2.236068d },
            new Object[] {new int[] {-5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15}, 5d, 6.204837d }
        };
    }

    private double round(final double in) {
        long l = (long) (in * 10000);
        return l / 10000d;
    }

}
