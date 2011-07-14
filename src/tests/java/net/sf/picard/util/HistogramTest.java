package net.sf.picard.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import static java.lang.Math.abs;

/**
 *
 */
public class HistogramTest {

    @Test(dataProvider = "histogramData")
    public void testHistogramFunctions(final int[] values, final double mean, final double stdev, final Integer trimByWidth) {
        final Histogram<Integer> histo = new Histogram<Integer>();
        for (int value : values) {
            histo.increment(value);
        }

        if (trimByWidth != null) histo.trimByWidth(trimByWidth);
        final double m = histo.getMean();
        final double sd = histo.getStandardDeviation();

        Assert.assertEquals(round(mean), round(m), "Means are not equal");
        Assert.assertEquals(round(stdev), round(sd), "Stdevs are not equal");
    }

    @DataProvider(name = "histogramData")
    public Object[][] histogramData() {
        return new Object[][] {
            new Object[] {new int[] {1,2,3,4,5,6,7,8,9,10} , 5.5d, 3.027650d, null },
            new Object[] {new int[] {1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9}, 6.333333d, 2.236068d, null  },
            new Object[] {new int[] {-5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15}, 5d, 6.204837d, null  },
                new Object[] {new int[] {1,2,3,4,5,6,7,8,9,10, 11, 11, 12, 100, 1000} , 5.5d, 3.027650d, 10 },
                new Object[] {new int[] {1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9, 20, 20, 21, 25, 25}, 6.333333d, 2.236068d, 11  },
                new Object[] {new int[] {-5, -4, -3, -2, -1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 101, 102, 103, 200, 2000}, 5d, 6.204837d, 20  }
        };
    }

    @Test
    public void testGeometricMean() {
        final int[] is = new int[] {4,4,4,4,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8};
        final Histogram<Integer> histo = new Histogram<Integer>();
        for (final int i : is) histo.increment(i);
        Assert.assertTrue(abs(histo.getGeometricMean() - 6.216797) < 0.00001);
    }

    @Test(dataProvider = "medianTestData")
    public void testMedian(final int [] values, final double median) {
        final Histogram<Integer> histo = new Histogram<Integer>();
        for (final int i : values) histo.increment(i);
        Assert.assertEquals(histo.getMedian(), median);
    }

    @DataProvider(name = "medianTestData")
    public Object[][] medianTestData() {
        return new Object[][] {
                new Object[] {new int[] {} , 0d},
                new Object[] {new int[] {999} , 999d},
                new Object[] {new int[] {1,2,3,4,5,6} , 3.5d},
                new Object[] {new int[] {5,5,5,5,5,6,6} , 5d},
                new Object[] {new int[] {5,5,5,5,5,6,6,6,6,6} , 5.5d},
        };
    }

    @Test
    public void testMad() {
        final int[] is = new int[] {4,4,4,4,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8};
        final Histogram<Integer> histo = new Histogram<Integer>();
        for (final int i : is) histo.increment(i);

        Assert.assertEquals(7d, histo.getMedian());
        Assert.assertEquals(1d, histo.getMedianAbsoluteDeviation());
        Assert.assertTrue(abs(histo.estimateSdViaMad() - 1.4826) < 0.0001);
    }


    private double round(final double in) {
        long l = (long) (in * 10000);
        return l / 10000d;
    }

}
