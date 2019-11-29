package picard.illumina;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.BarcodeEditDistanceQuery;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import static org.testng.Assert.assertEquals;

public class DistanceMetricTest {

    @DataProvider
    public Iterator<Object[]> distanceTestData() {
        final byte[][] barcode1      = new byte[][]{"ACGT".getBytes(), "AACCGT".getBytes()};
        final byte[][] read1         = new byte[][]{"ACGT".getBytes(), "AACCGT".getBytes()};
        final byte[][] readOneError  = new byte[][]{"AgGT".getBytes(), "AACCGT".getBytes()};
        final byte[][] readTwoError  = new byte[][]{"ACGT".getBytes(), "AACCtT".getBytes()};
        final byte[][] readTwoIns    = new byte[][]{"ACGT".getBytes(), "AACaCG".getBytes()};
        final byte[][] readOneDel    = new byte[][]{"ACTG".getBytes(), "AACCGT".getBytes()};
        final byte[][] readOneNoCall    = new byte[][]{"NN.n".getBytes(), "AACCGT".getBytes()};

        final byte[][] readOneNoisy  = new byte[][]{"AAAA".getBytes(), "AACCGT".getBytes()};

        final byte[][] goodQualities = new byte[][]{"IIII"     .getBytes(), "IIIIIII".getBytes()};
        final byte[][] badQuality    = new byte[][]{"I\u0001II".getBytes(), "IIIIIII".getBytes()};

        final List<Object[]> tests = new ArrayList<>();

        BarcodeEditDistanceQuery query = new BarcodeEditDistanceQuery(barcode1, read1, goodQualities, 2, 4);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 0});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 0});
        tests.add(new Object[]{query, DistanceMetric.FREE, 0});

        query = new BarcodeEditDistanceQuery(barcode1, readOneError, goodQualities, 2, 4);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 1});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 1});
        tests.add(new Object[]{query, DistanceMetric.FREE, 1});

        query = new BarcodeEditDistanceQuery(barcode1, readTwoError, goodQualities, 2, 4);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 1});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 1});
        tests.add(new Object[]{query, DistanceMetric.FREE, 1});

        query = new BarcodeEditDistanceQuery(barcode1, readTwoIns, goodQualities, 2, 4);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 3});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 3});
        tests.add(new Object[]{query, DistanceMetric.FREE, 1});

        query = new BarcodeEditDistanceQuery(barcode1, readOneDel, goodQualities, 2, 4);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 2});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 2});
        tests.add(new Object[]{query, DistanceMetric.FREE, 1});

        query = new BarcodeEditDistanceQuery(barcode1, readOneNoisy, goodQualities, 2, 4);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 3});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 3});
        tests.add(new Object[]{query, DistanceMetric.FREE, 3});

        query = new BarcodeEditDistanceQuery(barcode1, readOneNoisy, goodQualities, 2, 1);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 2});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 2});
        tests.add(new Object[]{query, DistanceMetric.FREE, 2});

        query = new BarcodeEditDistanceQuery(barcode1, read1, badQuality, 2, 3);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 0});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 1});
        tests.add(new Object[]{query, DistanceMetric.FREE, 0});

        query = new BarcodeEditDistanceQuery(barcode1, readOneError, badQuality, 2, 3);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 0});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 1});
        tests.add(new Object[]{query, DistanceMetric.FREE, 0});

        query = new BarcodeEditDistanceQuery(barcode1, readOneNoCall, goodQualities, 2, 3);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 0});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 0});
        tests.add(new Object[]{query, DistanceMetric.FREE, 0});

        query = new BarcodeEditDistanceQuery(barcode1, readOneNoisy, goodQualities, 2, 3);
        tests.add(new Object[]{query, DistanceMetric.LENIENT_HAMMING, 3});
        tests.add(new Object[]{query, DistanceMetric.HAMMING, 3});
        tests.add(new Object[]{query, DistanceMetric.FREE, 3});

        return tests.iterator();
    }

    @Test(dataProvider = "distanceTestData")
    public void testDistances(final BarcodeEditDistanceQuery query, final DistanceMetric metric, final int expectedDistance) {
        assertEquals(metric.distance(query), expectedDistance);
    }
}