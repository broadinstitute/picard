package picard.metrics;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.analysis.CollectRawWgsMetrics;
import picard.analysis.CollectWgsMetricsWithNonZeroCoverage;
import picard.analysis.directed.HsMetrics;
import picard.analysis.directed.TargetMetrics;
import picard.analysis.directed.TargetedPcrMetrics;

import java.lang.reflect.Field;
import java.util.HashSet;
import java.util.Set;


public class MetricBaseTest {

    @DataProvider(name= "testMetricClasses")
    public Object[][] getMetricClasses(){
        return new Object[][]{
                {TargetedPcrMetrics.class},
                {HsMetrics.class},
                {TargetMetrics.class},
                {CollectRawWgsMetrics.RawWgsMetrics.class},
                {CollectWgsMetricsWithNonZeroCoverage.WgsMetricsWithNonZeroCoverage.class}
            };
        }

    @Test(dataProvider = "testMetricClasses")
    public void testUniqueFields(Class metricClass) {

        Set<String> metricFields = new HashSet<>();

        for (final Field f : metricClass.getFields()) {
            Assert.assertTrue(metricFields.add(f.getName()), f.getName() + " has a naming collision");
        }
    }
}

