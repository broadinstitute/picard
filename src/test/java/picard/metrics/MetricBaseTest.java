package picard.metrics;

import htsjdk.samtools.metrics.MetricBase;
import org.broadinstitute.barclay.argparser.ClassFinder;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Field;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class MetricBaseTest {

    @DataProvider(name = "testMetricClasses")
    public Iterator<Object[]> getMetricClasses() {
        final ClassFinder classFinder = new ClassFinder();

        classFinder.find("picard", MetricBase.class);
        return classFinder.getClasses().stream()
                .map(c -> new Object[]{c})
                .iterator();
    }

    @Test(dataProvider = "testMetricClasses")
    public void testUniqueFields(final Class metricClass) {
        final Set<String> metricFields = new HashSet<>();

        for (final Field f : metricClass.getFields()) {
            Assert.assertFalse(!metricFields.add(f.getName()), metricClass + " has a naming collision in field " + f.getName());
        }
    }
}
