package picard;

import org.testng.Assert;
import org.testng.annotations.*;
import org.testng.annotations.Test;
import picard.util.TestNGUtil;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.*;
import java.util.stream.StreamSupport;

/**
 * This test is a mechanism to check that none of the data-providers fail to run.
 * It is needed because in the case that a data-provider fails (for some reason, perhaps a change in some other code
 * causes it to throw an exception) the tests that rely on it will be silently skipped.
 * The only mention of this will be in test logs but since we normally avoid reading these logs and rely on the
 * exit code, it will look like all the tests have passed.
 *
 * @author Yossi Farjoun
 */
public class TestDataProviders {

    @Test
    public void independentTestOfDataProviderTest() throws Exception {
        final Iterator<Object[]> data = TestNGUtil.getDataProviders("picard");

        Assert.assertTrue(data.hasNext(), "Found no data from testAllDataProvidersdata. Something is wrong.");

        Assert.assertEquals( StreamSupport.stream(Spliterators.spliteratorUnknownSize(data, 0), false)
                        .filter(c -> ((Method) c[0]).getName().equals("testAllDataProvidersData")).count(), 1,
                "getDataProviders didn't find testAllDataProvidersData, which is in this class. Something is wrong.");
    }

    @DataProvider(name = "DataprovidersThatDontTestThemselves")
    public Iterator<Object[]> testAllDataProvidersData() throws Exception {
        return TestNGUtil.getDataProviders("picard");
    }

    // runs all the @DataProviders it gets from DataProvidersThatDontTestThemselves.
    // runs the BeforeSuite and BeforeClass methods of the same class before it runs the provider itself.

    // @NoInjection annotations required according to this test:
    // https://github.com/cbeust/testng/blob/master/src/test/java/test/inject/NoInjectionTest.java
    @Test(dataProvider = "DataprovidersThatDontTestThemselves")
    public void testDataProviderswithDP(@NoInjection final Method method, final Class clazz) throws
            IllegalAccessException, InstantiationException {

        Object instance = clazz.newInstance();

        Set<Method> methodSet = new HashSet<>();
        methodSet.addAll(Arrays.asList(clazz.getDeclaredMethods()));
        methodSet.addAll(Arrays.asList(clazz.getMethods()));

        // Some tests assume that the @BeforeSuite methods will be called before the @DataProviders
        // However, it might be problematic to invoke @BeforeSuite twice, thus if a class has @BeforeSuite
        // method, we will skip testing for that @DataProvider.
        if (methodSet.stream().anyMatch(m->m.isAnnotationPresent(BeforeSuite.class))) return;

        for (final Method otherMethod : methodSet) {
            if (otherMethod.isAnnotationPresent(BeforeClass.class)) {
                try {
                    otherMethod.setAccessible(true);
                    otherMethod.invoke(instance);
                } catch (IllegalAccessException | InvocationTargetException e) {
                    throw new IllegalStateException(String.format("@BeforeClass threw an exception (%s::%s). Dependent tests will be skipped. Please fix.", clazz.getName(), method.getName()), e);
                }
            }
        }

        try {
            method.setAccessible(true);
            method.invoke(instance);
        } catch (IllegalAccessException | InvocationTargetException e) {
            throw new IllegalStateException(String.format("@DataProvider threw an exception (%s::%s). Dependent tests will be skipped. Please fix.", clazz.getName(), method.getName()), e);
        }
    }
}

