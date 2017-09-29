package picard;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.NoInjection;
import org.testng.annotations.Test;
import picard.cmdline.ClassFinder;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * This test is a mechanism to check that none of the data-providers fail to run.
 * It is needed because in the case that a data-provider fails (for some reason, perhaps a change in some other code
 * causes it to throw an exception) the tests that rely on it will be sliently skipped.
 * The only mention of this will be in test logs but since we normally avoid reading these logs and rely on the
 * exit code, it will look like all the tests have passed.
 *
 * @author Yossi Farjoun
 */
public class TestDataProviders {

    @Test
    public void IndependentTestOfDataProviderTest() throws IllegalAccessException, InvocationTargetException, InstantiationException {
        testAllDataProvidersdata();
    }

    @DataProvider(name = "DataprovidersThatDontTestThemselves")
    public Iterator<Object[]> testAllDataProvidersdata() throws IllegalAccessException, InstantiationException, InvocationTargetException {

        List<Object[]> data = new ArrayList<>();
        final ClassFinder classFinder = new ClassFinder();
        classFinder.find("picard", Object.class);

        for (final Class<?> testClass : classFinder.getClasses()) {
            if (Modifier.isAbstract(testClass.getModifiers())) continue;
            for (final Method method : testClass.getMethods()) {
                if (method.isAnnotationPresent(DataProvider.class)) {
                    data.add(new Object[]{method, testClass});
                }
            }
        }
        Assert.assertTrue(data.size() > 1);
        Assert.assertEquals(data.stream().filter(c -> ((Method) c[0]).getName().equals("testAllDataProvidersdata")).count(), 1);

        return data.iterator();
    }

    // @NoInjection annotations required according to this test:
    // https://github.com/cbeust/testng/blob/master/src/test/java/test/inject/NoInjectionTest.java
    @Test(dataProvider = "DataprovidersThatDontTestThemselves")
    public void testDataProviderswithDP(@NoInjection final Method method, final Class clazz) throws IllegalAccessException, InstantiationException, InvocationTargetException {
        System.err.println("Method: " + method + " Class: " + clazz.getName());

        Object instance = clazz.newInstance();

        // Some tests assume that the @BeforeSuite methods will be called before the @DataProviders
        for (final Method otherMethod : clazz.getMethods()) {
            if (otherMethod.isAnnotationPresent(BeforeSuite.class)) {
                otherMethod.invoke(instance);
            }
        }

        method.invoke(instance);
    }
}
