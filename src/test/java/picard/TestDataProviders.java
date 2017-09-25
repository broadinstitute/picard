package picard;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.ClassFinder;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;

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
    public void testAllDataProviders() throws IllegalAccessException, InstantiationException, InvocationTargetException {
        int i = 0;
        final ClassFinder classFinder = new ClassFinder();
        classFinder.find("picard", Object.class);
        for (Class testClass : classFinder.getClasses()) {
            if (Modifier.isAbstract(testClass.getModifiers())) continue;
            for (final Method method : testClass.getMethods()) {
                if (method.isAnnotationPresent(DataProvider.class)) {
                    System.err.println("Method: " + testClass.getName() + ":" + method.getName());
                    method.invoke(testClass.newInstance());
                    i++;
                }
            }
        }
        Assert.assertNotSame(i,0);
        System.err.println("Found: "+ i + " @DataProviders.");
    }
}
