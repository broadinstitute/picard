package picard.util;

import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.*;

import org.testng.Assert;
import org.testng.collections.Sets;
import picard.cmdline.ClassFinder;

import static java.lang.Math.abs;

/**
 * Small class implementing some utility functions that are useful for test and interfacing with the TestNG framework.
 */
public class TestNGUtil {

    static final double EPSILON = 1e-300; //a constant near the smallest possible positive value representable in double,

     // not actually the smallest possible value on purpose, since that would be indistinguishable from 0 and then useless.
    // This is small enough to be meaningless, but representable.


    public static abstract class SingleTestUnitTest<TestClazz extends TestNGParameterizable> {
        @DataProvider(name = "testcases")
        Object[][] getParams() {
            return TestNGUtil.generateDataProvider(getTestcases());
        }

        @Test(dataProvider = "testcases")
        public void doMetaTest(final TestClazz testcase) {
            // Delegate to another method to avoid overwriting of annotations.
            doTest(testcase);
        }

        abstract public Iterable<TestClazz> getTestcases();

        abstract public void doTest(TestClazz testcase);
    }

    /**
     * Interface for exposing an implementation for converting an object into a Object[] appropriate for positional
     * parameter-passing as performed by TestNG.  Useful for simplifying data sources.
     */
    public static class TestNGParameterizable {
        public Object[] toObjectArray() {
            return new Object[]{this};
        }
    }

    public static Object[][] generateDataProvider(final TestNGParameterizable[] testcases) {
        return generateDataProvider(Arrays.asList(testcases));
    }

    public static Object[][] generateDataProvider(final Iterable<? extends TestNGParameterizable> testcases) {
        final Iterator<? extends TestNGParameterizable> i = testcases.iterator();
        final List<Object[]> parameterList = new LinkedList<Object[]>();
        while (i.hasNext()) {
            parameterList.add(i.next().toObjectArray());
        }
        Object[][] parameterArray = new Object[parameterList.size()][];
        for (int j = 0; j < parameterList.size(); j++) {
            parameterArray[j] = parameterList.get(j);
        }
        return parameterArray;
    }


    //TestNG doesn't have utility functions with these signatures

    /**
     * Small utility function for determining if two doubles are within a _relative_ accuracy of each other.
     *
     * @param lhs first number
     * @param rhs second number
     * @param accuracy maximal allowed relative difference between the two numbers
     * @return true if numbers are within the relative tolerance of each other
     */
    public static boolean compareDoubleWithAccuracy(final double lhs, final double rhs, final double accuracy) {
        if (accuracy <= 0) throw new IllegalArgumentException("Accuracy must be positive.");
        return abs(lhs - rhs) / (abs(lhs) + abs(rhs) + EPSILON) < accuracy;
    }

    public static void assertEqualDoubleArrays(final double[] lhs, final double[] rhs, final double accuracy) {
        Assert.assertNotNull(lhs);
        Assert.assertNotNull(rhs);

        if (accuracy <= 0) throw new IllegalArgumentException("Accuracy must be positive.");
        Assert.assertEquals(lhs.length, rhs.length, "Arrays not same length: " + lhs.length + " vs. " + rhs.length);

        for (int i = 0; i < lhs.length; ++i) {
            Assert.assertTrue(compareDoubleWithAccuracy(lhs[i], rhs[i], accuracy), "Arrays disagree at position " + i + ":  " + lhs[i] + " vs. " + rhs[i] + ". ");
        }
    }

    /** A Method that returns all the Methods that are annotated with @DataProvider
     * in a given package. Should be moved to htsjdk and used from there
     *
     * @param packageName the package under which to look for classes and methods
     * @return an iterator to collection of Object[]'s consisting of {Method method, Class clazz} pair.
     * where method has the @DataProviderAnnotation and is a member of clazz.
     */
    public static Iterator<Object[]> getDataProviders(final String packageName) {
        List<Object[]> data = new ArrayList<>();
        final ClassFinder classFinder = new ClassFinder();
        classFinder.find(packageName, Object.class);

        for (final Class<?> testClass : classFinder.getClasses()) {
            final int modifiers;
            final Method[] declaredMethods;
            final Method[] methods;

            try {
                modifiers = testClass.getModifiers();
                declaredMethods = testClass.getDeclaredMethods();
                methods = testClass.getMethods();
            } catch (final  NoClassDefFoundError e){
                continue;
            }
            if (Modifier.isAbstract(modifiers) || Modifier.isInterface(modifiers))
                continue;
            Set<Method> methodSet = Sets.newHashSet();
            methodSet.addAll(Arrays.asList(declaredMethods));
            methodSet.addAll(Arrays.asList(methods));

            for (final Method method : methodSet) {
                if (method.isAnnotationPresent(DataProvider.class)) {
                    data.add(new Object[]{method, testClass});
                }
            }
        }

        return data.iterator();
    }

    static public void assertGreaterThan(final double lhs, final double rhs) {
        Assert.assertTrue(lhs > rhs, String.format("Expected inequality is not true: %g > %g", lhs, rhs));
    }

    static public void assertLessThan(final double lhs, final double rhs) {
        Assert.assertTrue(lhs < rhs, String.format("Expected inequality is not true: %g < %g", lhs, rhs));
    }
}
