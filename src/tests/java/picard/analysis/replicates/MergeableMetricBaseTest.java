package picard.analysis.replicates;

import org.testng.Assert;
import org.testng.annotations.Test;

public class MergeableMetricBaseTest {

    class TestMergeableMetric extends MergeableMetricBase {
        @MergeByAdding
        Integer boxedInt = 1;
        @MergeByAdding
        int unboxedInt = 2;

        @MergeByAdding
        Double boxedDouble = 3D;
        @MergeByAdding
        double unboxedDouble = 4D;

        @MergeByAdding
        Long boxedLong = 5L;
        @MergeByAdding
        long unboxedLong = 6L;

        @MergeByAdding
        Float boxedFloat = 7F;
        @MergeByAdding
        float unboxedFloat = 8F;

        @MergeByAdding
        Short boxedShort = 9;
        @MergeByAdding
        short unboxedShort = 10;

        @MergeByAdding
        Byte boxedByte = 11;
        @MergeByAdding
        byte unboxedByte = 12;

        @MergeByAssertEquals
        String mustBeEqualString = "hello";

        @MergeByAssertEquals
        Double mustBeEqualDouble = 0.5;

        @MergeByAssertEquals
        boolean mustBeEqualUnboxedBoolean = false;

        @NoMergingIsDerived
        double ratioIntValues;

        @Override
        public void calculateDerivedFields() {
            ratioIntValues = boxedInt / (double) unboxedInt;
        }
    }

    @Test
    public void testMerging() {
        final TestMergeableMetric metric1 = new TestMergeableMetric(), metric2 = new TestMergeableMetric();
        metric1.merge(metric2);

        Assert.assertEquals(metric1.boxedInt, (Integer) 2);
        Assert.assertEquals(metric1.unboxedInt, 4);

        Assert.assertEquals(metric1.boxedDouble, 6D);
        Assert.assertEquals(metric1.unboxedDouble, 8D);

        Assert.assertEquals(metric1.boxedLong, (Long) 10L);
        Assert.assertEquals(metric1.unboxedLong, 12L);

        Assert.assertEquals(metric1.boxedFloat, 14F);
        Assert.assertEquals(metric1.unboxedFloat, 16F);

        Assert.assertEquals(metric1.boxedShort, (Short) (short) 18);
        Assert.assertEquals(metric1.unboxedShort, (short) 20);

        Assert.assertEquals(metric1.boxedByte, (Byte) (byte) 22);
        Assert.assertEquals(metric1.unboxedByte, 24);

        Assert.assertEquals(metric1.mustBeEqualDouble, metric2.mustBeEqualDouble);
        Assert.assertEquals(metric1.mustBeEqualString, metric2.mustBeEqualString);
        Assert.assertEquals(metric1.mustBeEqualUnboxedBoolean, metric2.mustBeEqualUnboxedBoolean);

        metric1.calculateDerivedFields();

        Assert.assertEquals(metric1.ratioIntValues, 0.5D);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMergingUnequalString() {

        final TestMergeableMetric metric1 = new TestMergeableMetric(), metric2 = new TestMergeableMetric();
        metric1.mustBeEqualString = "goodbye";

        Assert.assertFalse(metric1.canMerge(metric2));
        metric1.merge(metric2);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMergingUnequalDouble() {

        final TestMergeableMetric metric1 = new TestMergeableMetric(), metric2 = new TestMergeableMetric();
        metric1.mustBeEqualDouble = 1D;

        Assert.assertFalse(metric1.canMerge(metric2));
        metric1.merge(metric2);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMergingUnequalBoolean() {

        final TestMergeableMetric metric1 = new TestMergeableMetric(), metric2 = new TestMergeableMetric();
        metric1.mustBeEqualUnboxedBoolean = true;

        Assert.assertFalse(metric1.canMerge(metric2));
        metric1.merge(metric2);
    }

    private class TestMergeableMericIllegal extends MergeableMetricBase {
        Integer undecorated = 0;
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testIllegalClass() {
        final TestMergeableMericIllegal illegal1 = new TestMergeableMericIllegal(), illegal2 = new TestMergeableMericIllegal();

        illegal1.merge(illegal2);
    }

    private class TestDerivedMergableMetric extends TestMergeableMetric {
        @MergeByAdding
        Integer anotherBoxed = 1;
    }

    @Test
    public void TestMergingDerivedClass() {
        final TestMergeableMetric instance1 = new TestMergeableMetric();
        final TestDerivedMergableMetric instance2 = new TestDerivedMergableMetric();

        instance1.merge(instance2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void TestMergingSuperClass() {
        final TestMergeableMetric instance1 = new TestMergeableMetric();
        final TestDerivedMergableMetric instance2 = new TestDerivedMergableMetric();

        instance2.merge(instance1);
    }

    @Test
    public void TestCanMerge() {
        final TestMergeableMetric instance1 = new TestMergeableMetric();
        instance1.unboxedInt=1;
        final TestDerivedMergableMetric instance2 = new TestDerivedMergableMetric();
        instance2.unboxedInt=2;

        instance1.merge(instance2);
        Assert.assertEquals(instance1.unboxedInt, 3);
    }
}
