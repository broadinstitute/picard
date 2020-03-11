/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package picard.analysis;

import org.testng.Assert;
import org.testng.annotations.Test;

public class MergeableMetricBaseTest {

    static class TestMergeableMetric extends MergeableMetricBase {
        @MergeByAdding
        public Integer boxedInt = 1;
        @MergeByAdding
        public int unboxedInt = 2;

        @MergeByAdding
        public Double boxedDouble = 3D;
        @MergeByAdding
        public double unboxedDouble = 4D;

        @MergeByAdding
        public Long boxedLong = 5L;
        @MergeByAdding
        public long unboxedLong = 6L;

        @MergeByAdding
        public Float boxedFloat = 7F;
        @MergeByAdding
        public float unboxedFloat = 8F;

        @MergeByAdding
        public Short boxedShort = 9;
        @MergeByAdding
        public short unboxedShort = 10;

        @MergeByAdding
        public Byte boxedByte = 11;
        @MergeByAdding
        public byte unboxedByte = 12;

        @MergeByAssertEquals
        public String mustBeEqualString = "hello";

        @MergeByAssertEquals
        public Double mustBeEqualDouble = 0.5;

        @MergeByAssertEquals
        public boolean mustBeEqualUnboxedBoolean = false;

        @MergeByAdding
        protected Integer protectedIntBase = 1;

        @MergeByAdding
        private Integer privateIntBase = 1;

        @NoMergingIsDerived
        protected double ratioIntValues;

        @Override
        public void calculateDerivedFields() {
            ratioIntValues = boxedInt / (double) unboxedInt;
        }
    }

     static class TestMergeableMetricWithRestrictedMembers extends TestMergeableMetric {

        @MergeByAdding
        private Integer privateIntDerived = 1;

        @MergeByAdding
        protected Integer protectedIntDerived = 1;

        @MergeByAssertEquals
        private String privateString = "hi!";

        @Override
        public void calculateDerivedFields() {}
    }


    static class TestMergeableMetricWithStaticMembers extends TestMergeableMetric {

        public static int staticInteger = 1;

        protected static String staticString = "hi!";
        @MergeByAssertEquals
        protected static String staticEqualsString = "hi!";
        @NoMergingIsDerived
        protected static String staticDerivedString = "hi!";
        @MergingIsManual
        protected static String staticManualString = "hi!";

        @MergeByAdding
        private Integer privateIntDerived = 1;

        @MergeByAdding
        protected Integer protectedIntDerived = 1;

        @MergeByAssertEquals
        private String privateString = "hi!";

        @Override
        public void calculateDerivedFields() {}
    }



    @Test
    public void testMergingWithStaticMembers() {
        final TestMergeableMetricWithStaticMembers metric1 = new TestMergeableMetricWithStaticMembers();
        final TestMergeableMetricWithStaticMembers metric2 = new TestMergeableMetricWithStaticMembers();

        metric1.merge(metric2);
        Assert.assertEquals(metric1.boxedInt, (Integer) 2);
        Assert.assertEquals(metric1.unboxedInt, 4);
        Assert.assertEquals(metric1.privateIntDerived, (Integer) 2);
        Assert.assertEquals(metric1.protectedIntDerived, (Integer) 2);
        Assert.assertEquals(metric1.protectedIntBase, (Integer) 2);
        Assert.assertEquals(((TestMergeableMetric)metric1).privateIntBase, (Integer) 2);
        Assert.assertEquals(metric1.privateString, "hi!");

        Assert.assertEquals(metric1.staticInteger,1);
        Assert.assertEquals(metric1.staticString, "hi!");
    }

    static class TestMergeableMetricWithStaticAnnotatedMembers extends TestMergeableMetric {
        // should fail
        @MergeByAdding
        public static int staticInteger = 1;
    }


    @Test(expectedExceptions = IllegalStateException.class)
    public void testMergingWithStaticAnnotatedMembers() {
        final TestMergeableMetricWithStaticAnnotatedMembers metric1 = new TestMergeableMetricWithStaticAnnotatedMembers();
        final TestMergeableMetricWithStaticAnnotatedMembers metric2 = new TestMergeableMetricWithStaticAnnotatedMembers();

        metric1.merge(metric2);

    }


    @Test
    public void testMergingWithRestrictedMembers() {
        final TestMergeableMetricWithRestrictedMembers metric1 = new TestMergeableMetricWithRestrictedMembers();
        final TestMergeableMetricWithRestrictedMembers metric2 = new TestMergeableMetricWithRestrictedMembers();
        metric1.merge(metric2);

        Assert.assertEquals(metric1.boxedInt, (Integer) 2);
        Assert.assertEquals(metric1.unboxedInt, 4);
        Assert.assertEquals(metric1.privateIntDerived, (Integer) 2);
        Assert.assertEquals(metric1.protectedIntDerived, (Integer) 2);
        Assert.assertEquals(metric1.protectedIntBase, (Integer) 2);
        Assert.assertEquals(((TestMergeableMetric)metric1).privateIntBase, (Integer) 2);
        Assert.assertEquals(metric1.privateString, "hi!");
    }

    @Test
    public void testMergingWithRestrictedMembersUsingBaseClass() {
        final MergeableMetricBase metric1 = new TestMergeableMetricWithRestrictedMembers();
        final MergeableMetricBase metric2 = new TestMergeableMetricWithRestrictedMembers();
        metric1.merge(metric2);

        Assert.assertEquals(((TestMergeableMetricWithRestrictedMembers)metric1).boxedInt, (Integer) 2);
        Assert.assertEquals(((TestMergeableMetricWithRestrictedMembers)metric1).unboxedInt, 4);
        Assert.assertEquals(((TestMergeableMetricWithRestrictedMembers)metric1).privateIntDerived, (Integer) 2);
        Assert.assertEquals(((TestMergeableMetricWithRestrictedMembers)metric1).privateString, "hi!");
    }


    @Test
    public void testMerging() {
        final TestMergeableMetric metric1 = new TestMergeableMetric();
        final TestMergeableMetric metric2 = new TestMergeableMetric();
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

        final TestMergeableMetric metric1 = new TestMergeableMetric();
        final TestMergeableMetric metric2 = new TestMergeableMetric();
        metric1.mustBeEqualString = "goodbye";

        Assert.assertFalse(metric1.canMerge(metric2));
        metric1.merge(metric2);
    }

    @Test
    public void testMergingANull() {

        final TestMergeableMetric metric1 = new TestMergeableMetric();
        final TestMergeableMetric metric2 = new TestMergeableMetric();
        metric1.mustBeEqualString = "goodbye";
        metric2.mustBeEqualString = null;

        Assert.assertTrue(metric1.canMerge(metric2));
        metric1.merge(metric2);
        Assert.assertEquals(metric1.mustBeEqualString, "goodbye");

        Assert.assertTrue(metric2.canMerge(metric1));
        metric2.merge(metric1);
        Assert.assertEquals(metric2.mustBeEqualString, "goodbye");
    }


    @Test(expectedExceptions = IllegalStateException.class)
    public void testMergingUnequalDouble() {

        final TestMergeableMetric metric1 = new TestMergeableMetric();
        final TestMergeableMetric metric2 = new TestMergeableMetric();
        metric1.mustBeEqualDouble = 1D;

        Assert.assertFalse(metric1.canMerge(metric2));
        metric1.merge(metric2);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testMergingUnequalBoolean() {

        final TestMergeableMetric metric1 = new TestMergeableMetric();
        final TestMergeableMetric metric2 = new TestMergeableMetric();
        metric1.mustBeEqualUnboxedBoolean = true;

        Assert.assertFalse(metric1.canMerge(metric2));
        metric1.merge(metric2);
    }

    private class TestMergeableMetricIllegal extends MergeableMetricBase {
        Integer undecorated = 0;

        @Override
        public void calculateDerivedFields() {}
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testIllegalClass() {
        final TestMergeableMetricIllegal illegal1 = new TestMergeableMetricIllegal();
        final TestMergeableMetricIllegal illegal2 = new TestMergeableMetricIllegal();

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
        instance1.unboxedInt = 1;
        final TestDerivedMergableMetric instance2 = new TestDerivedMergableMetric();
        instance2.unboxedInt = 2;

        instance1.merge(instance2);
        Assert.assertEquals(instance1.unboxedInt, 3);
    }
}
