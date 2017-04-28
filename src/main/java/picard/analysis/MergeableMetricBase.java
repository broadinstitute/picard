/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.metrics.MetricBase;

import java.lang.annotation.Annotation;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * An extension of MetricBase that knows how to merge-by-adding fields that are appropriately annotated. It also provides an interface
 * for calculating derived fields (and an annotation that informs that said fields are derived). Finally, it also allows for an annotation
 * that suggests that a field will be used as an ID and thus merging will simply require that these fields are equal.
 *
 * merge-by-adding is only enabled for the following types: int, Integer, float, Float, double, Double, short, Short, long, Long, byte, Byte.
 * Overflow will be detected (for the short, and byte types) and an exception thrown.
 *
 * @author Yossi Farjoun
 */
abstract public class MergeableMetricBase extends MetricBase {

    /** Metrics whose values can be merged by adding. */
    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface MergeByAdding {}

    /** Metrics whose values should be equal when merging. */
    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface MergeByAssertEquals {}

    /** Metrics that are not merged, but are subsequently derived from other metrics, for example by
     * {@link #calculateDerivedFields()}. */
    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface NoMergingIsDerived {}

    /** Metrics that are merged manually in the {@link #merge(MergeableMetricBase)} ()}. Typically these metrics need
     * access to both metrics being merged. */
    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface MergingIsManual {}

    /** Metrics that are not merged. */
    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface NoMergingKeepsValue {}

    /** checks if this instance can be merged with another
     *
     * Other must have all the fields that this instance has, and
     * the fields that are annotated as MergeByAssertEquals must contain the same value
     *
     * @param other metric that will be merged into this one.
     * @return true if the other metric can be merged into this one.
     */
    public boolean canMerge(final MergeableMetricBase other) {

        try {
            for (final Field field : this.getClass().getDeclaredFields()) {
                if (field.isSynthetic()) continue;

                //try to get field from other, will throw exception if other instance doesn't have the
                field.get(other);

                final Annotation[] equalAnnotations = field.getAnnotationsByType(MergeByAssertEquals.class);
                if (equalAnnotations.length != 0) {
                    if (field.get(this) == null) return true;
                    if (field.get(other) == null) return true;
                    if (!field.get(this).equals(field.get(other))) return false;
                }
            }
        } catch (final Exception e) {
            return false;
        }
        return true;
    }

    /**
     * Merges another MergableMetricBase if possible
     *
     * @param other another MergableMetricBase instance to merge, must of the same class as this.
     * @return true if the other metric can be merged into this one.
     */
    public boolean mergeIfCan(final MergeableMetricBase other) {

        if (canMerge(other)) {
            merge(other);
            return true;
        } else {
            return false;
        }
    }

    /**
     * for a collection of MergeableMetricBase, merge them all into "this" one.
     *
     * @param many a Collection of MergeableMetricBase
     */
    public MergeableMetricBase merge(final Collection<? extends MergeableMetricBase> many) {
        many.stream().forEach(this::merge);
        calculateDerivedFields();
        return this;
    }

    /**
     * Merge another metric into this one
     *
     * @param other metric to merge into this one.
     */
    public MergeableMetricBase merge(final MergeableMetricBase other) {

        for (final Field field : getAllFields(this.getClass())) {
            if (field.isSynthetic()) continue;

            if (field.getAnnotationsByType(MergeByAdding.class).length +
                    field.getAnnotationsByType(MergeByAssertEquals.class).length +
                    field.getAnnotationsByType(NoMergingIsDerived.class).length +
                    field.getAnnotationsByType(MergingIsManual.class).length +
                    field.getAnnotationsByType(NoMergingKeepsValue.class).length == 0) {
                throw new IllegalStateException("All fields of this class must be annotated with @MergeByAdding, @NoMergingIsDerived, or @MergeByAssertEquals. " +
                        "Field " + field.getName() + " isn't annotated.");
            }

            final Annotation[] summableAnnotations = field.getAnnotationsByType(MergeByAdding.class);
            field.setAccessible(true);

            if (summableAnnotations.length != 0) {
                try {
                    if (field.getType() == Integer.class) {
                        field.set(this, (Integer) field.get(this) + (Integer) field.get(other));
                    } else if (field.getType() == int.class) {
                        field.set(this, (int) field.get(this) + (int) field.get(other));
                    } else if (field.getType() == Float.class) {
                        field.set(this, (Float) field.get(this) + (Float) field.get(other));
                    } else if (field.getType() == float.class) {
                        field.set(this, (float) field.get(this) + (float) field.get(other));
                    } else if (field.getType() == Double.class) {
                        field.set(this, (Double) field.get(this) + (Double) field.get(other));
                    } else if (field.getType() == double.class) {
                        field.set(this, (double) field.get(this) + (double) field.get(other));
                    } else if (field.getType() == Long.class) {
                        field.set(this, (Long) field.get(this) + (Long) field.get(other));
                    } else if (field.getType() == long.class) {
                        field.set(this, (long) field.get(this) + (long) field.get(other));
                    } else if (field.getType() == Byte.class) {
                        final Integer result = (Byte) field.get(this) + (Byte) field.get(other);
                        if (result > Byte.MAX_VALUE)
                            throw new IllegalArgumentException("Overflow detected in adding " + field.get(this) + " to " + field.get(other));
                        field.set(this, (byte) (int) result);
                    } else if (field.getType() == byte.class) {
                        final int result = (byte) field.get(this) + (byte) field.get(other);
                        if (result > Byte.MAX_VALUE)
                            throw new IllegalArgumentException("Overflow detected in adding " + field.get(this) + " to " + field.get(other));
                        field.set(this, (byte) result);
                    } else if (field.getType() == Short.class) {
                        final Integer result = (Short) field.get(this) + (Short) field.get(other);
                        if (result > Short.MAX_VALUE)
                            throw new IllegalArgumentException("Overflow detected in adding " + field.get(this) + " to " + field.get(other));
                        field.set(this, (Short) (short) (int) result);
                    } else if (field.getType() == short.class) {
                        final Integer result = (short) field.get(this) + (short) field.get(other);
                        if (result > Short.MAX_VALUE)
                            throw new IllegalArgumentException("Overflow detected in adding " + field.get(this) + " to " + field.get(other));
                        field.set(this, (short) (int) result);
                    } else
                        throw new IllegalArgumentException("I don't know how to MergeByAdding type " + field.getDeclaringClass().getCanonicalName() +
                                " of field " + field.getName() + "please teach me!");
                } catch (IllegalAccessException e) {
                    e.printStackTrace();
                }
            }

            final Annotation[] equalAnnotations = field.getAnnotationsByType(MergeByAssertEquals.class);
            if (equalAnnotations.length != 0) {
                try {
                    if (field.get(this) == null) {
                        field.set(this, field.get(other));
                    } else if (field.get(other) != null && !field.get(this).equals(field.get(other))) {
                        throw new IllegalStateException("Field " + field.getName() +
                                " is annotated as @MergeByAssertEquals, but found two different values: " + field.get(this) + " & " + field.get(other));
                    }
                } catch (IllegalAccessException e) {
                    e.printStackTrace();
                }
            }
        }
        return this;
    }

    private static List<Field> getAllFields(Class clazz){
        final List<Field> fields = new ArrayList<>();
        fields.addAll(Arrays.asList(clazz.getDeclaredFields()));
        final Class superClass = clazz.getSuperclass();

        if (superClass != null) fields.addAll(getAllFields(superClass));

        return fields;
    }

    /**
     * placeholder method that will calculate the derived fields from the other ones. classes that are derived from non-trivial base classes
     * should consider calling super.calculateDerivedFields() as well.
     */
    public void calculateDerivedFields() {}
}
