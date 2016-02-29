package picard.analysis.replicates;

import htsjdk.samtools.metrics.MetricBase;

import java.lang.annotation.Annotation;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.lang.reflect.Field;

/**
 * An extension of MetricBase that knows how to merge-by-adding fields that are appropriately annotated. It also provides an interface
 * for calculating derived fields (and an annotation that informs that said fields are derived). Finally, it also allows for an annotation
 * that suggests that a field will be used as an ID and thus merging will simply require that these fields are equal.
 *
 * merge-by-adding is only enabled for the following types: int, Integer, float, Float, double, Double, short, Short, long, Long, byte, Byte.
 * Overflow will be detected (for the short, and byte types) and an exception thrown.
 */
public class MergeableMetricBase extends MetricBase {

    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface MergeByAdding {
    }

    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface MergeByAssertEquals {
    }

    @Retention(RetentionPolicy.RUNTIME)
    @Target(ElementType.FIELD)
    protected @interface NoMergingIsDerived {
    }

    /** checks if this instance can be merged with another
     *
     * Other must have all the fields that this instance has, and
     * the fields that are annotated as MergeByAssertEquals must contain the same value
     *
     * @param other
     * @return
     */
    public boolean canMerge(final MergeableMetricBase other) {


        try {
            for (final Field field : this.getClass().getDeclaredFields()) {
                if (field.isSynthetic()) continue;

                //try to get field from other, will throw exception if other instance doesn't have the
                field.get(other);

                final Annotation[] equalAnnotations = field.getAnnotationsByType(MergeByAssertEquals.class);
                if (equalAnnotations.length != 0) {
                    if (!field.get(this).equals(field.get(other))) {
                        return false;
                    }
                }
            }
        } catch (final Exception e) {
            return false;
        }
        return true;
    }

    /** Merges another MergableMetricBase if possible
     *
     * @param other another MergableMetricBase instance to merge, must of the same class as this.
     * @return
     */
    public boolean mergeIfCan(final MergeableMetricBase other) {
        if(canMerge(other)) {
            merge(other);
            return true;
        }
        else {
            return false;
        }
    }

    public void merge(final MergeableMetricBase other) {


        for (final Field field : this.getClass().getDeclaredFields()) {
            if(field.isSynthetic()) continue;

            if (field.getAnnotationsByType(MergeByAdding.class).length +
                    field.getAnnotationsByType(MergeByAssertEquals.class).length +
                    field.getAnnotationsByType(NoMergingIsDerived.class).length == 0) {
                throw new IllegalStateException("All fields of this class must be annotated with @MergeByAdding, @NoMergingIsDerived, or @MergeByAssertEquals. " +
                        "Field " + field.getName() + " isn't annotated.");
            }

            final Annotation[] summableAnnotations = field.getAnnotationsByType(MergeByAdding.class);
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
                    if (!field.get(this).equals(field.get(other))) {
                        throw new IllegalStateException("Field " + field.getName() +
                                " is annotated as @MergeByAssertEquals, but found two different values: " + field.get(this) + " & " + field.get(other));
                    }
                } catch (IllegalAccessException e) {
                    e.printStackTrace();
                }

            }
        }
    }

     public void calculateDerivedFields(){}
}
