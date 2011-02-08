/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
 */

package net.sf.picard.metrics;

import net.sf.picard.PicardException;
import net.sf.picard.util.FormatUtil;

import java.lang.reflect.Field;

/**
 * A base class from which all Metric classes should inherit.
 *
 * @author Tim Fennell
 */
public class MetricBase {
    /**
     * An equals method that checks equality by asserting that the classes are of the exact
     * same type and that all public fields are equal.
     *
     * @param o an instance to compare to
     * @return true if they are equal, false otherwise
     */
    public boolean equals(final Object o) {
        if (o == null) return false;
        if (o.getClass() != getClass()) return false;

        // Loop through all the fields and check that they are either
        // null in both objects or equal in both objects
        for (final Field f : getClass().getFields()) {
            try {
                final Object lhs = f.get(this);
                final Object rhs = f.get(o);

                if (lhs == null) {
                    if (rhs == null) {
                        // keep going
                    }
                    else if (rhs != null) {
                        return false;
                    }
                }
                else {
                    if (lhs.equals(rhs)) {
                        // keep going
                    }
                    else {
                        return false;
                    }
                }
            }
            catch (IllegalAccessException iae) {
                throw new PicardException("Could not read field " + f.getName() + " from a " + getClass().getSimpleName());
            }
        }

        // If we got this far all the fields are equal
        return true;
    }

    public int hashCode() {
        int result = 0;
        for (final Field f : getClass().getFields()) {
            try {
                result = 31 * result + f.get(this).hashCode();
            } catch (IllegalAccessException e) {
                throw new PicardException("Could not read field " + f.getName() + " from a " + getClass().getSimpleName());
            }
        }
        return result;
    }

    /** Converts the metric class to a human readable string. */
    public String toString() {
        final StringBuilder buffer = new StringBuilder();
        final FormatUtil formatter = new FormatUtil();

        for (final Field f : getClass().getFields()) {
            try {
                buffer.append(f.getName());
                buffer.append("\t");
                buffer.append(formatter.format(f.get(this)));
                buffer.append("\n");
            }
            catch (IllegalAccessException iae) {
                throw new PicardException("Could not read field " + f.getName() + " from a " + getClass().getSimpleName());
            }
        }

        return buffer.toString();
    }

    public boolean equals(MetricBase that) {
        for (Field field : this.getClass().getFields()) {
            try {
                if (!field.get(this).equals(field.get(that))) {
                    return false;
                }
            }
            catch (IllegalAccessException ex) {
                return false;
            }
        }
        return true;

    }
}
