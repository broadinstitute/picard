/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import java.util.*;

/**
 * Small utility methods for dealing with collection classes.
 */
public class CollectionUtil {

    public static <T> List<T> makeList (final T... list) {
        final List<T> result = new ArrayList<T>();
        for (final T item : list) {
            result.add(item);
        }
        return result;
    }

    public static <T> Set<T> makeSet (final T... list) {
        final Set<T> result = new HashSet<T>();
        for (final T item : list) {
            result.add(item);
        }
        return result;
    }

    /** Construct a string by toString()ing each item in the collection with inBetween between each item. */
    public static String join(final Collection<?> items, final String inBetween) {
        final StringBuilder builder = new StringBuilder();
        for (final Object item : items) {
            if (builder.length() > 0) builder.append(inBetween);
            builder.append(item);
        }

        return builder.toString();
    }
}
