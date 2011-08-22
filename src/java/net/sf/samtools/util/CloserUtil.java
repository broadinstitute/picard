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
package net.sf.samtools.util;

import net.sf.samtools.util.CloseableIterator;

import java.util.List;
import java.util.Arrays;
import java.io.Closeable;
import java.io.IOException;

/**
 * Utility to close things that implement Closeable
 * WARNING: This should only be used for Closeable things open for read, because it ignores exceptions, and
 * the caller will probably want to know about exceptions when closing a file being written to, because
 * this may indicate a failure to flush.
 *
 * @author Kathleen Tibbetts
 */
public class CloserUtil {

    /**
     * Calls close() on <code>obj</code> if it implements Closeable
     *
     * @param obj   The potentially closeable object
     */
    public static void close(Object obj) {
        if (obj != null) {
            close(Arrays.asList(obj));
        }
    }

    /**
     * Calls close() on all elements of <code>objs</code> that implement Closeable
     *
     * @param objs   A list of potentially closeable objects
     */
    public static void close(List<Object> objs) {
        for (Object o : objs) {
            if (o instanceof Closeable) {
                try {
                    ((Closeable)o).close();
                }
                catch (IOException ioe) {
                    // Do nothing 
                }
            } else if (o instanceof CloseableIterator) {
                ((CloseableIterator)o).close();
            }
            else {
                try {
                    java.lang.reflect.Method m = o.getClass().getMethod("close");
                    m.invoke(o);
                }
                catch (Exception e) { /** Ignore */ }
            }
        }
    }
}
