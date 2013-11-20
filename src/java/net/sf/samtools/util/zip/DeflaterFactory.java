/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package net.sf.samtools.util.zip;

import net.sf.samtools.Defaults;
import net.sf.samtools.SAMException;

import java.lang.reflect.Constructor;
import java.util.zip.Deflater;

/**
 * Create zlib-based Deflater if JNI library and other require libraries are available, otherwise create standard
 * JDK Deflater.
 * Java 7 has its own Deflater implementation (libzip.so).  This is almost as fast as a zlib-based Deflater, so in general
 * there isn't a compelling reason to use zlib.  However, Intel has created a hardware-assisted zlib implementation
 * as part of their IPP (Integrated Performance Primitives) package that can run significantly faster on some Intel
 * hardware.  We have seen compression times reduced by 13% to 33% depending on particular hardware, and hope that
 * newer Intel processors will be even better.
 *
 * Note that this class will no longer be necessary once Java 8 is required, because JDK 8 will use zlib instead
 * of libzip implementation.
 */
public class DeflaterFactory {

    private static Constructor<IntelDeflater> intelDeflaterConstructor;

    static {
        try {
            if (Defaults.TRY_USE_INTEL_DEFLATER) {
                final Class<IntelDeflater> clazz = (Class<IntelDeflater>) Class.forName("net.sf.samtools.util.zip.IntelDeflater");
                intelDeflaterConstructor = clazz.getConstructor(Integer.TYPE, Boolean.TYPE);
            }
        } catch (ClassNotFoundException e) {
            intelDeflaterConstructor = null;
        } catch (NoSuchMethodException e) {
            intelDeflaterConstructor = null;
        } catch (UnsatisfiedLinkError e) {
            intelDeflaterConstructor = null;
        }
    }

    public static Deflater makeDeflater(final int compressionLevel, final boolean nowrap) {
        if (intelDeflaterConstructor != null) {
            try {
                return intelDeflaterConstructor.newInstance(compressionLevel, nowrap);
            } catch (Exception e) {
                throw new SAMException("Exception constructing IntelDeflater", e);
            }
        } else {
            return new Deflater(compressionLevel, nowrap);
        }
    }

    public static boolean usingIntelDeflater() {
        return intelDeflaterConstructor != null;
    }
}
