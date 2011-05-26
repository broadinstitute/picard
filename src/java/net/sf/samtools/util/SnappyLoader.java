/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

import net.sf.samtools.SAMException;
import org.xerial.snappy.LoadSnappy;
import org.xerial.snappy.SnappyInputStream;
import org.xerial.snappy.SnappyOutputStream;
import sun.font.TrueTypeFont;

import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

/**
 * If Snappy is available, obtain single-arg ctors for SnappyInputStream and SnappyOutputStream.
 */
public class SnappyLoader {
    private final Constructor<InputStream> SnappyInputStreamCtor;
    private final Constructor<OutputStream> SnappyOutputStreamCtor;
    public final boolean SnappyAvailable;

    // Force bcel to load Snappy.
    //private static final Class SnappyClass = SnappyInputStream.class;

    private static final boolean DefaultVerbosity = Boolean.valueOf(System.getProperty("snappy.loader.verbosity", "false"));

    public SnappyLoader() {
        this(DefaultVerbosity);
    }

    public SnappyLoader(boolean verbose) {
        Constructor<InputStream> inputStreamCtor = null;
        Constructor<OutputStream> outputStreamCtor = null;
        if (java.lang.Boolean.valueOf(System.getProperty("snappy.disable", "false"))) {
            System.err.println("Snappy is disabled via system property.");
        } else {
            try {
                final Class<InputStream> snappyInputStreamClass = (Class<InputStream>)Class.forName("org.xerial.snappy.SnappyInputStream");
                final Class<OutputStream> snappyOutputStreamClass = (Class<OutputStream>)Class.forName("org.xerial.snappy.SnappyOutputStream");
                inputStreamCtor = snappyInputStreamClass.getConstructor(InputStream.class);
                outputStreamCtor = snappyOutputStreamClass.getConstructor(OutputStream.class);
            } catch (NoSuchMethodException e) {
            } catch (ClassNotFoundException e) {
            }
        }
        SnappyInputStreamCtor = inputStreamCtor;
        SnappyOutputStreamCtor = outputStreamCtor;
        if (SnappyInputStreamCtor != null && SnappyOutputStreamCtor != null) {
            // Don't try to call any Snappy code until classes have been found via reflection above.
            if (!LoadSnappy.load()) {
                if (verbose) System.err.println("Snappy dll failed to load.");
                SnappyAvailable = false;
            } else {
                SnappyAvailable = true;
                if (verbose) System.err.println("Snappy stream classes loaded.");
            }
        } else {
            SnappyAvailable = false;
            if (verbose) System.err.println("Snappy stream classes not loaded.");
        }
    }

    public InputStream wrapInputStream(InputStream inputStream) {
        try {
            return SnappyInputStreamCtor.newInstance(inputStream);
        } catch (Exception e) {
            throw new SAMException("Error instantiating SnappyInputStream", e);
        }
    }

    public OutputStream wrapOutputStream(OutputStream outputStream) {
        try {
            return SnappyOutputStreamCtor.newInstance(outputStream);
        } catch (Exception e) {
            throw new SAMException("Error instantiating SnappyOutputStream", e);
        }
    }
}
