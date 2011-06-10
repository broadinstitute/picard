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

import java.io.InputStream;
import java.io.OutputStream;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;

/**
 * If Snappy is available, obtain single-arg ctors for SnappyInputStream and SnappyOutputStream.
 */
public class SnappyLoader {
    private static final int SNAPPY_BLOCK_SIZE = 32768;  // keep this as small as can be without hurting compression ratio.
    private final Constructor<InputStream> SnappyInputStreamCtor;
    private final Constructor<OutputStream> SnappyOutputStreamCtor;
    public final boolean SnappyAvailable;

    // Force Snappy-java code to be loaded into executable jars.
    private final SnappyInputStream ignoreMe = null;

    // Force bcel to load Snappy.
    //private static final Class SnappyClass = SnappyInputStream.class;

    private static final boolean DefaultVerbosity = Boolean.valueOf(System.getProperty("snappy.loader.verbosity", "false"));

    public SnappyLoader() {
        this(DefaultVerbosity);
    }

    /**
     * Constructs a new SnappyLoader which will check to see if snappy is available in the JVM/library path.
     * @param verbose if true output a small number of debug messages to System.err
     */
    public SnappyLoader(final boolean verbose) {
        Constructor<InputStream> inputStreamCtor = null;
        Constructor<OutputStream> outputStreamCtor = null;

        if (java.lang.Boolean.valueOf(System.getProperty("snappy.disable", "false"))) {
            System.err.println("Snappy is disabled via system property.");
        }
        else {
            try {
                final Class<InputStream> snappyInputStreamClass = (Class<InputStream>)Class.forName("org.xerial.snappy.SnappyInputStream");
                final Class<OutputStream> snappyOutputStreamClass = (Class<OutputStream>)Class.forName("org.xerial.snappy.SnappyOutputStream");
                inputStreamCtor = snappyInputStreamClass.getConstructor(InputStream.class);
                outputStreamCtor = snappyOutputStreamClass.getConstructor(OutputStream.class, Integer.TYPE);
            }
            catch (NoSuchMethodException e) { /* Do nothing. */ }
            catch (ClassNotFoundException e) { /* Do nothing. */ }
        }

        this.SnappyInputStreamCtor = inputStreamCtor;
        this.SnappyOutputStreamCtor = outputStreamCtor;

        if (this.SnappyInputStreamCtor != null && this.SnappyOutputStreamCtor != null) {
            // Don't try to call any Snappy code until classes have been found via reflection above.
            if (!LoadSnappy.load()) {
                if (verbose) System.err.println("Snappy dll failed to load.");
                SnappyAvailable = false;
            }
            else {
                if (verbose) System.err.println("Snappy stream classes loaded.");
                SnappyAvailable = true;
            }
        }
        else {
            if (verbose) System.err.println("Snappy stream classes not loaded.");
            SnappyAvailable = false;
        }
    }

    /** Wrap an InputStream in a SnappyInputStream. If Snappy is not available will throw an exception. */
    public InputStream wrapInputStream(final InputStream inputStream) {
        try {
            return SnappyInputStreamCtor.newInstance(inputStream);
        } catch (Exception e) {
            throw new SAMException("Error instantiating SnappyInputStream", e);
        }
    }

    /** Wrap an InputStream in a SnappyInputStream. If Snappy is not available will throw an exception. */
    public OutputStream wrapOutputStream(final OutputStream outputStream) {
        try {
            return SnappyOutputStreamCtor.newInstance(outputStream, SNAPPY_BLOCK_SIZE);
        } catch (Exception e) {
            throw new SAMException("Error instantiating SnappyOutputStream", e);
        }
    }
}
