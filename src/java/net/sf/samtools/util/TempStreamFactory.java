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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * Factory class for wrapping input and output streams for temporary files.  If available, Snappy is used to
 * compress output files.  Therefore, if a temporary output file is written with an output stream obtained
 * from this class, it must be read by an input stream created by this class, otherwise a file written with
 * compression will not be read with decompression.
 */
public class TempStreamFactory {
    private static SnappyLoader snappyLoader = null;

    private static synchronized SnappyLoader getSnappyLoader() {
        if (snappyLoader == null) snappyLoader = new SnappyLoader();
        return snappyLoader;
    }

    /**
     * Wrap the given InputStream in a SnappyInputStream if available.
     * @return If Snappy is available, a SnappyInputStream wrapping inputStream.
     * If not, and bufferSize > 0, a BufferedInputStream.
     * Otherwise inputStream is returned.
     */
    public InputStream wrapTempInputStream(final InputStream inputStream, final int bufferSize) {
        if (getSnappyLoader().SnappyAvailable) {
            try {
                return getSnappyLoader().wrapInputStream(inputStream);
            } catch (Exception e) {
                throw new SAMException("Error creating SnappyInputStream", e);
            }
        } else if (bufferSize > 0) {
            return new BufferedInputStream(inputStream, bufferSize);
        } else {
            return inputStream;
        }
    }

    /**
     * Wrap the given OutputStream in a SnappyOutputStream if available.
     * @return If Snappy is available, a SnappyOutputStream wrapping outputStream.
     * If not, and bufferSize > 0, a BufferedOutputStream.
     * Otherwise outputStream is returned.
     */
    public OutputStream wrapTempOutputStream(final OutputStream outputStream, final int bufferSize) {
        if (getSnappyLoader().SnappyAvailable) {
            try {
                return getSnappyLoader().wrapOutputStream(outputStream);
            } catch (Exception e) {
                throw new SAMException("Error creating SnappyOutputStream", e);
            }
        } else if (bufferSize > 0) {
            return new BufferedOutputStream(outputStream, bufferSize);
        } else {
            return outputStream;
        }
    }
}
