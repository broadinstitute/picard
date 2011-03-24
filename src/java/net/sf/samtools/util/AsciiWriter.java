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

import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;

/**
 * Fast (I hope) buffered Writer that converts char to byte merely by casting, rather than charset conversion.
 */
public class AsciiWriter extends Writer {

    private final OutputStream os;
    // Buffer size has not been tuned.
    private final byte[] buffer = new byte[IOUtil.STANDARD_BUFFER_SIZE];
    private int numBytes;

    /**
     * @param os need not be buffered as this class buffers
     */
    public AsciiWriter(final OutputStream os) {
        this.os = os;
        numBytes = 0;
    }

    /**
     * flushes and closes underlying OutputStream.
     */
    public void close() throws IOException {
        flush();
        os.close();
    }

    /**
     * flushes underlying OutputStream
     */
    public void flush() throws IOException {
        os.write(buffer, 0, numBytes);
        numBytes = 0;
        os.flush();
    }

    /**
     * All other Writer methods vector through this, so this is the only one that must be overridden.
     */
    public void write(final char[] chars, int offset, int length) throws IOException {
        while (length > 0) {
            final int charsToConvert = Math.min(length, buffer.length - numBytes);
            StringUtil.charsToBytes(chars, offset, charsToConvert, buffer, numBytes);
            numBytes += charsToConvert;
            offset += charsToConvert;
            length -= charsToConvert;
            if (numBytes == buffer.length) {
                os.write(buffer, 0, numBytes);
                numBytes = 0;
            }
        }
    }
}
