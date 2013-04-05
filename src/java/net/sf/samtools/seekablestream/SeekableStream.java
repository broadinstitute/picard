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
package net.sf.samtools.seekablestream;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

public abstract class SeekableStream extends InputStream {

    public abstract long length();

    public abstract long position() throws IOException;

    public abstract void seek(long position) throws IOException;

    public abstract int read(byte[] buffer, int offset, int length) throws IOException;

    public abstract void close() throws IOException;

    public abstract boolean eof() throws IOException;

    /**
     * @return String representation of source (e.g. URL, file path, etc.), or null if not available.
     */
    public abstract String getSource();

    /**
     * Read enough bytes to fill the input buffer.
     * @param b
     * @throws EOFException If EOF is reached before buffer is filled
     */
    public void readFully(byte b[]) throws IOException {
        int len = b.length;
        if (len < 0){
            throw new IndexOutOfBoundsException();
        }
        int n = 0;
        while (n < len) {
            int count = read(b, n, len - n);
            if (count < 0){
                throw new EOFException();
            }
            n += count;
        }
    }

}
