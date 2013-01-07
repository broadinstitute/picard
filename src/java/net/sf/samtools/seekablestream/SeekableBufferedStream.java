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

package net.sf.samtools.seekablestream;

import java.io.BufferedInputStream;
import java.io.IOException;

/**
 * A wrapper class to provide buffered read access to a SeekableStream.  Just wrapping such a stream with
 * a BufferedInputStream will not work as it does not support seeking.  In this implementation a
 * seek call is delegated to the wrapped stream, and the buffer reset. 
 */
public class SeekableBufferedStream extends SeekableStream {

    public static final int DEFAULT_BUFFER_SIZE = 512000;

    final private int bufferSize;
    final SeekableStream wrappedStream;
    BufferedInputStream bufferedStream;
    long position;

    public SeekableBufferedStream(SeekableStream httpStream, int bufferSize) {
        this.bufferSize = bufferSize;
        this.wrappedStream = httpStream;
        this.position = 0;
        bufferedStream = new BufferedInputStream(wrappedStream, bufferSize);
    }
    public SeekableBufferedStream(SeekableStream httpStream) {
        this(httpStream, DEFAULT_BUFFER_SIZE);
    }

    public long length() {
        return wrappedStream.length();
    }

    public void seek(long position) throws IOException {
        this.position = position;
        wrappedStream.seek(position);
        bufferedStream = new BufferedInputStream(wrappedStream, bufferSize);
    }

    public int read() throws IOException {
        int b = bufferedStream.read();
        position++;
        return b;
    }

    public int read(byte[] buffer, int offset, int length) throws IOException {
        int nBytesRead = bufferedStream.read(buffer, offset, length);
        if (nBytesRead > 0) {
            position += nBytesRead;
        }
        return nBytesRead;
    }

    public void close() throws IOException {
        wrappedStream.close();
    }

    public boolean eof() throws IOException {
        return position >= wrappedStream.length();
    }

    @Override
    public String getSource() {
        return wrappedStream.getSource();
    }

    @Override
    public long position() throws IOException {
        return position;
    }
}
