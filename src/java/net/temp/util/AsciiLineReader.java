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

import java.io.InputStream;

/**
 * LineReader implementation that detects the OS and selects the preferred implementation for that
 * OS and delegates to it. 
 *
 * @author alecw
 */
public class AsciiLineReader implements LineReader {

    private final LineReader readerImpl;

    private static boolean useAsciiLineReaderImpl = !isLinux();

    private static boolean isMacOs() {
        return "Mac OS X".equals(System.getProperty("os.name"));
    }

    private static boolean isLinux() {
        return "Linux".equals(System.getProperty("os.name"));
    }

    /**
     * Prepare to read lines from the given input stream, with default buffer size
     * @param is need not be buffered, because this class does buffered reading
     */
    public AsciiLineReader(final InputStream is) {
        if (useAsciiLineReaderImpl) {
            readerImpl = new AsciiLineReaderImpl(is);
        } else {
            readerImpl = new BufferedLineReader(is);
        }
    }

    /**
     * Prepare to read lines from the given input stream
     * @param is need not be buffered, because this class does buffered reading
     * @param bufferSize
     */
    public AsciiLineReader(final InputStream is, final int bufferSize) {
        if (useAsciiLineReaderImpl) {
            readerImpl = new AsciiLineReaderImpl(is, bufferSize);
        } else {
            readerImpl = new BufferedLineReader(is, bufferSize);
        }
    }

    /**
     * Read a line and remove the line terminator
     *
     * @return the line read, or null if EOF has been reached.
     */
    public String readLine() {
        return readerImpl.readLine();
    }

    /**
     * @return 1-based number of line most recently read
     */
    public int getLineNumber() {
        return readerImpl.getLineNumber();
    }

    /**
     * Non-destructive one-character look-ahead.
     *
     * @return If not eof, the next character that would be read.  If eof, -1.
     */
    public int peek() {
        return readerImpl.peek();
    }

    public void close() {
        readerImpl.close();
    }
}

