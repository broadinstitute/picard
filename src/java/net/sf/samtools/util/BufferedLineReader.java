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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;

/**
 * Implementation of LineReader that is a thin wrapper around BufferedReader.  On Linux, this is faster
 * than AsciiLineReaderImpl.  If you use AsciiLineReader rather than this class, it will detect the OS
 * and delegate to the preferred implementation.
 *
 * @author alecw@broadinstitute.org
 */
public class BufferedLineReader implements LineReader {

    private final BufferedReader reader;
    private int lineNumber = 0;
    private String peekedLine;

    public BufferedLineReader(final InputStream is) {
        this(is, IOUtil.STANDARD_BUFFER_SIZE);
    }

    public BufferedLineReader(final InputStream is, final int bufferSize) {
        reader = new BufferedReader(new InputStreamReader(is, Charset.forName("UTF-8")), bufferSize);
    }

    /**
     * Read a line and remove the line terminator
     *
     * @return the line read, or null if EOF has been reached.
     */
    public String readLine() {
        ++lineNumber;
        try {
            final String ret;
            if (peekedLine != null) {
                ret = peekedLine;
                peekedLine = null;
            } else {
            ret = reader.readLine();
            }
            return ret;
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * @return 1-based number of line most recently read
     */
    public int getLineNumber() {
        return lineNumber;
    }

    /**
     * Non-destructive one-character look-ahead.
     *
     * @return If not eof, the next character that would be read.  If eof, -1.
     */
    public int peek() {
        if (peekedLine == null) {
            try {
                peekedLine = reader.readLine();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }
        if (peekedLine == null) {
            return -1;
        }
        if (peekedLine.length() == 0) {
            return '\n';
        }
        return peekedLine.charAt(0);
    }

    public void close() {
        peekedLine = null;
        try {
            reader.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }
}
