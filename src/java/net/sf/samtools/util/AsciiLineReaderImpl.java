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
import java.io.InputStream;

/**
 * This is the old AsciiLineReader, except that lineBuffer is stored as chars rather than bytes.  It appears
 * that the String ctor that takes byte[] is significantly slower than the one that takes char[].
 *
 * Note that in general in Picard, conversion between ASCII bytes and chars is encapsulated in StringUtil.
 * This class is an exception because there are significant performance implications.  It simply casts bytes
 * to chars.
 *
 * On some platforms (e.g. Linux), BufferedLineReader is faster than this class.  If you use AsciiLineReader
 * instead of this class, it will detect the OS and delegate to the preferred implementation for that platform.
 *
 * @author alecw@broadinstitute.org
 */
public class AsciiLineReaderImpl implements LineReader {
    private static final char LINEFEED = '\n';
    private static final char CARRIAGE_RETURN = '\r';

    private final InputStream is;
    private byte[] buffer;
    private int nextChar;
    private int nChars;
    // Allocate this only once, despite the fact that it is essentially a local variable of readLine()
    private char[] lineBuffer = new char[10000];

    private int lineNumber = 0;

    /**
     * Prepare to read lines from the given input stream, with default buffer size
     * @param is need not be buffered, because this class does buffered reading
     */
    public AsciiLineReaderImpl(final InputStream is) {
        this(is, 512000);
    }

    /**
     * Prepare to read lines from the given input stream
     * @param is need not be buffered, because this class does buffered reading
     * @param bufferSize
     */
    public AsciiLineReaderImpl(final InputStream is, final int bufferSize) {
        this.is = is;
        buffer = new byte[bufferSize];
        nextChar = nChars = 0;
    }

    /**
     * Read a line of text.  A line is considered to be terminated by any one
     * of a line feed ('\n'), a carriage return ('\r'), or a carriage return
     * followed immediately by a linefeed.
     *
     * @return     A String containing the contents of the line or null if the
     *             end of the stream has been reached.  Terminator is not included in
     *             the returned string.
     */
    public String readLine() {
        return readLine(false);
    }

    /**
     * Read a line of text.  A line is considered to be terminated by any one
     * of a line feed ('\n'), a carriage return ('\r'), or a carriage return
     * followed immediately by a linefeed.
     *
     * @param      includeTerminators  If true, the line-termination characters
     *             are included in the returned string.
     *
     * @return     A String containing the contents of the line or null if the
     *             end of the stream has been reached
     */
    public String readLine(final boolean includeTerminators){
        int linePosition = 0;

        while (true)
        {
            if (nChars == -1)
            {
                return null;
            }

            // Refill buffer if neccessary
            if (nextChar == nChars)
            {
                fill();
                if (nextChar == nChars || nChars == -1)
                {
                    // eof reached.  Return the last line, or null if this is a new line
                    if (linePosition > 0)
                    {
                        ++lineNumber;
                        return new String(lineBuffer, 0, linePosition);
                    } else
                    {
                        return null;
                    }
                }
            }


            final char c = (char) (buffer[nextChar++] & 0xFF);
            if (c == LINEFEED || c == CARRIAGE_RETURN)
            {

                if (includeTerminators)
                {
                    lineBuffer[linePosition++] = c;
                    if (c == CARRIAGE_RETURN && peek() == LINEFEED)
                    {
                        lineBuffer[linePosition++] = c;
                        nextChar++; // <= to account for the '\n' we just ate
                    }
                }
                else {
                    if (c == CARRIAGE_RETURN && peek() == LINEFEED)
                    {
                        nextChar++; // <= skip the trailing \n in case of \r\n termination
                    }

                }
                ++lineNumber;
                return new String(lineBuffer, 0, linePosition);
            } else
            {
                // Expand line buffer size if neccessary.  Reservce at least 2 characters
                // for potential line-terminators in return string

                if (linePosition > (lineBuffer.length - 3))
                {
                    final char[] temp = new char[lineBuffer.length + 1000];
                    System.arraycopy(lineBuffer, 0, temp, 0, lineBuffer.length);
                    lineBuffer = temp;
                }

                lineBuffer[linePosition++] = c;
            }
        }
    }

    /**
     * @return the line number of the most recently read line.  After reading the first line, this method
     * return 1.
     */
    public int getLineNumber() {
        return lineNumber;
    }

    /**
     * Peek ahead one character, filling from the underlying stream if neccessary.
     *
     * @return peeked character.
     * @throws java.io.IOException
     */
    public int peek(){
        // Refill buffer if neccessary
        if (nextChar == nChars)
        {
            fill();
            if (nextChar == nChars)
            {
                // eof reached.
                return -1;
            }
        }
        return buffer[nextChar];

    }

    private void fill() {
        try {
            nChars = is.read(buffer);
            nextChar = 0;
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    public void close()  {
        try {
            is.close();
        } catch (IOException e) {
            // Ignore exception
        }
    }
}
