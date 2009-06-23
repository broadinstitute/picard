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
 * Fast replacement for BufferedReader that assumes that bytes can be converted to chars simply by casting.
 *
 * @author jrobinso
 */
public class AsciiLineReader implements LineReader {
    private static final byte LINEFEED = (byte)('\n' & 0xff);
    private static final byte CARRIAGE_RETURN = (byte)('\r' & 0xff);

    private final InputStream is;
    private byte[] buffer;
    private int nextChar;
    private int nChars;
    // Allocate this only once, despite the fact that it is essentially a local variable of readLine()
    private byte[] lineBuffer = new byte[1000];

    private int lineNumber = 0;

    /**
     * Prepare to read lines from the given input stream, with default buffer size
     * @param is need not be buffered, because this class does buffered reading
     */
    public AsciiLineReader(final InputStream is) {
        this(is, 512000);
    }

    /**
     * Prepare to read lines from the given input stream
     * @param is need not be buffered, because this class does buffered reading
     * @param bufferSize
     */
    public AsciiLineReader(final InputStream is, final int bufferSize) {
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
                        return StringUtil.bytesToString(lineBuffer, 0, linePosition);
                    } else
                    {
                        return null;
                    }
                }
            }


            final byte b = buffer[nextChar++];
            if (b == LINEFEED || b == CARRIAGE_RETURN)
            {

                if (includeTerminators)
                {
                    lineBuffer[linePosition++] = b;
                    if (b == CARRIAGE_RETURN && peek() == LINEFEED)
                    {
                        lineBuffer[linePosition++] = b;
                        nextChar++; // <= to account for the '\n' we just ate
                    }
                }
                else {
                    if (b == CARRIAGE_RETURN && peek() == LINEFEED)
                    {
                        nextChar++; // <= skip the trailing \n in case of \r\n termination
                    }
                    
                }
                ++lineNumber;
                return StringUtil.bytesToString(lineBuffer, 0, linePosition);
            } else
            {
                // Expand line buffer size if neccessary.  Reservce at least 2 characters
                // for potential line-terminators in return string

                if (linePosition > (lineBuffer.length - 3))
                {
                    final byte[] temp = new byte[lineBuffer.length + 100];
                    System.arraycopy(lineBuffer, 0, temp, 0, lineBuffer.length);
                    lineBuffer = temp;
                }

                lineBuffer[linePosition++] = b;
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

