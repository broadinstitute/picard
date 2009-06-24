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
package net.sf.picard.io;

import net.sf.picard.PicardException;
import net.sf.samtools.util.CloserUtil;

import java.io.InputStream;
import java.io.IOException;

/**
 * Line-oriented InputStream reader that uses one buffer for disk buffering and line-termination-finding,
 * in order to improve performance.
 *
 * Implementation detail: All public methods must leave the input buffer in a non-empty state, unless at EOF.
 *
 * @author alecw@broadinstitute.org
 */
public class FastLineReader {
    private InputStream in;
    private byte[] fileBuffer = new byte[512000];
    // Next byte to read in fileBuffer
    private int nextByte = 0;
    // Number of bytes in fileBuffer
    private int numBytes = 0;
    private boolean atEof;

    public FastLineReader(final InputStream in) {
        this.in = in;
        ensureBufferNotEmpty();
    }

    /**
     * @return true if input is exhausted
     */
    public boolean eof() {
        return atEof;
    }

    /**
     * @return peeks at the next byte in the stream and returns true if it is CR or LF.  Returns false if EOF or
     * next byte is not CR or LF.
     */
    public boolean atEoln() {
        return ensureBufferNotEmpty() && (fileBuffer[nextByte] == '\n' || fileBuffer[nextByte] == '\r');
    }

    /**
     * Advance over any EOLN chars (CR or LF)
     * @return true if saw one or more CR or LFs
     */
    public boolean skipNewlines() {
        boolean sawEoln = false;
        while (atEoln()) {
            sawEoln = true;
            ++nextByte;
        }
        return sawEoln;
    }

    public void close() {
        CloserUtil.close(in);
        in = null;
        fileBuffer = null;
    }

    /**
     * @return Next byte from the input.  Do not call if at EOF.
     */
    public byte getByte() {
        final byte ret = peekByte(); 
        ++nextByte;
        ensureBufferNotEmpty();
        return ret;
    }

    /**
     * @return Next byte from the input, without advancing over that byte.  Do not call if at EOF.
     */
    public byte peekByte() {
        if (eof()) {
            throw new IllegalStateException("Cannot getByte() if EOF.");
        }
        return fileBuffer[nextByte];
    }

    /**
     * Read from input until input is exhausted, EOLN is seen, or output buffer is filled
     * @param outputBuffer where to put bytes read
     * @param startOutputIndex where to start putting bytes read
     * @return number of bytes read
     */
    public int readToEndOfOutputBufferOrEoln(final byte[] outputBuffer, final int startOutputIndex) {
        boolean sawNewline;
        int totalGrabbed = 0;
        do {
            if (!ensureBufferNotEmpty()) {
                break;
            }
            final int startInputIndex = nextByte;
            sawNewline = advanceToEobOrEoln();
            int lengthOfChunk = nextByte - startInputIndex;

            // Roll back if went past the amount that can be stored in the output buffer.
            // Assumption is that lines are relatively short so this won't happen very often.
            if (lengthOfChunk > outputBuffer.length - startOutputIndex) {
                lengthOfChunk = outputBuffer.length - startOutputIndex;
                nextByte = startInputIndex + lengthOfChunk;
            }
            System.arraycopy(fileBuffer, startInputIndex, outputBuffer, startOutputIndex + totalGrabbed, lengthOfChunk);
            totalGrabbed += lengthOfChunk;
        } while (!sawNewline && totalGrabbed < outputBuffer.length - startOutputIndex);
        ensureBufferNotEmpty();
        return totalGrabbed;
    }

    /**
     * Advance nextByte to end of currently-buffered input or to line terminator
     * @return true if saw a line terminator
     */
    private boolean advanceToEobOrEoln() {
        while (nextByte < numBytes) {
            if (atEoln()) {
                return true;
            }
            ++nextByte;
        }
        return false;
    }

    /**
     * Ensure that fileBuffer has at least one byte available in it.  Potentially wipes out
     * what is in fileBuffer so everything from fileBuffer[0..nextByte] should already have been pulled out.
     * @return false if EOF, else true
     */
    private boolean ensureBufferNotEmpty() {
        try {
            if (nextByte < numBytes) {
                return true;
            }
            nextByte = 0;
            numBytes = in.read(fileBuffer);
            atEof = (numBytes < 1);
            return !atEof;
        } catch (IOException e) {
            throw new PicardException("Exception reading InputStream", e);
        }
    }

}
