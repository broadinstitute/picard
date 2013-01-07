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
package org.broad.tribble.readers;

import java.io.IOException;

/**
 * User: depristo
 * Date: Oct 7, 2010
 * Time: 10:53:20 AM
 *
 * Minimal interface for an object at support getting the current position in the stream / writer / file.
 *
 * The constrain here is simple.  If you are a output stream / writer, and you've written 50 bytes to the stream,
 * then getFilePointer() should return 50 bytes.  If you are an input stream or file reader, and you've read
 * 25 bytes from the object, then getFilePointer() should return 25.
 */
public interface Positional {
    /**
     * @return the current offset, in bytes, in the stream / writer / file.
     */
    public long getPosition();

    /**
     * Is the stream done?  Equivalent to ! hasNext() for an iterator?
     * @return true if the stream has reached EOF, false otherwise
     */
    public boolean isDone() throws IOException;

    /**
     * Skip the next nBytes in the stream.
     * @param nBytes to skip, must be >= 0
     * @return the number of bytes actually skippped.
     * @throws IOException
     */
    public long skip(long nBytes) throws IOException;

    /**
     * Return the next byte in the first, without actually reading it from the stream.
     *
     * Has the same output as read()
     *
     * @return the next byte, or -1 if EOF encountered
     * @throws IOException
     */
    public int peek() throws IOException;
}
