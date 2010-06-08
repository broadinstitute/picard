/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.samtools;

/**
 * Represents the origin of a SAM record.
 *
 * @author mhanna
 * @version 0.1
 */
public class SAMFileSource {
    /**
     * The reader originating this SAM record.
     */
    private SAMFileReader mReader;

    /**
     * The point on disk from which a record originates.
     */
    private SAMFileSpan mFilePointer;

    /**
     * Create a new SAMFileSource with the given reader and file pointer.
     * @param reader reader.
     * @param filePointer File pointer.
     */
    public SAMFileSource(final SAMFileReader reader, final SAMFileSpan filePointer) {
        this.mReader = reader;
        this.mFilePointer = filePointer;
    }

    /**
     * Retrieves the reader from which this read was initially retrieved.
     * @return The reader.
     */
    public SAMFileReader getReader() {
        return mReader;
    }

    /**
     * A pointer to the region on disk from which the read originated.
     * @return A pointer within the file.
     */
    public SAMFileSpan getFilePointer() {
        return mFilePointer;
    }
}
