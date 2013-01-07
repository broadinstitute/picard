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
package org.broad.tribble;

/**
 * A class to represent a header of a feature containing file.  Specific to a codec.  All
 * codecs must return a non-null value now for their header, but the header can be the
 * empty header object or the end case be set to NO_HEADER_END.
 *
 * Note that if the headerEnd value is > 0 the readers will skip the header for the codec upfront,
 * so that decode() doesn't have to deal with potentially seeing header records in the inputstream.
 *
 * @author Mark DePristo
 * @since 5/2/12
 */
public class FeatureCodecHeader {
    /** The value of the headerEnd field when there's no header */
    public final static long NO_HEADER_END = 0;

    /** An public instance representing no header */
    public final static FeatureCodecHeader EMPTY_HEADER = new FeatureCodecHeader(null, NO_HEADER_END);

    private final Object headerValue;
    private final long headerEnd;

    /**
     * Create a FeatureCodecHeader indicating the contents of the header (can be null)
     * and the byte position in the file where the header ends (not inclusive).  headerEnd
     * should be NO_HEADER_END when no header is present.
     *
     * @param headerValue the header data read by the codec
     * @param headerEnd the position (not inclusive) of the end of the header.  1 would
     *                  mean just the first byte of the file is the header.  Zero indicates
     *                  there's no header at all
     */
    public FeatureCodecHeader(final Object headerValue, final long headerEnd) {
        if ( headerEnd < 0 ) throw new TribbleException("Header end < 0");
        this.headerValue = headerValue;
        this.headerEnd = headerEnd;
    }

    /**
     * @return the header value provided by the codec for this file
     */
    public Object getHeaderValue() {
        return headerValue;
    }

    /**
     * @return to position, not inclusive, where the header ends.  Must be >= 0
     */
    public long getHeaderEnd() {
        return headerEnd;
    }

    /**
     * @return true if we should skip some bytes to skip this header
     */
    public boolean skipHeaderBytes() {
        return getHeaderEnd() != NO_HEADER_END;
    }
}
