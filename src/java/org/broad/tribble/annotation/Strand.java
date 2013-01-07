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
package org.broad.tribble.annotation;

/**
 * a public enum for Strand encoding
 */
public enum Strand {
    POSITIVE("+"), NEGATIVE("-"), NONE("!");  // not really sure what we should do for the NONE Enum

    // how we encode the strand information as text
    private String encoding;
    Strand(String str) {
        encoding = str;
    }

    /**
     * provide a way to take an encoding string, and produce a Strand
     * @param encoding the encoding string
     * @return a Strand object, if an appropriate one cannot be located an IllegalArg exception
     */
    public static Strand toStrand(String encoding) {
        for (Strand st : Strand.values())
            if (st.encoding.equals(encoding))
                return st;
        throw new IllegalArgumentException("Unable to match encoding to Strand enum for encoding string " + encoding);
    }

}