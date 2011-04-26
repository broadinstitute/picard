/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.util;

import java.io.File;
import java.io.InputStream;

/**
 * Parser for tab-delimited files
 *
 * @author Kathleen Tibbetts
 */
public class TabbedInputParser extends BasicInputParser {

    /**
     * Constructor
     *
     * @param stream  The input stream(s) to parse
     */
    public TabbedInputParser(boolean treatGroupedDelimitersAsOne, InputStream... stream) {
        super(treatGroupedDelimitersAsOne, stream);
    }

    /**
     * Constructor
     *
     * @param file  The file(s) to parse
     */
    public TabbedInputParser(boolean treatGroupedDelimitersAsOne, File... file) {
        super(treatGroupedDelimitersAsOne, file);
    }

    /**
     * Determines whether a given character is a delimiter
     *
     * @param b the character to evaluate
     * @return  true if <code>b</code> is a delimiter; otherwise false
     */
    protected boolean isDelimiter(byte b) {
        return b == '\t';
    }
}
