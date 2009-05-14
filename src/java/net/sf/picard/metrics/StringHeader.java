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

package net.sf.picard.metrics;

import net.sf.samtools.util.StringUtil;

/**
 * A simple header who's data type is a single String.  Should not be used for anything other
 * than comments or descriptive text.
 *
 * @author Tim Fennell
 */
public class StringHeader implements Header {
    private String value;

    /** Default constructor. */
    public StringHeader() {}

    /** Constructor that uses the supplied value as the value of the header. */
    public StringHeader(String value) {
        setValue(value);
    }

    public void parse(String in) { value = in.trim(); }
    public String toString() { return value; }

    public String getValue() { return value; }
    public void setValue(String value) { this.value = StringUtil.assertCharactersNotInString(value, '\n'); }

    /** Checks equality on the value of the header. */
    public boolean equals(Object o) {
        if (o != null && o instanceof StringHeader) {
            StringHeader that = (StringHeader) o;
            if (this.value == null) {
                return that.value == null;
            }
            else {
                return this.value.equals(that.value);
            }
        }
        else {
            return false;
        }
    }
}
