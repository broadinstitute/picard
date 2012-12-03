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

import java.util.regex.Pattern;

/**
 * Implementation of LineReader that gets its input from a String.  No charset conversion
 * is necessary because the String is in unicode.  Handles CR, LF or CRLF line termination,
 * but if asked to return the line terminator, it always comes back as LF.
 */
public class StringLineReader implements LineReader {
    private static final Pattern CRLF = Pattern.compile("\r\n");
    private final String theString;
    private int curPos = 0;
    private int lineNumber = 0;

    public StringLineReader(final String s) {
        // Simplify later processing by replacing crlf with just lf, and replacing solo cr with lf.
        // Note that String.replace(String, String) causes a regex to be used, so precompilation should be
        // the best we can do short of handling the string directly.
        this.theString = CRLF.matcher(s).replaceAll("\n").replace('\r', '\n');
    }

    /**
     * Read a line and remove the line terminator
     */
    public String readLine() {
        return readLine(false);
    }

    /**
     * Read a line and optionally include the line terminator
     *
     * @param includeTerminators
     * @return the next line from the input, with \n terminator if present and requested, or null if no more input.
     */
    private String readLine(final boolean includeTerminators) {
        if (curPos == theString.length()) {
            return null;
        }
        final int nextLfIndex = theString.indexOf('\n', curPos);
        if (nextLfIndex == -1) {
            final int startPos = curPos;
            curPos = theString.length();
            ++lineNumber;
            return theString.substring(startPos);
        }
        final int startPos = curPos;
        final int endPos = nextLfIndex + (includeTerminators? 1: 0);
        curPos = nextLfIndex + 1;
        ++lineNumber;
        return theString.substring(startPos, endPos);
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
        if (curPos == theString.length()) {
            return -1;
        }
        return theString.charAt(curPos);
    }

    public void close() {
        curPos = theString.length();
    }
}
