/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
package picard.sam.util;

/**
 * Common functions for quickly parsing strings. Used for parsing the tile and coordinates from the read names
 */
public class ReadNameParsingUtils {

    /**
     * Single pass method to parse the read name for the default regex.  Examines the last three fields as split by the delimiter.
     */
    public static int getRapidDefaultReadNameRegexSplit(final String readName, final char delim, final int[] tokens) {
        int tokensIdx = 0;
        int prevIdx = 0;
        int numFields = 1;
        for (int i = 0; i < readName.length(); i++) {
            if (readName.charAt(i) == delim) numFields++;
        }
        int startOffset = numFields - 2 - 1; // zero-based (ex. 7 -> 4, 5 -> 2)
        if (startOffset < 0) return -1;
        int endOffset = startOffset + 2; // zero-based
        for (int i = 0; i < readName.length(); i++) {
            if (readName.charAt(i) == delim) {
                if (startOffset <= tokensIdx && tokensIdx <= endOffset)
                    tokens[tokensIdx - startOffset] = rapidParseInt(readName.substring(prevIdx, i)); // only fill in 2-4 inclusive for 5 fields
                tokensIdx++;
                prevIdx = i + 1;
            }
        }
        if (prevIdx < readName.length()) {
            if (startOffset <= tokensIdx && tokensIdx <= endOffset)
                tokens[tokensIdx - startOffset] = rapidParseInt(readName.substring(prevIdx, readName.length())); // only fill in 2-4 inclusive
            tokensIdx++;
        }
        return tokensIdx;
    }

    /**
     * Very specialized method to rapidly parse a sequence of digits from a String up until the first
     * non-digit character.
     */
    public static int rapidParseInt(final String input) {
        final int len = input.length();
        int val = 0;
        int i = 0;
        boolean isNegative = false;

        if (0 < len && '-' == input.charAt(0)) {
            i = 1;
            isNegative = true;
        }

        for (; i < len; ++i) {
            final char ch = input.charAt(i);
            if (Character.isDigit(ch)) {
                val = (val * 10) + (ch - 48);
            } else {
                break;
            }
        }

        if (isNegative) val = -val;

        return val;
    }
}
