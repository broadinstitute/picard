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

/**
 * Grab-bag of stateless String-oriented utilities.
 */
public class StringUtil {
    /**
     *
     * @param separator String to interject between each string in strings arg
     * @param strings List of strings to be joined.
     * @return String that concatenates each item of strings arg, with separator btw each of them.
     */
    public static String join(final String separator, final String... strings) {
        if (strings.length == 0) {
            return "";
        }
        final StringBuilder ret = new StringBuilder(strings[0]);
        for (int i = 1; i < strings.length; ++i) {
            ret.append(separator);
            ret.append(strings[i]);
        }
        return ret.toString();
    }

    /**
     * Split the string into tokesn separated by the given delimiter.  Profiling has
     * revealed that the standard string.split() method typically takes > 1/2
     * the total time when used for parsing ascii files.
     * Note that if tokens arg is not large enough to all the tokens in the string, excess tokens are discarded.
     *
     * @param aString  the string to split
     * @param tokens an array to hold the parsed tokens
     * @param delim  character that delimits tokens
     * @return the number of tokens parsed
     */
    public static int split(final String aString, final String[] tokens, final char delim) {

        final int maxTokens = tokens.length;
        int nTokens = 0;
        int start = 0;
        int end = aString.indexOf(delim);
        if(end < 0) {
            tokens[nTokens++] = aString;
            return nTokens;
        }
        while ((end > 0) && (nTokens < maxTokens))
        {
            tokens[nTokens++] = aString.substring(start, end);
            start = end + 1;
            end = aString.indexOf(delim, start);

        }
        // Add the trailing string,  if there is room and if it is not empty.
        if (nTokens < maxTokens)
        {
            final String trailingString = aString.substring(start);
            if (trailingString.length() > 0)
            {
                tokens[nTokens++] = trailingString;
            }
        }
        return nTokens;
    }

    /**
     * @param b ASCII character
     * @return uppercase version of arg if it was lowercase, otherwise returns arg 
     */
    public static byte toUpperCase(final byte b) {
        if (b < 'a' || b > 'z') {
            return b;
        }
        return (byte)(b + ('A' - 'a'));
    }


    /**
     * Checks that a String doesn't contain one or more characters of interest.
     *
     * @param illegalChars the String to check
     * @param chars the characters to check for
     * @return String the input String for convenience
     * @throws IllegalArgumentException if the String contains one or more of the characters
     */
    public static String assertCharactersNotInString(final String illegalChars, final char... chars) {
        for (final char illegalChar : illegalChars.toCharArray()) {
            for (final char ch: chars) {
                if (illegalChar == ch) {
                    throw new IllegalArgumentException("Supplied String contains illegal character '" + illegalChar + "'.");
                }
            }
        }

        return illegalChars;
    }

    /**
     * Return input string with newlines inserted to ensure that all lines
     * have length <= maxLineLength.  if a word is too long, it is simply broken
     * at maxLineLength.  Does not handle tabs intelligently (due to implementer laziness).
     */
    public static String wordWrap(final String s, final int maxLineLength) {
        final String[] lines = s.split("\n");
        final StringBuilder sb = new StringBuilder();
        for (final String line: lines) {
            if (sb.length() > 0) {
                sb.append("\n");
            }
            sb.append(wordWrapSingleLine(line, maxLineLength));
        }
        if (s.endsWith("\n")) {
            sb.append("\n");
        }
        return sb.toString();
    }

    public static String wordWrapSingleLine(final String s, final int maxLineLength) {
        if (s.length() <= maxLineLength) {
            return s;
        }
        final StringBuilder sb = new StringBuilder();
        int startCopyFrom = 0;
        while (startCopyFrom < s.length()) {
            int lastSpaceIndex = startCopyFrom;
            int i;
            // Find break point (if it exists)
            for (i = startCopyFrom; i < s.length() && i - startCopyFrom < maxLineLength; ++i) {
                if (Character.isWhitespace(s.charAt(i))) {
                    lastSpaceIndex = i;
                }
            }
            if (i - startCopyFrom < maxLineLength) {
                lastSpaceIndex = i;
            }
            // Include any trailing whitespace
            for (; lastSpaceIndex < s.length() && Character.isWhitespace(s.charAt(lastSpaceIndex)); ++lastSpaceIndex) {}
            if (sb.length() > 0) {
                sb.append("\n");
            }
            // Handle situation in which there is no word break.  Just break the word in the middle.
            if (lastSpaceIndex == startCopyFrom) {
                lastSpaceIndex = i;
            }
            sb.append(s.substring(startCopyFrom, lastSpaceIndex));
            startCopyFrom = lastSpaceIndex;
        }
        return sb.toString();
    }

    ////////////////////////////////////////////////////////////////////
    // The following methods all convert btw bytes and Strings, without
    // using the Java character set mechanism.
    ////////////////////////////////////////////////////////////////////

    public static String bytesToString(final byte[] data) {
        if (data == null) {
            return null;
        }
        return bytesToString(data, 0, data.length);
    }

    @SuppressWarnings("deprecation")
    public static String bytesToString(final byte[] buffer, final int offset, final int length) {
/*
        The non-deprecated way, that requires allocating char[]
        final char[] charBuffer = new char[length];
        for (int i = 0; i < length; ++i) {
            charBuffer[i] = (char)buffer[i+offset];
        }
        return new String(charBuffer);
*/
        return new String(buffer, 0, offset, length);
    }

    @SuppressWarnings("deprecation")
    public static byte[] stringToBytes(final String s) {
/*
        The non-deprecated way, that requires allocating char[]
        final byte[] byteBuffer = new byte[s.length()];
        final char[] charBuffer = s.toCharArray();
        for (int i = 0; i < charBuffer.length; ++i) {
            byteBuffer[i] = (byte)(charBuffer[i] & 0xff);
        }
        return byteBuffer;
*/
        final byte[] byteBuffer = new byte[s.length()];
        s.getBytes(0, byteBuffer.length, byteBuffer, 0);
        return byteBuffer;
    }

    // This method might more appropriately live in BinaryCodec, but all the byte <=> char conversion
    // should be in the same place.
    public static String readNullTerminatedString(final BinaryCodec binaryCodec) {
        final StringBuilder ret = new StringBuilder();
        for (byte b = binaryCodec.readByte(); b != 0; b = binaryCodec.readByte()) {
            ret.append((char)(b & 0xff));
        }
        return ret.toString();
    }

    /**
     * Convert chars to bytes merely by casting
     * @param chars input chars
     * @param charOffset where to start converting from chars array
     * @param length how many chars to convert
     * @param bytes where to put the converted output
     * @param byteOffset where to start writing the converted output.
     */
    public static void charsToBytes(final char[] chars, final int charOffset, final int length,
                                    final byte[] bytes, final int byteOffset) {
        for (int i = 0; i < length; ++i) {
            bytes[byteOffset + i] = (byte)chars[charOffset + i];
        }
    }

}
