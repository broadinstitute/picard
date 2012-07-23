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

import java.util.Collection;
import java.util.List;
import java.util.Arrays;

/**
 * Grab-bag of stateless String-oriented utilities.
 */
public class StringUtil {
    private static final byte UPPER_CASE_OFFSET = 'A' - 'a';

    /**
     * @param separator String to interject between each string in strings arg
     * @param objs List of objs to be joined
     * @return String that concatenates the result of each item's to String method for all items in objs, with separator between each of them.
     */
    public static <T> String join(final String separator, final Collection<T> objs) {
        if (objs.size() == 0) {
            return "";
        }

        boolean notFirst = false;

        final StringBuilder ret = new StringBuilder();
        for (final Object obj : objs) {
            if(notFirst) {
                ret.append(separator);
            }
            ret.append(obj.toString());
            notFirst = true;
        }
        return ret.toString();
    }

    public static <T> String join(final String separator, final T... objs) {
        final List<T> values = Arrays.asList(objs);
        return join(separator, values);
    }


    /**
     * Split the string into tokens separated by the given delimiter.  Profiling has
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
     * Split the string into tokens separated by the given delimiter.  Profiling has
     * revealed that the standard string.split() method typically takes > 1/2
     * the total time when used for parsing ascii files.
     * Note that the string is split into no more elements than tokens arg will hold, so the final tokenized
     * element may contain delimiter chars.
     *
     * @param aString  the string to split
     * @param tokens an array to hold the parsed tokens
     * @param delim  character that delimits tokens
     * @return the number of tokens parsed
     */
    public static int splitConcatenateExcessTokens(final String aString, final String[] tokens, final char delim) {

        final int maxTokens = tokens.length;
        int nTokens = 0;
        int start = 0;
        int end = aString.indexOf(delim);
        if(end < 0) {
            tokens[nTokens++] = aString;
            return nTokens;
        }
        while ((end > 0) && (nTokens < maxTokens - 1))
        {
            tokens[nTokens++] = aString.substring(start, end);
            start = end + 1;
            end = aString.indexOf(delim, start);

        }
        // Add the trailing string,  if it is not empty.
        final String trailingString = aString.substring(start);
        if (trailingString.length() > 0)
        {
            tokens[nTokens++] = trailingString;
        }
        return nTokens;
    }

    /**
     * @param b ASCII character
     * @return lowercase version of arg if it was uppercase, otherwise returns arg
     */
    public static byte toLowerCase(final byte b) {
        if (b < 'A' || b > 'Z') {
            return b;
        }
        return (byte)(b - UPPER_CASE_OFFSET);
    }

    /**
     * @param b ASCII character
     * @return uppercase version of arg if it was lowercase, otherwise returns arg
     */
    public static byte toUpperCase(final byte b) {
        if (b < 'a' || b > 'z') {
            return b;
        }
        return (byte)(b + UPPER_CASE_OFFSET);
    }

    /**
     * Converts in place all lower case letters to upper case in the byte array provided.
     */
    public static void toUpperCase(final byte[] bytes) {
        final int length = bytes.length;
        for (int i=0; i<length; ++i) {
            if (bytes[i] >= 'a' && bytes[i] <= 'z') {
                bytes[i] = (byte) (bytes[i] + UPPER_CASE_OFFSET);
            }
        }
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


    public static String intValuesToString(final int[] intVals) {
        final StringBuilder sb = new StringBuilder(intVals.length);
        if(intVals.length > 0) {
            sb.append(String.valueOf(intVals[0]));
            for(int i = 1; i < intVals.length; i++) {
                sb.append(", ");
                sb.append(String.valueOf(intVals[i]));
            }
        }

        return sb.toString();
    }

    public static String intValuesToString(final short[] shortVals) {
        final StringBuilder sb = new StringBuilder(shortVals.length);
        if(shortVals.length > 0) {
            sb.append(String.valueOf(shortVals[0]));
            for(int i = 1; i < shortVals.length; i++) {
                sb.append(", ");
                sb.append(String.valueOf(shortVals[i]));
            }
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

    @SuppressWarnings("deprecation")
    public static byte[] stringToBytes(final String s, final int offset, final int length) {
        final byte[] byteBuffer = new byte[length];
        s.getBytes(offset, offset + length, byteBuffer, 0);
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

    /**
     * Convert ASCII char to byte.
     */
    public static byte charToByte(final char c) {
        return (byte)c;
    }

    /**
     * Convert ASCII byte to ASCII char.
     */
    public static char byteToChar(final byte b) {
        return (char)(b & 0xff);
    }

    /**
     * Convert a byte array into a String hex representation.
     * @param data Input to be converted.
     * @return String twice as long as data.length with hex representation of data.
     */
    public static String bytesToHexString(final byte[] data) {
        final char[] chars = new char[2 * data.length];
        for (int i = 0; i < data.length; i++) {
            final byte b = data[i];
            chars[2*i] = toHexDigit((b >> 4) & 0xF);
            chars[2*i+1] = toHexDigit(b & 0xF);
        }
        return new String(chars);
    }

    /**
     * Convert a String containing hex characters into an array of bytes with the binary representation
     * of the hex string
     * @param s Hex string.  Length must be even because each pair of hex chars is converted into a byte.
     * @return byte array with binary representation of hex string.
     * @throws NumberFormatException
     */
    public static byte[] hexStringToBytes(final String s)  throws NumberFormatException {
        if (s.length() % 2 != 0) {
            throw new NumberFormatException("Hex representation of byte string does not have even number of hex chars: " + s);
        }
        final byte[] ret = new byte[s.length() / 2];
        for (int i = 0; i < ret.length; ++i) {
            ret[i] = (byte) ((fromHexDigit(s.charAt(i * 2)) << 4) | fromHexDigit(s.charAt(i * 2 + 1)));
        }
        return ret;
    }

    public static char toHexDigit(final int value) {
        return (char) ((value < 10) ? ('0' + value) : ('A' + value - 10));
    }

    public static int fromHexDigit(final char c) throws NumberFormatException {
        final int ret = Character.digit(c, 16);
        if (ret == -1) {
            throw new NumberFormatException("Not a valid hex digit: " + c);
        }
        return ret;
    }

    /**
     * Reverse the given string.  Does not check for null.
     * @param s String to be reversed.
     * @return New string that is the reverse of the input string.
     */
    public static String reverseString(final String s) {
        final StringBuilder sb = new StringBuilder(s);
        sb.reverse();
        return sb.toString();
    }

    /**
     * <p>Checks if a String is whitespace, empty ("") or null.</p>
     *
     * <pre>
     * StringUtils.isBlank(null)      = true
     * StringUtils.isBlank("")        = true
     * StringUtils.isBlank(" ")       = true
     * StringUtils.isBlank("sam")     = false
     * StringUtils.isBlank("  sam  ") = false
     * </pre>
     *
     * @param str  the String to check, may be null
     * @return <code>true</code> if the String is null, empty or whitespace
     */
    public static boolean isBlank(String str) {
        int strLen;
        if (str == null || (strLen = str.length()) == 0) {
            return true;
        }
        for (int i = 0; i < strLen; i++) {
            if (!Character.isWhitespace(str.charAt(i)) ) {
                return false;
            }
        }
        return true;
    }

     /* <p>Generates a string of one character to a specified length</p>
     *
     * @param c  the Character to repeat
     * @param repeatNumber the number of times to repeat the character
     * @return String with the character c repeated repeatNumber times
     */
    public static String repeatCharNTimes(char c, int repeatNumber) {
        char[] output = new char[repeatNumber];
        Arrays.fill(output, c);
        return String.valueOf(output);
    }

}
