/*
 * The MIT License
 *
 * Copyright (c) 2019 The Broad Institute
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

package picard.arrays.illumina;

import org.apache.commons.io.IOUtils;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * A class to provide methods for accessing Illumina Infinium Data Files.
 */
abstract class InfiniumDataFile {

    private String identifier;
    private int numberOfEntries;
    private int fileVersion;

    final DataInputStream stream;

    InfiniumDataFile(final DataInputStream stream, final boolean cacheStream) throws IOException {
        if (cacheStream) {
            final byte[] data = readStreamIntoByteArray(stream);

            // Don't need to buffer this one because it is sitting in memory
            this.stream = new DataInputStream(new ByteArrayInputStream(data));
        } else {
            this.stream = stream;
        }
    }

    /**
     * Utility method for reading a string of data. (Reads from the current offset)
     *
     * @return The parsed string.
     * @throws java.io.IOException is thrown when there is a problem reading the stream.
     */
    String parseString() throws IOException {
        final String dataString;
        final byte strLen = stream.readByte();
        if (strLen != 0) {
            final byte[] stringBytes = new byte[strLen];
            final int bytesRead = stream.read(stringBytes);
            if (bytesRead != stringBytes.length) {
                throw new IOException("Did not fully read string. Read " + bytesRead + " out of "
                        + stringBytes.length + ".");
            }
            dataString = new String(byteArrayToCharArray(stringBytes));
        } else {
            dataString = "";
        }
        return dataString;
    }

    private char[] byteArrayToCharArray(final byte[] stringBytes) {
        final char[] chars = new char[stringBytes.length];
        for (int i = 0; i < chars.length; i++) {
            chars[i] = (char) stringBytes[i];
        }
        return chars;
    }

    /**
     * Utility method for parsing an array of byte values.
     *
     * @param toc The table of content record for parsing the byte values.
     * @return An array of byte values for the given TOC.
     * @throws java.io.IOException is thrown when there is a problem reading the stream.
     */
    byte[] parseByteArray(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        final int arrayLen = Integer.reverseBytes(stream.readInt());
        final byte[] byteArray = new byte[arrayLen];
        for (int i = 0; i < arrayLen; i++) {
            byteArray[i] = stream.readByte();
        }
        return byteArray;
    }

    /**
     * Utility method to convert an unsigned short to an int.
     *
     * @param bytes The byte array representing the unsigned short.
     *              (Java has no unsigned values which is why we promote it to an int)
     * @return The converted int.
     */
    private int unsignedShortToInt(final byte[] bytes) {
        int integer = 0;
        integer |= bytes[1] & 0xFF;
        integer <<= 8;
        integer |= bytes[0] & 0xFF;
        return integer;
    }

    /**
     * Utility method to convert a byte array to a float value.
     *
     * @param bytes The byte array representing the float value.
     * @return The converted float.
     */
    private float byteArrayToFloat(final byte[] bytes) {
        int tempInt = ((0xff & bytes[0])
                | ((0xff & bytes[1]) << 8)
                | ((0xff & bytes[2]) << 16)
                | ((0xff & bytes[3]) << 24));
        return Float.intBitsToFloat(tempInt);
    }

    private static final int FLOAT_BYTES_LENGTH = 4;


    /**
     * Utility method for parsing a float value. (Reads from current offset)
     *
     * @return The parsed float value.
     * @throws java.io.IOException is thrown when there is a problem reading the stream.
     */
    float parseFloat() throws IOException {
        final byte[] floatBytes = new byte[FLOAT_BYTES_LENGTH];
        stream.readFully(floatBytes);
        return byteArrayToFloat(floatBytes);
    }

    /**
     * Utility method for parsing a float value.
     *
     * @param toc The table of contents record to parse the float from.
     * @return The parsed float value.
     * @throws java.io.IOException is thrown when there is a problem reading the stream.
     */
    float parseFloat(final InfiniumFileTOC toc) throws IOException {
        final byte[] floatBytes = new byte[FLOAT_BYTES_LENGTH];
        stream.skipBytes(toc.getOffset());
        stream.readFully(floatBytes);
        return byteArrayToFloat(floatBytes);
    }

    /**
     * Utility method for parsing an array of unsigned short values.
     *
     * @param toc The table of content record for parsing the unsigned short values.
     * @return An array of unsigned short values for the given TOC.
     * (Java has no unsigned values which is why we promote it to an int)
     * @throws java.io.IOException is thrown when there is a problem reading the stream.
     */
    int[] parseUnsignedShortArray(final InfiniumFileTOC toc)
            throws IOException {
        final byte[] shortBytes = new byte[2];
        stream.skipBytes(toc.getOffset());
        final int arrayLen = Integer.reverseBytes(stream.readInt());
        int[] unsignedShortArray = new int[arrayLen];
        for (int i = 0; i < arrayLen; i++) {
            stream.readFully(shortBytes);
            unsignedShortArray[i] = unsignedShortToInt(shortBytes);
        }
        return unsignedShortArray;
    }

    int parseShort(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        return readShort();
    }

    int readShort() throws IOException {
        final byte[] shortBytes = new byte[2];
        stream.readFully(shortBytes);
        return unsignedShortToInt(shortBytes);
    }

    int parseInt(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        return Integer.reverseBytes(stream.readInt());
    }

    /**
     * Utility method for parsing a string.
     *
     * @param toc The table of contents information for this string.
     * @return The parsed string from the given table of contents.
     * @throws java.io.IOException thrown when there is an error reading the data stream.
     */
    String parseString(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        return parseString();
    }

    /**
     * Utility method for parsing an array of float values.
     *
     * @param toc The table of content record for parsing the float values.
     * @return An array of float values for the given TOC.
     * @throws java.io.IOException is thrown when there is a problem reading the stream.
     */
    float[] parseFloatArray(final InfiniumFileTOC toc) throws IOException {
        stream.skipBytes(toc.getOffset());
        final int arrayLen = Integer.reverseBytes(stream.readInt());

        final float[] floatArray = new float[arrayLen];
        for (int i = 0; i < arrayLen; i++) {
            floatArray[i] = parseFloat();
        }
        return floatArray;
    }

    /**
     * This method is used to avoid a null assignment that was being done for some reason in the caller. By limiting
     * the scope of the stream, it will accomplish this without the need for a null.
     *
     * @param streamToCache The input stream that will be cached.
     * @return The byte array
     * @throws java.io.IOException Errors reading the stream
     */
    private byte[] readStreamIntoByteArray(final InputStream streamToCache) throws IOException {
        try (final ByteArrayOutputStream outStream = new ByteArrayOutputStream()) {
            //read the entire inputstream into memory
            IOUtils.copy(streamToCache, outStream);
            return outStream.toByteArray();
        }
    }

    public String getIdentifier() {
        return identifier;
    }

    public void setIdentifier(final String identifier) {
        this.identifier = identifier;
    }

    private int getNumberOfEntries() {
        return numberOfEntries;
    }

    void setNumberOfEntries(int numberOfEntries) {
        this.numberOfEntries = numberOfEntries;
    }

    int getFileVersion() {
        return fileVersion;
    }

    void setFileVersion(final int fileVersion) {
        this.fileVersion = fileVersion;
    }

    /**
     * Parse the table of contents.
     *
     * @return The TOC
     * @throws java.io.IOException Any errors reading in the TOC
     */
    InfiniumFileTOC[] getTableOfContents() throws IOException {

        final InfiniumFileTOC[] tableOfContents = new InfiniumFileTOC[getNumberOfEntries()];

        //read in the table of contents... order them by offset so that we
        //only have to traverse the input stream once.
        for (int i = 0; i < getNumberOfEntries(); i++) {
            final InfiniumFileTOC toc = new InfiniumFileTOC();
            toc.setTableOfContentsId(Short.reverseBytes(stream.readShort()));
            toc.setOffset(Integer.reverseBytes(stream.readInt()));
            tableOfContents[i] = toc;
        }

        return tableOfContents;
    }

    int parseInt() throws IOException {
        return Integer.reverseBytes(stream.readInt());
    }

    void skipFloats(int numFloats) throws IOException {
        stream.skipBytes(numFloats * FLOAT_BYTES_LENGTH);
    }

    void skipFloat() throws IOException {
        skipFloats(1);
    }

    void skipBoolean() throws IOException {
        stream.skipBytes(1);
    }

    void skipString() throws IOException {
        byte strLen = stream.readByte();
        stream.skipBytes(strLen);
    }
}
