/*
 * Copyright (c) 2007-2010 by The Broad Institute, Inc. and the Massachusetts Institute of Technology.
 * All Rights Reserved.
 *
 * This software is licensed under the terms of the GNU Lesser General Public License (LGPL), Version 2.1 which
 * is available at http://www.opensource.org/licenses/lgpl-2.1.php.
 *
 * THE SOFTWARE IS PROVIDED "AS IS." THE BROAD AND MIT MAKE NO REPRESENTATIONS OR WARRANTIES OF
 * ANY KIND CONCERNING THE SOFTWARE, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT
 * OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE.  IN NO EVENT SHALL THE BROAD OR MIT, OR THEIR
 * RESPECTIVE TRUSTEES, DIRECTORS, OFFICERS, EMPLOYEES, AND AFFILIATES BE LIABLE FOR ANY DAMAGES OF
 * ANY KIND, INCLUDING, WITHOUT LIMITATION, INCIDENTAL OR CONSEQUENTIAL DAMAGES, ECONOMIC
 * DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER THE BROAD OR MIT SHALL
 * BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE
 * FOREGOING.
 */

package org.broad.tribble;

import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.LineReaderUtil;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broad.tribble.util.ParsingUtils;

import java.io.IOException;
import java.util.regex.Pattern;

/**
 * A convenience base class for codecs that want to read in features from ASCII lines.
 * <p/>
 * This class overrides the general decode locs for streams and presents instead
 * Strings to decode(String) and readHeader(LineReader) functions.
 *
 * @param <T> The feature type this codec reads
 */
public abstract class AsciiFeatureCodec<T extends Feature> extends AbstractFeatureCodec<T> {
    /**
     * regex used to identify what separates fields
     */
    protected Pattern splitPattern = Pattern.compile("\\t");

    protected AsciiFeatureCodec(final Class<T> myClass) {
        super(myClass);
    }

    @Override
    public Feature decodeLoc(final PositionalBufferedStream stream) throws IOException {
        String line = readLine(stream);
        try {
            return decodeLoc(line);
        } catch (RuntimeException e) {
            String msg = "\nLine: " + line;
            throw new RuntimeException(msg, e);
        }
    }

    @Override
    public T decode(final PositionalBufferedStream stream) throws IOException {
        String line = readLine(stream);
        try {
            return decode(line);
        } catch (RuntimeException e) {
            String msg = "\nLine: " + line;
            throw new RuntimeException(msg, e);
        }
    }

    @Override
    public FeatureCodecHeader readHeader(final PositionalBufferedStream stream) throws IOException {
        final LineReader br = LineReaderUtil.fromBufferedStream(stream, LineReaderUtil.LineReaderOption.SYNCHRONOUS);
        // TODO -- track header end here
        return new FeatureCodecHeader(readHeader(br), FeatureCodecHeader.NO_HEADER_END);
    }

    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     *
     * @param line the input line to decode
     * @return Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    /**
     * Decode a set of tokens as a Feature.
     * For backwards compatibility, the
     * default implementation joins by tabs, and calls {@link #decode(String)}.
     * <p/>
     * It is recommended that you override {@link #decode(String[])}
     * as well as {@link #decode(String)}
     *
     * @param tokens
     * @return
     */
    public T decode(String[] tokens) {
        String line = ParsingUtils.join("\t", tokens);
        return decode(line);
    }

    /**
     * Decode a line as a Feature.
     *
     * @param line the input line to decode
     * @return Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    abstract public T decode(String line);

    /**
     * Read and return the header, or null if there is no header.
     *
     * @return the actual header data in the file, or null if none is available
     */
    public Object readHeader(LineReader reader) {
        return null;
    }

    /**
     * Read a line of text.  A line is considered to be terminated by any one
     * of a line feed ('\n'), a carriage return ('\r'), or a carriage return
     * followed immediately by a linefeed.
     *
     * @param stream the stream to read the next line from
     * @return A String containing the contents of the line or null if the
     *         end of the stream has been reached
     */
    private String readLine(final PositionalBufferedStream stream) throws IOException {
        int linePosition = 0;

        while (true) {
            final int b = stream.read();

            if (b == -1) {
                // eof reached.  Return the last line, or null if this is a new line
                if (linePosition > 0) {
                    return new String(lineBuffer, 0, linePosition);
                } else {
                    return null;
                }
            }

            final char c = (char) (b & 0xFF);
            if (c == LINEFEED || c == CARRIAGE_RETURN) {
                if (c == CARRIAGE_RETURN && stream.peek() == LINEFEED) {
                    stream.read(); // <= skip the trailing \n in case of \r\n termination
                }

                return new String(lineBuffer, 0, linePosition);
            } else {
                // Expand line buffer size if neccessary.  Reserve at least 2 characters
                // for potential line-terminators in return string

                if (linePosition > (lineBuffer.length - 3)) {
                    char[] temp = new char[BUFFER_OVERFLOW_INCREASE_FACTOR * lineBuffer.length];
                    System.arraycopy(lineBuffer, 0, temp, 0, lineBuffer.length);
                    lineBuffer = temp;
                }

                lineBuffer[linePosition++] = c;
            }
        }
    }

    private static final int BUFFER_OVERFLOW_INCREASE_FACTOR = 2;
    private static final byte LINEFEED = (byte) ('\n' & 0xff);
    private static final byte CARRIAGE_RETURN = (byte) ('\r' & 0xff);
    private char[] lineBuffer = new char[10000];
}
