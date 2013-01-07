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

import org.broad.tribble.readers.AsciiLineReader;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.PositionalBufferedStream;
import org.broad.tribble.util.ParsingUtils;

import java.io.IOException;
import java.util.regex.Pattern;

/**
 * A convenience base class for codecs that want to read in features from ASCII lines.
 *
 * This class overrides the general decode locs for streams and presents instead
 * Strings to decode(String) and readHeader(LineReader) functions.
 *
 * @param <T> The feature type this codec reads
 */
public abstract class AsciiFeatureCodec<T extends Feature> extends AbstractFeatureCodec<T> {
    /** A cached line reader we will use for decode and decodeLoc() */
    private final AsciiLineReader lineReader = new AsciiLineReader();

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
        try{
            return decodeLoc(line);
        }catch (RuntimeException e){
            String msg = "\nLine: " + line;
            throw new RuntimeException(msg, e);
        }
    }

    @Override
    public T decode(final PositionalBufferedStream stream) throws IOException {
        String line = readLine(stream);
        try{
            return decode(line);
        }catch (RuntimeException e){
            String msg = "\nLine: " + line;
            throw new RuntimeException(msg, e);
        }
    }

    @Override
    public FeatureCodecHeader readHeader(final PositionalBufferedStream stream) throws IOException {
        final AsciiLineReader br = new AsciiLineReader(stream);
        // TODO -- track header end here
        return new FeatureCodecHeader(readHeader(br), FeatureCodecHeader.NO_HEADER_END);
    }

    private final String readLine(final PositionalBufferedStream stream) throws IOException {
        return lineReader.readLine(stream);
    }

    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     *
     * @param line the input line to decode
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decodeLoc(String line) {
        return decode(line);
    }

    /**
     * Decode a set of tokens as a Feature.
     * Default implementation joins by tabs, and calls
     * {@link #decode(String)}, for backwards compatibility.
     *
     * It is recommended that you override {@link #decode(String[])}
     * as well as {@link #decode(String)}
     * @param tokens
     * @return
     */
    public T decode(String[] tokens){
        String line = ParsingUtils.join("\t", tokens);
        return decode(line);
    }

    /**
     * Decode a line as a Feature.
     *
     * @param line the input line to decode
     * @return  Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     * a comment)
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

}
