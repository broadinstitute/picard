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

import org.broad.tribble.readers.PositionalBufferedStream;

import java.io.IOException;

/**
 * the base interface for classes that read in features.
 *
 * FeatureCodecs have to implement two key methods:
 *
 * readHeader() => starting from the first line of the file, read the full header, if any, and
 *   return it, as well as the position in the file where the header stops and records begin.  The
 *   contract with the readers is that the header and decoders see fresh streams, so you can
 *   safely read into the first record of the file looking for the header.  It's always why
 *   you need to return the file position (via getPosition() in the stream) of the last
 *   header record.
 *
 * decode(stream) => parse out the next record, and return it.  Decode is always called on a
 *   fresh stream after the header is read.
 *
 * Note that it's not safe to carry state about the PositionalBufferedStream arguments here.  There's
 * no guarentee on the state of the stream between calls.
 *
 * The canDecode is used to determine if a file can be decoded by this codec.  Just open up the
 * file and check if it can be decoded with this codec.
 *
 * @param <T> The feature type this codec reads
 */
public interface FeatureCodec<T extends Feature> {
    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     *
     *
     * @param stream the input stream from which to decode the next record
     * @return  Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public Feature decodeLoc(final PositionalBufferedStream stream) throws IOException;

    /**
     * Decode a line as a Feature.
     *
     *
     * @param stream the input stream from which to decode the next record
     * @return  Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     * a comment)
     */
    public T decode(final PositionalBufferedStream stream) throws IOException;

    /**
     * Like {@link #decodeLoc(org.broad.tribble.readers.PositionalBufferedStream)}, but the source is a string representing a line.
     * @param line
     * @return
     */
    public Feature decodeLoc(final String line);

    /**
     * Like {@link #decode(org.broad.tribble.readers.PositionalBufferedStream)}, but the source is a string representing a line.
     * @param line
     * @return
     */
    public T decode(final String line);
    
    /**
     * Read and return the header, or null if there is no header.
     *
     *
     *
     * @param stream the input stream from which to decode the header
     * @return header object
     */
    public FeatureCodecHeader readHeader(final PositionalBufferedStream stream) throws IOException;

    /**
     * This function returns the object the codec generates.  This is allowed to be Feature in the case where
     * conditionally different types are generated.  Be as specific as you can though.
     *
     * This function is used by reflections based tools, so we can know the underlying type
     *
     * @return the feature type this codec generates.
     */
    public Class<T> getFeatureType();

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     *
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     *
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(final String path);
}
