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

import org.broad.tribble.readers.LocationAware;

import java.io.IOException;
import java.io.InputStream;

/**
 * The base interface for classes that read in features.
 * <p/>
 * FeatureCodecs have to implement two key methods:
 * <p/>
 * {@link #readHeader(SOURCE)} - Reads the header, provided a {@link SOURCE} pointing at the beginning of the source input.
 * {@link #decode(SOURCE)} - Reads a {@link Feature} record, provided a {@link SOURCE} pointing at the beginning of a record within the 
 * source input.
 * <p/>
 * Note that it's not safe to carry state about the {@link SOURCE} within the codec.  There's no guarantee about its  state between calls.
 *
 * @param <FEATURE_TYPE> The type of {@link Feature} this codec generates
 * @param <SOURCE> The type of the data source this codec reads from
 */
public interface FeatureCodec<FEATURE_TYPE extends Feature, SOURCE> {
    /**
     * Decode a line to obtain just its FeatureLoc for indexing -- contig, start, and stop.
     *
     * @param source the input stream from which to decode the next record
     * @return Return the FeatureLoc encoded by the line, or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    public Feature decodeLoc(final SOURCE source) throws IOException;

    /**
     * Decode a single {@link Feature} from the {@link SOURCE}, reading no further in the underlying source than beyond that feature.
     *
     * @param source the input stream from which to decode the next record
     * @return Return the Feature encoded by the line,  or null if the line does not represent a feature (e.g. is
     *         a comment)
     */
    public FEATURE_TYPE decode(final SOURCE source) throws IOException;

    /**
     * Read and return the header, or null if there is no header.
     * 
     * Note: Implementers of this method must be careful to read exactly as much from {@link SOURCE} as needed to parse the header, and no 
     * more. Otherwise, data that might otherwise be fed into parsing a {@link Feature} may be lost.
     *
     * @param source the source from which to decode the header
     * @return header object
     */
    public FeatureCodecHeader readHeader(final SOURCE source) throws IOException;

    /**
     * This function returns the object the codec generates.  This is allowed to be Feature in the case where
     * conditionally different types are generated.  Be as specific as you can though.
     * <p/>
     * This function is used by reflections based tools, so we can know the underlying type
     *
     * @return the feature type this codec generates.
     */
    public Class<FEATURE_TYPE> getFeatureType();

    /**
     * Generates a reader of type {@link SOURCE} appropriate for use by this codec from the generic input stream.  Implementers should
     * assume the stream is buffered.
     */
    public SOURCE makeSourceFromStream(final InputStream bufferedInputStream);

    /**
     * Generates a {@link LocationAware} reader of type {@link SOURCE}.  Like {@link #makeSourceFromStream(java.io.InputStream)}, except
     * the {@link LocationAware} compatibility is required for creating indexes.
     * 
     * Implementers of this method must return a type that is both {@link LocationAware} as well as {@link SOURCE}.  Note that this 
     * requirement cannot be enforced via the method signature due to limitations in Java's generic typing system.  Instead, consumers
     * should cast the call result into a {@link SOURCE} when applicable.
     */
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream);

    /** Adapter method that assesses whether the provided {@link SOURCE} has more data. True if it does, false otherwise. */
    public boolean isDone(final SOURCE source);

    /** Adapter method that closes the provided {@link SOURCE}. */
    public void close(final SOURCE source);

    /**
     * This function returns true iff the File potentialInput can be parsed by this
     * codec.
     * <p/>
     * There is an assumption that there's never a situation where two different Codecs
     * return true for the same file.  If this occurs, the recommendation would be to error out.
     * <p/>
     * Note this function must never throw an error.  All errors should be trapped
     * and false returned.
     *
     * @param path the file to test for parsability with this codec
     * @return true if potentialInput can be parsed, false otherwise
     */
    public boolean canDecode(final String path);
}
