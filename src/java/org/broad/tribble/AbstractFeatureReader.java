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

import org.broad.tribble.index.Index;
import org.broad.tribble.util.ParsingUtils;

import java.io.IOException;
import java.util.Iterator;

/**
 * jrobinso
 * <p/>
 * the feature reader class, which uses indices and codecs to read in Tribble file formats.
 */
public abstract class AbstractFeatureReader<T extends Feature> implements FeatureReader<T> {
    // the logging destination for this source
    //private final static Logger log = Logger.getLogger("BasicFeatureSource");

    // the path to underlying data source
    String path;

    // the query source, codec, and header
    // protected final QuerySource querySource;
    protected final FeatureCodec codec;
    protected FeatureCodecHeader header;

    // A hook for the future, when we might allow clients to specify this.


    public static final AbstractFeatureReader getFeatureReader(String featureFile, FeatureCodec codec) throws TribbleException {
        return getFeatureReader(featureFile, codec, true);
    }

    /**
     * factory for unknown file type,  could be ascii, binary, or could be tabix, or something else.
     *
     * @param featureResource the feature file to create from
     * @param codec           the codec to read features with
     */
    public static final AbstractFeatureReader getFeatureReader(String featureResource, FeatureCodec codec, boolean requireIndex) throws TribbleException {

        try {
            // Test for tabix index
            if (featureResource.endsWith(".gz") &&
                    ParsingUtils.resourceExists(featureResource + ".tbi")) {
                if ( ! (codec instanceof AsciiFeatureCodec) )
                    throw new TribbleException("Tabix indexed files only work with ASCII codecs, but received non-Ascii codec " + codec.getClass().getSimpleName());
                return new TabixFeatureReader(featureResource, (AsciiFeatureCodec)codec);
            }
            // Not tabix => tribble index file (might be gzipped, but not block gzipped)
            else {
                return new TribbleIndexedFeatureReader(featureResource, codec, requireIndex);
            }
        } catch (IOException e) {
            throw new TribbleException.MalformedFeatureFile("Unable to create BasicFeatureReader using feature file ", featureResource, e);
        } catch (TribbleException e) {
            e.setSource(featureResource);
            throw e;
        }
    }

    /**
     * Return a reader with a supplied index.
     *
     * @param featureResource the path to the source file containing the features
     * @param codec used to decode the features
     * @param index index of featureResource
     * @return a reader for this data
     * @throws TribbleException
     */
    public static final AbstractFeatureReader getFeatureReader(String featureResource, FeatureCodec codec, Index index) throws TribbleException {
        try {
            return new TribbleIndexedFeatureReader(featureResource, codec, index);
        } catch (IOException e) {
            throw new TribbleException.MalformedFeatureFile("Unable to create AbstractFeatureReader using feature file ", featureResource, e);
        }

    }

    protected AbstractFeatureReader(String path, FeatureCodec codec) {
        this.path = path;
        this.codec = codec;
    }


    /**
     * get the header
     *
     * @return the header object we've read-in
     */
    public Object getHeader() {
        return header.getHeaderValue();
    }

    static class EmptyIterator<T extends Feature> implements CloseableTribbleIterator {
        public Iterator iterator() { return this; }
        public boolean hasNext() { return false; }
        public Object next() { return null; }
        public void remove() { }
        @Override public void close() { }
    }
}
