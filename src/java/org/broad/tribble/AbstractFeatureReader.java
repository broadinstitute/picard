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
import org.broad.tribble.util.TabixUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * jrobinso
 * <p/>
 * the feature reader class, which uses indices and codecs to read in Tribble file formats.
 */
public abstract class AbstractFeatureReader<T extends Feature, SOURCE> implements FeatureReader<T> {
    // the logging destination for this source
    //private final static Logger log = Logger.getLogger("BasicFeatureSource");

    // the path to underlying data source
    String path;

    // the query source, codec, and header
    // protected final QuerySource querySource;
    protected final FeatureCodec<T, SOURCE> codec;
    protected FeatureCodecHeader header;

    private static ComponentMethods methods = new ComponentMethods();

    public static final Set<String> BLOCK_COMPRESSED_EXTENSIONS = Collections.unmodifiableSet(new HashSet<String>(Arrays.asList(".gz", ".gzip", ".bgz", ".bgzf")));

    /**
     * Calls {@link #getFeatureReader(String, FeatureCodec, boolean)} with {@code requireIndex} = true
     */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureFile, final FeatureCodec<FEATURE, SOURCE> codec) throws TribbleException {
        return getFeatureReader(featureFile, codec, true);
    }

    /**
     * {@link #getFeatureReader(String, String, FeatureCodec, boolean)} with {@code null} for indexResource
     * @throws TribbleException
     */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureResource, final FeatureCodec<FEATURE, SOURCE> codec, final boolean requireIndex) throws TribbleException {
        return getFeatureReader(featureResource, null, codec, requireIndex);
    }

    /**
     *
     * @param featureResource the feature file to create from
     * @param indexResource   the index for the feature file. If null, will auto-generate (if necessary)
     * @param codec
     * @param requireIndex    whether an index is required for this file
     * @return
     * @throws TribbleException
     */
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureResource, String indexResource, final FeatureCodec<FEATURE, SOURCE> codec, final boolean requireIndex) throws TribbleException {

        try {
            // Test for tabix index
            if (methods.isTabix(featureResource, indexResource)) {
                if ( ! (codec instanceof AsciiFeatureCodec) )
                    throw new TribbleException("Tabix indexed files only work with ASCII codecs, but received non-Ascii codec " + codec.getClass().getSimpleName());
                return new TabixFeatureReader<FEATURE, SOURCE>(featureResource, (AsciiFeatureCodec) codec);
            }
            // Not tabix => tribble index file (might be gzipped, but not block gzipped)
            else {
                return new TribbleIndexedFeatureReader<FEATURE, SOURCE>(featureResource, codec, requireIndex);
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
    public static <FEATURE extends Feature, SOURCE> AbstractFeatureReader<FEATURE, SOURCE> getFeatureReader(final String featureResource, final FeatureCodec<FEATURE, SOURCE>  codec, final Index index) throws TribbleException {
        try {
            return new TribbleIndexedFeatureReader<FEATURE, SOURCE>(featureResource, codec, index);
        } catch (IOException e) {
            throw new TribbleException.MalformedFeatureFile("Unable to create AbstractFeatureReader using feature file ", featureResource, e);
        }

    }

    protected AbstractFeatureReader(final String path, final FeatureCodec<T, SOURCE> codec) {
        this.path = path;
        this.codec = codec;
    }

    /**
     * Whether the reader has an index or not
     * Default implementation returns false
     * @return
     */
    public boolean hasIndex(){
        return false;
    }

    public static void setComponentMethods(ComponentMethods methods){
        AbstractFeatureReader.methods = methods;
    }

    /**
     * Whether a filename ends in one of the BLOCK_COMPRESSED_EXTENSIONS
     * @param fileName
     * @return
     */
    public static boolean hasBlockCompressedExtension (final String fileName) {
        for (final String extension : BLOCK_COMPRESSED_EXTENSIONS) {
            if (fileName.toLowerCase().endsWith(extension))
                return true;
        }
        return false;
    }

    /**
     * Whether the name of a file ends in one of the BLOCK_COMPRESSED_EXTENSIONS
     * @param file
     * @return
     */
    public static boolean hasBlockCompressedExtension (final File file) {
        return hasBlockCompressedExtension(file.getName());
    }

    /**
     * get the header
     *
     * @return the header object we've read-in
     */
    public Object getHeader() {
        return header.getHeaderValue();
    }

    static class EmptyIterator<T extends Feature> implements CloseableTribbleIterator<T> {
        public Iterator iterator() { return this; }
        public boolean hasNext() { return false; }
        public T next() { return null; }
        public void remove() { }
        @Override public void close() { }
    }

    public static class ComponentMethods{

        public boolean isTabix(String resourcePath, String indexPath) throws IOException{
            if(indexPath == null){
                indexPath = ParsingUtils.appendToPath(resourcePath, TabixUtils.STANDARD_INDEX_EXTENSION);
            }
            return hasBlockCompressedExtension(resourcePath) && ParsingUtils.resourceExists(indexPath);
        }
    }
}
