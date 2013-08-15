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

import net.sf.samtools.util.CloserUtil;
import org.broad.tribble.readers.*;

import java.io.IOException;
import java.io.InputStream;

/**
 * A convenience base class for codecs that want to read in features from ASCII lines.
 * <p/>
 * This class overrides the general decode locs for streams and presents instead
 * Strings to decode(String) and readHeader(LineReader) functions.
 *
 * @param <T> The feature type this codec reads
 */
public abstract class AsciiFeatureCodec<T extends Feature> extends AbstractFeatureCodec<T, LineIterator> {
    protected AsciiFeatureCodec(final Class<T> myClass) {
        super(myClass);
    }
    
    @Override
    public void close(final LineIterator lineIterator) {
        CloserUtil.close(lineIterator);
    }

    @Override
    public boolean isDone(final LineIterator lineIterator) {
        return !lineIterator.hasNext();
    }

    @Override
    public LocationAware makeIndexableSourceFromStream(final InputStream bufferedInputStream) {
        final PositionalBufferedStream pbs;
        if (bufferedInputStream instanceof PositionalBufferedStream) {
            pbs = (PositionalBufferedStream) bufferedInputStream;
        } else {
            pbs = new PositionalBufferedStream(bufferedInputStream);
        }
        return new AsciiLineReaderIterator(new AsciiLineReader(pbs));
    }

    @Override
    public LineIterator makeSourceFromStream(final InputStream bufferedInputStream) {
        return new LineIteratorImpl(LineReaderUtil.fromBufferedStream(bufferedInputStream));
    }

    /** 
     * Convenience method.  Decoding in ASCII files operates line-by-line, so obviate the need to call 
     * {@link org.broad.tribble.readers.LineIterator#next()} in implementing classes and, instead, have them implement
     * {@link AsciiFeatureCodec#decode(String)}.
     */
    @Override
    public T decode(final LineIterator lineIterator) {
        return decode(lineIterator.next());
    }

    /** @see {@link AsciiFeatureCodec#decode(org.broad.tribble.readers.LineIterator)} */
    public abstract T decode(String s);

    @Override
    public FeatureCodecHeader readHeader(final LineIterator lineIterator) throws IOException {
        // TODO: Track header end here, rather than assuming there isn't one.
        return new FeatureCodecHeader(readActualHeader(lineIterator), FeatureCodecHeader.NO_HEADER_END);
    }

    /**
     * Read and return the header, or null if there is no header.
     *
     * @return the actual header data in the file, or null if none is available
     */
    abstract public Object readActualHeader(final LineIterator reader);
}
