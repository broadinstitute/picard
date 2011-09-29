/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

package net.sf.picard.illumina.parser;

import net.sf.picard.PicardException;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Parse various formats and versions of Illumina Basecall files, and use them the to populate
 * ClusterData objects.  Use IlluminaDataProviderFactory to create one of these.
 *
 * @author alecw@broadinstitute.org
 */
public class IlluminaDataProvider {

    /**
     * Underlying parsers to use.
     */
    private final List<IlluminaParser> parsers;
    private final boolean pairedEnd;
    private final boolean barcoded;

    // These are here just to generate error messages
    private final File basecallDirectory;
    private final int lane;

    /**
     * Create an IlluminaDataProvider given a list of parsers for particular file formats
     *
     * @param barcoded True if this is a barcoded lane.
     * @param pairedEnd True if this is a paired-end lane.
     * @param parsers List of parsers that will populate the ClusterData streams.
     * @param basecallDirectory For error reporting only.
     * @param lane For error reporting only.
     */
    IlluminaDataProvider(final boolean barcoded, final boolean pairedEnd,
                         final List<IlluminaParser> parsers,
                         final File basecallDirectory, final int lane) {
        this.barcoded = barcoded;
        this.basecallDirectory = basecallDirectory;
        this.lane = lane;
        this.pairedEnd = pairedEnd;
        this.parsers = parsers;
    }

    public boolean hasNext() {
        if (parsers.isEmpty()) {
            return false;
        }
        final boolean ret = parsers.get(0).hasNext();
        for (int i = 1; i < parsers.size(); ++i) {
            if (parsers.get(i).hasNext() != ret) {
                throw new PicardException("Unequal length basecall files in " + basecallDirectory + ", lane " + lane);
            }
        }
        return ret;
    }

    public ClusterData next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        final ClusterData cluster = new ClusterData();
        cluster.setFirstEnd(new ReadData());
        if (pairedEnd) {
            cluster.setSecondEnd(new ReadData());
        }
        if (barcoded) {
            cluster.setBarcodeRead(new ReadData());
        }
        for (final IlluminaParser parser: parsers) {
            parser.next(cluster);
        }
        return cluster;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /** Jump so that the next record returned will be from the specified tile. */
    public void seekToTile(final int oneBasedTileNumber) {
        for (final IlluminaParser parser: parsers) {
            parser.seekToTile(oneBasedTileNumber);
        }
    }

}
