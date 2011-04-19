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

import net.sf.samtools.util.CloserUtil;
import net.sf.picard.PicardException;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.MappedByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.Arrays;

/**
 * Read a .cnf (binary noise) or .cif (binary intensity) file.  A file in this format contains
 * 1 or more cycles of data for a set of clusters, with 4 values per cycle, one for each channel.
 * A file can store its values in either a byte or a short per value, but the API treats them all as shorts.
 * This class does not distinguish btw CIF and CNF files.
 *
 * @author alecw@broadinstitute.org
 */
public class ClusterIntensityFileReader {

    /**
     * The two supporte file types.
     */
    public enum FileType {cif, cnf}

    private static final byte[] IDENTIFIER = StringUtil.stringToBytes("CIF");
    private static final byte FILE_VERSION = 1;
    private static final int HEADER_SIZE = 13;
    private static final int NUM_CHANNELS = IntensityChannel.values().length;

    // Just for error reporting
    private final File file;

    /**
     * The entire file is mmapped
     */
    private final MappedByteBuffer buf;
    private final int elementSize;
    private final int firstCycle;
    private final int numCycles;
    private final int numClusters;

    // Precomputed for speed, I hope.
    private final int cycleSize;
    private final int channelSize;


    /**
     * Prepare to parse a CIF or CNF file.
     * @param file The file to be parsed.
     */
    public ClusterIntensityFileReader(final File file) {
        try {
            this.file = file;
            final FileInputStream is = new FileInputStream(this.file);
            final FileChannel channel = is.getChannel();
            final long fileSize = channel.size();
            buf = channel.map(FileChannel.MapMode.READ_ONLY, 0, fileSize);
            buf.order(ByteOrder.LITTLE_ENDIAN);
            CloserUtil.close(channel);
            CloserUtil.close(is);
            final byte[] identifierBuf = new byte[IDENTIFIER.length];
            buf.get(identifierBuf);
            if (!Arrays.equals(identifierBuf, IDENTIFIER)) {
                throw new PicardException("Cluster intensity file " + file + " contains unexpected header: " +
                        StringUtil.bytesToString(identifierBuf));
            }
            final byte fileVersion = buf.get();
            if (fileVersion != FILE_VERSION) {
                throw new PicardException("Cluster intensity file " + file + " contains unexpected version: " + fileVersion);
            }
            elementSize = buf.get();
            if (elementSize < 1 || elementSize > 2) {
                throw new PicardException("Cluster intensity file " + file + " contains unexpected element size: " + elementSize);
            }
            // convert these to unsigned
            firstCycle = buf.getShort() & 0xffff;
            numCycles = buf.getShort() & 0xffff;
            if (numCycles == 0) {
                throw new PicardException("Cluster intensity file " + file + " has zero cycles.");
            }
            numClusters = buf.getInt();
            if (numClusters < 0) {
                // It is possible for there to be no clusters in a tile.
                throw new PicardException("Cluster intensity file " + file + " has negative number of clusters: " +numClusters);
            }
        } catch (IOException e) {
            throw new PicardException("IOException opening cluster intensity file " + file, e);
        }
        cycleSize = NUM_CHANNELS * numClusters * elementSize;
        channelSize = numClusters * elementSize;
    }

    /**
     * Get the value for the given args.  Value is returned as a signed short regardless of whether storage is
     * in bytes or shorts.
     * @param cluster 0-based cluster number.
     * @param channel Which channel is desired.
     * @param cycle Absolute cycle number.  E.g. if the first cycle in the file is N, then the first value that can
     * be fetched is cycle=N
     * @return Intensity or noise (depending on whether this is a CIF or CNF file).
     */
    public short getValue(final int cluster, final IntensityChannel channel, final int cycle) {
        if (cycle < firstCycle || cycle >= firstCycle + numCycles) {
            throw new IllegalArgumentException("Requested cycle (" + cycle + ") number out of range.  First cycle=" +
                    firstCycle + "; numCycles=" + numCycles);
        }
        if (cluster < 0 || cluster >= numClusters) {
            throw new IllegalArgumentException("Requested cluster (" + cluster + ") number out of range. numClusters=" + numClusters);
        }
        final int relativeCycle = cycle - firstCycle;
        final int position = HEADER_SIZE + relativeCycle * cycleSize + channel.ordinal() * channelSize + cluster * elementSize;
        buf.position(position);
        if (elementSize == 1) {
            return buf.get();
        } else {
            return buf.getShort();
        }
    }

    public File getFile() {
        return file;
    }

    /**
     * @return The first (one-based) cycle stored in this file.
     */
    public int getFirstCycle() {
        return firstCycle;
    }

    /**
     * @return Number of clusters stored in this file.
     */
    public int getNumClusters() {
        return numClusters;
    }

    /**
     * @return Number of cycles stored in this file.
     */
    public int getNumCycles() {
        return numCycles;
    }
}
