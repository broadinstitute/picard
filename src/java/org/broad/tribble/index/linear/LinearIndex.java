/*
 * Copyright (c) 2009-2010 by The Broad Institute, Inc.
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

package org.broad.tribble.index.linear;

import org.broad.tribble.TribbleException;
import org.broad.tribble.index.AbstractIndex;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Index defined by dividing the genome by chromosome, then each chromosome into bins of fixed width (in
 * genomic coordinates).   Features are allocated to bins by start position.  The longest feature in each
 * recorded and used to adjust the start position of a query to include all bins that might have a feature
 * that overlaps the query interval.  This works well for feature sets of approximately homogeneous length,
 * or whose longest feature is on the order of the bin width or less.
 * <p/>
 * magicNumber      integer
 * type             integer
 * version          integer
 * filename         null terminated character array
 * filesize         long
 * lastModified     long
 * md5              String
 * flags            integer
 * <p/>
 * ------  LINEAR INDEX
 * nChromosomes     integer
 */
public class LinearIndex extends AbstractIndex {

    // NOTE: To debug uncomment the System.getProperty and recompile.
    public static final double MAX_FEATURES_PER_BIN = Double.valueOf(System.getProperty("MAX_FEATURES_PER_BIN", "100"));
    public static final int INDEX_TYPE = IndexType.LINEAR.fileHeaderTypeIdentifier;

    private final static int MAX_BIN_WIDTH = 1 * 1000 * 1000 * 1000; //  widths must be less than 1 billion

    // 1MB: we will no merge bins with any features in them beyond this size, no matter how sparse, per chromosome
    private static final long MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX = Long.valueOf(System.getProperty("MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX", "1024000"));

    public static boolean enableAdaptiveIndexing = true;

    /**
     * Initialize using the specified {@code indices}
     * @param indices
     * @param featureFile
     */
    public LinearIndex(final List<ChrIndex> indices, final File featureFile) {
        super(featureFile.getAbsolutePath());
        for (final ChrIndex index : indices)
            chrIndices.put(index.getName(), index);
    }

    private LinearIndex(final LinearIndex parent, final List<ChrIndex> indices) {
        super(parent);
        for (final ChrIndex index : indices)
            chrIndices.put(index.getName(), index);
    }

    /**
     * Initialize with default parameters
     * @param featureFile File for which this is an index
     */
    public LinearIndex(final String featureFile) {
        super(featureFile);
    }

    /**
     * Load from file.
     * @param inputStream This method assumes that the input stream is already buffered as appropriate.
     */
    public LinearIndex(final InputStream inputStream) throws IOException {
        final LittleEndianInputStream dis = new LittleEndianInputStream(inputStream);
        validateIndexHeader(INDEX_TYPE, dis);
        read(dis);
    }

    public boolean isCurrentVersion() {
        if (!super.isCurrentVersion()) return false;

        // todo fixme nasty hack to determine if this is an old style V3 linear index (without nFeaturesPerBin)
        for (final org.broad.tribble.index.ChrIndex chrIndex : chrIndices.values())
            if (((ChrIndex) chrIndex).OLD_V3_INDEX)
                return false;

        return true;
    }

    @Override
    protected int getType() {
        return INDEX_TYPE;
    }

    public List<String> getSequenceNames() {
        return (chrIndices == null ? Collections.EMPTY_LIST :
                Collections.unmodifiableList(new ArrayList<String>(chrIndices.keySet())));
    }

    @Override
    public Class getChrIndexClass() {
        return ChrIndex.class;
    }


    /**
     * Blocks are organized as a simple flat list:
     * <p/>
     * Block 0
     * Block 1
     * Block 2
     * <p/>
     * There's a constant bin width, so that each block corresponds to a specific interval
     * over the genome based on its index, as in:
     * <p/>
     * Block 0: (0 - binWidth]
     * Block 1: (binWidth - 2 * binWidth]
     * Block 2: (2 * binWidth - 3 * binWidth]
     * <p/>
     * Note that covered regions are open on the left ( and closed on the right ].
     * <p/>
     * In general, if block i is the ith block (starting from 0), then block i
     * contains all records that have starting position > (i * binWidth) and
     * <= ((i + 1) * binWidth))
     */
    public static class ChrIndex implements org.broad.tribble.index.ChrIndex {
        private String name = "";
        private int binWidth;
        private int longestFeature;
        private int nFeatures;
        private List<Block> blocks;

        private boolean OLD_V3_INDEX = false;

        /**
         * Default constructor needed for factory methods -- DO NOT REMOVE
         */
        public ChrIndex() {

        }

        ChrIndex(final String name, final int binWidth) {
            this.name = name;
            this.binWidth = binWidth;
            this.blocks = new ArrayList<Block>(100);
            this.longestFeature = 0;
            //this.largestBlockSize = 0;
            this.nFeatures = 0;
        }

        public String getName() {
            return name;
        }

        void addBlock(final Block block) {
            blocks.add(block);
            //largestBlockSize = Math.max(largestBlockSize, block.getSize());
        }

        public int getNBlocks() {
            return blocks.size();
        }

        public List<Block> getBlocks() {
            return blocks;
        }

        public List<Block> getBlocks(final int start, final int end) {
            if (blocks.isEmpty()) {
                return Collections.emptyList();
            } else {
                // Adjust position for the longest feature in this chromosome.  This insures we get
                // features that start before the bin but extend into it
                final int adjustedPosition = Math.max(start - longestFeature, 0);
                final int startBinNumber = adjustedPosition / binWidth;
                if (startBinNumber >= blocks.size()) // are we off the end of the bin list, so return nothing
                    return Collections.emptyList();
                else {
                    final int endBinNumber = Math.min((end - 1) / binWidth, blocks.size() - 1);

                    // By definition blocks are adjacent for the liner index.  Combine them into one merged block

                    final long startPos = blocks.get(startBinNumber).getStartPosition();
                    final long endPos = blocks.get(endBinNumber).getStartPosition() + blocks.get(endBinNumber).getSize();
                    final long size = endPos - startPos;
                    if (size == 0) {
                        return Collections.EMPTY_LIST;
                    } else {
                        final Block mergedBlock = new Block(startPos, size);
                        return Arrays.asList(mergedBlock);
                    }
                }
            }
        }


        public void updateLongestFeature(final int featureLength) {
            longestFeature = Math.max(longestFeature, featureLength);
        }

        public int getNFeatures() {
            return this.nFeatures;
        }

        public void incrementFeatureCount() {
            this.nFeatures++;
        }

        public void write(final LittleEndianOutputStream dos) throws IOException {

            // Chr name, binSize,  # bins,  longest feature
            dos.writeString(name);
            dos.writeInt(binWidth);
            dos.writeInt(blocks.size());
            dos.writeInt(longestFeature);
            dos.writeInt(0);    // no longer used
            //dos.writeInt(largestBlockSize);
            dos.writeInt(nFeatures);

            long pos = 0;
            long size = 0;
            for (final Block block : blocks) {
                pos = block.getStartPosition();
                size = block.getSize();
                dos.writeLong(pos);
            }
            // End of last block for this chromosome
            dos.writeLong(pos + size);
        }

        public void read(final LittleEndianInputStream dis) throws IOException {
            name = dis.readString();
            binWidth = dis.readInt();
            final int nBins = dis.readInt();
            longestFeature = dis.readInt();
            //largestBlockSize = dis.readInt();
            // largestBlockSize and totalBlockSize are old V3 index values.  largest block size should be 0 for
            // all newer V3 block.  This is a nasty hack that should be removed when we go to V4 (XML!) indices
            OLD_V3_INDEX = dis.readInt() > 0;
            nFeatures = dis.readInt();

            // note the code below accounts for > 60% of the total time to read an index
            blocks = new ArrayList<Block>(nBins);
            long pos = dis.readLong();
            for (int binNumber = 0; binNumber < nBins; binNumber++) {
                final long nextPos = dis.readLong();
                final long size = nextPos - pos;
                blocks.add(new Block(pos, size));
                pos = nextPos;
            }
        }

        public boolean equals(final Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof ChrIndex)) return false;
            final ChrIndex other = (ChrIndex) obj;
            return binWidth == other.binWidth
                    && longestFeature == other.longestFeature
                    //&& largestBlockSize == other.largestBlockSize
                    && nFeatures == other.nFeatures
                    && name.equals(other.name)
                    && blocks.equals(other.blocks);
        }

        /**
         * @return  Total size of all blocks
         */
        public long getTotalSize() {
            long n = 0;
            for (final Block b : getBlocks())
                n += b.getSize();
            return n;
        }

        public double getAverageFeatureSize() {
            return (1.0 * getTotalSize()) / getNFeatures();
        }

        public double getFeaturesPerBlock() {
            return (1.0 * getNFeatures()) / getNBlocks();
        }

        private double getNFeaturesOfMostDenseBlock(final double featureSize) {
            double m = -1;
            for (final Block b : getBlocks()) {
                final double n = b.getSize() / featureSize;
                if (m == -1 || n > m) m = n;
            }
            return m;
        }

        private double optimizeScore() {
            return getNFeaturesOfMostDenseBlock(getAverageFeatureSize());
        }

        public ChrIndex optimize(final double threshold) {
            return optimize(this, threshold, 0);
        }

        private static boolean badBinWidth(final ChrIndex idx) {
            if (idx.binWidth > MAX_BIN_WIDTH || idx.binWidth < 0) // an overflow occurred
                return true;
            else if (MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX != 0 && idx.getNFeatures() > 1 && idx.binWidth > MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX) {
                return true;
            } else {
                return false;
            }
        }

        private static ChrIndex optimize(ChrIndex idx, final double threshold, int level) {
            ChrIndex best = idx;

            while (true) {
                final double score = idx.optimizeScore();

                if (score > threshold || idx.getNBlocks() == 1 || badBinWidth(idx))
                    break;
                else {
                    best = idx; // remember the last best option

                    // try to make a better one
                    idx = mergeBlocks(idx);
                    level++;
                }

                if (level > 30) throw new IllegalStateException("Too many iterations");
            }

            return best;
        }

        private static ChrIndex mergeBlocks(final ChrIndex idx) {
            final ChrIndex merged = new ChrIndex(idx.name, idx.binWidth * 2); // increasing width by 2 each time
            merged.longestFeature = idx.longestFeature;
            merged.nFeatures = idx.nFeatures;

            final Iterator<Block> blocks = idx.getBlocks().iterator();
            if (!blocks.hasNext())
                throw new IllegalStateException("Block iterator cannot be empty at the start for " + idx.getName());

            // extremely simple merging algorithm.  Walk left to right, joining up blocks adjacent blocks.
            while (blocks.hasNext()) {
                final Block b1 = blocks.next();
                final Block b2 = blocks.hasNext() ? blocks.next() : null;

                if (b2 == null)
                    merged.addBlock(b1);
                else
                    // the new block is simply the start of the first block and the size of both together
                    merged.addBlock(new Block(b1.getStartPosition(), b1.getSize() + b2.getSize()));
            }

            return merged;
        }

        private static String dupString(final char c, final int nCopies) {
            final char[] chars = new char[nCopies];
            Arrays.fill(chars, c);
            return new String(chars);
        }
    }

    /**
     * Adapative optimization of the linear index
     * @param threshold threshold to use for optimizing each constituent {@code chrIndex}
     * @return The new optimized index
     */
    public Index optimize(final double threshold) {
        if (enableAdaptiveIndexing) {

            final List<ChrIndex> newIndices = new ArrayList<ChrIndex>(this.chrIndices.size());
            for (final String name : chrIndices.keySet()) {
                final LinearIndex.ChrIndex oldIdx = (LinearIndex.ChrIndex) chrIndices.get(name);
                final LinearIndex.ChrIndex newIdx = oldIdx.optimize(threshold);
                newIndices.add(newIdx);
            }
            return new LinearIndex(this, newIndices);
        } else {
            return this;
        }
    }

    public Index optimize() {
        return optimize(MAX_FEATURES_PER_BIN);
    }

    /**
     * Code to convert linear index to a text table for analysis
     * @param out Stream to which to write out table to
     */
    public void writeTable(final PrintStream out) {
        out.printf("chr binWidth avg.feature.size nFeatures.total block.id start.pos size nFeatures%n");
        for (final String name : chrIndices.keySet()) {
            final LinearIndex.ChrIndex chrIdx = (LinearIndex.ChrIndex) chrIndices.get(name);
            int blockCount = 0;
            for (final Block b : chrIdx.getBlocks()) {
                out.printf("%s %d %.2f %d %d %d %d %d%n", name, chrIdx.binWidth, chrIdx.getAverageFeatureSize(), chrIdx.getNFeatures(), blockCount,
                        blockCount * chrIdx.binWidth, b.getSize(), (int) (b.getSize() / chrIdx.getAverageFeatureSize()));
                blockCount++;
            }
        }
    }

    // purely for testing purposes
    protected final void setTS(final long ts) {
        this.indexedFileTS = ts;
    }
}

