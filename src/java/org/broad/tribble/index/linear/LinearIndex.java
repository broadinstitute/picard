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

import org.broad.tribble.index.AbstractIndex;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.File;
import java.io.IOException;
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

public class LinearIndex extends AbstractIndex implements Index {

    // NOTE: To debug uncomment the System.getProperty and recompile.
    public static final double MAX_FEATURES_PER_BIN = Double.valueOf(System.getProperty("MAX_FEATURES_PER_BIN", "100"));

    private final static int MAX_BIN_WIDTH = 1 * 1000 * 1000 * 1000; //  widths must be less than 1 billion

    // 1MB: we will no merge bins with any features in them beyond this size, no matter how sparse, per chromosome
    private static final long MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX = Long.valueOf(System.getProperty("MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX", "1024000"));

    public static boolean enableAdaptiveIndexing = true;

    /**
     * Default constructor -- used by factory methods.  Do not remove.
     */
    public LinearIndex() {}

    public LinearIndex(List<ChrIndex> indices, File inputFile) {
        super(inputFile.getAbsolutePath());
        for (ChrIndex index : indices)
            chrIndices.put(index.getName(), index);
    }

    private LinearIndex(LinearIndex parent, List<ChrIndex> indices) {
        super(parent);
        for (ChrIndex index : indices)
            chrIndices.put(index.getName(), index);
    }

    public LinearIndex(String featureFile) {
        super(featureFile);
    }


    public boolean isCurrentVersion() {
        if (!super.isCurrentVersion()) return false;

        // todo fixme nasty hack to determine if this is an old style V3 linear index (without nFeaturesPerBin)
        for (org.broad.tribble.index.ChrIndex chrIndex : chrIndices.values())
            if (((ChrIndex) chrIndex).OLD_V3_INDEX)
                return false;

        return true;
    }

    @Override
    protected int getType() {
        return IndexFactory.IndexType.LINEAR.getHeaderValue();
    }

    public LinkedHashSet<String> getSequenceNames() {
        return (chrIndices == null ? new LinkedHashSet<String>() : new LinkedHashSet<String>(chrIndices.keySet()));
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

        ChrIndex(String name, int binWidth) {
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

        void addBlock(Block block) {
            blocks.add(block);
            //largestBlockSize = Math.max(largestBlockSize, block.getSize());
        }

        public int getNBlocks() {
            return blocks.size();
        }

        public List<Block> getBlocks() {
            return blocks;
        }

        /**
         * @param start the start position, one based
         * @param end   the end position, one based
         * @return a list of blocks that include the region defined from start to stop.  Can never return null
         */
        public List<Block> getBlocks(int start, int end) {
            if (blocks.isEmpty()) {
                return Collections.emptyList();
            } else {
                // Adjust position for the longest feature in this chromosome.  This insures we get
                // features that start before the bin but extend into it
                int adjustedPosition = Math.max(start - longestFeature, 0);
                int startBinNumber = adjustedPosition / binWidth;
                if (startBinNumber >= blocks.size()) // are we off the end of the bin list, so return nothing
                    return Collections.emptyList();
                else {
                    int endBinNumber = Math.min((end - 1) / binWidth, blocks.size() - 1);

                    // By definition blocks are adjacent for the liner index.  Combine them into one merged block

                    long startPos = blocks.get(startBinNumber).getStartPosition();
                    long endPos = blocks.get(endBinNumber).getStartPosition() + blocks.get(endBinNumber).getSize();
                    long size = endPos - startPos;
                    if (size == 0) {
                        return Collections.EMPTY_LIST;
                    } else {
                        Block mergedBlock = new Block(startPos, size);
                        return Arrays.asList(mergedBlock);
                    }
                }
            }
        }


        public void updateLongestFeature(int featureLength) {
            longestFeature = Math.max(longestFeature, featureLength);
        }

        public int getNFeatures() {
            return this.nFeatures;
        }

        public void incrementFeatureCount() {
            this.nFeatures++;
        }

        public void write(LittleEndianOutputStream dos) throws IOException {

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
            for (Block block : blocks) {
                pos = block.getStartPosition();
                size = block.getSize();
                dos.writeLong(pos);
            }
            // End of last block for this chromosome
            dos.writeLong(pos + size);
        }

        public void read(LittleEndianInputStream dis) throws IOException {
            name = dis.readString();
            binWidth = dis.readInt();
            int nBins = dis.readInt();
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
                long nextPos = dis.readLong();
                long size = nextPos - pos;
                blocks.add(new Block(pos, size));
                pos = nextPos;
            }
        }

        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof ChrIndex)) return false;
            ChrIndex other = (ChrIndex) obj;
            return binWidth == other.binWidth
                    && longestFeature == other.longestFeature
                    //&& largestBlockSize == other.largestBlockSize
                    && nFeatures == other.nFeatures
                    && name.equals(other.name)
                    && blocks.equals(other.blocks);
        }

        public long getTotalSize() {
            long n = 0;
            for (Block b : getBlocks())
                n += b.getSize();
            return n;
        }

        public double getAverageFeatureSize() {
            return (1.0 * getTotalSize()) / getNFeatures();
        }

        public double getFeaturesPerBlock() {
            return (1.0 * getNFeatures()) / getNBlocks();
        }

        private double getNFeaturesOfMostDenseBlock(double featureSize) {
            double m = -1;
            for (Block b : getBlocks()) {
                double n = b.getSize() / featureSize;
                if (m == -1 || n > m) m = n;
            }
            return m;
        }

        private double optimizeScore() {
            return getNFeaturesOfMostDenseBlock(getAverageFeatureSize());
        }

        public ChrIndex optimize(double threshold) {
            return optimize(this, threshold, 0);
        }

        private static boolean badBinWidth(ChrIndex idx) {
            if (idx.binWidth > MAX_BIN_WIDTH || idx.binWidth < 0) // an overflow occurred
                return true;
            else if (MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX != 0 && idx.getNFeatures() > 1 && idx.binWidth > MAX_BIN_WIDTH_FOR_OCCUPIED_CHR_INDEX) {
                return true;
            } else {
                return false;
            }
        }

        private static ChrIndex optimize(ChrIndex idx, double threshold, int level) {
            ChrIndex best = idx;

            while (true) {
                double score = idx.optimizeScore();

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

        private static ChrIndex mergeBlocks(ChrIndex idx) {
            ChrIndex merged = new ChrIndex(idx.name, idx.binWidth * 2); // increasing width by 2 each time
            merged.longestFeature = idx.longestFeature;
            merged.nFeatures = idx.nFeatures;

            Iterator<Block> blocks = idx.getBlocks().iterator();
            if (!blocks.hasNext())
                throw new IllegalStateException("Block iterator cannot be empty at the start for " + idx.getName());

            // extremely simple merging algorithm.  Walk left to right, joining up blocks adjacent blocks.
            while (blocks.hasNext()) {
                Block b1 = blocks.next();
                Block b2 = blocks.hasNext() ? blocks.next() : null;

                if (b2 == null)
                    merged.addBlock(b1);
                else
                    // the new block is simply the start of the first block and the size of both together
                    merged.addBlock(new Block(b1.getStartPosition(), b1.getSize() + b2.getSize()));
            }

            return merged;
        }

        private static String dupString(char c, int nCopies) {
            char[] chars = new char[nCopies];
            Arrays.fill(chars, c);
            return new String(chars);
        }
    }

    // ----------------------------------------------------------------------------------------------------
    //
    // Adapative optimization of the linear index
    //
    // ----------------------------------------------------------------------------------------------------
    public Index optimize(double threshold) {
        if (enableAdaptiveIndexing) {

            List<ChrIndex> newIndices = new ArrayList<ChrIndex>(this.chrIndices.size());
            for (String name : chrIndices.keySet()) {
                LinearIndex.ChrIndex oldIdx = (LinearIndex.ChrIndex) chrIndices.get(name);
                LinearIndex.ChrIndex newIdx = oldIdx.optimize(threshold);
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


    // ----------------------------------------------------------------------------------------------------
    //
    // Code to convert linear index to a text table for analysis
    //
    // ----------------------------------------------------------------------------------------------------
    public void writeTable(PrintStream out) {
        out.printf("chr binWidth avg.feature.size nFeatures.total block.id start.pos size nFeatures%n");
        for (String name : chrIndices.keySet()) {
            LinearIndex.ChrIndex chrIdx = (LinearIndex.ChrIndex) chrIndices.get(name);
            int blockCount = 0;
            for (Block b : chrIdx.getBlocks()) {
                out.printf("%s %d %.2f %d %d %d %d %d%n", name, chrIdx.binWidth, chrIdx.getAverageFeatureSize(), chrIdx.getNFeatures(), blockCount,
                        blockCount * chrIdx.binWidth, b.getSize(), (int) (b.getSize() / chrIdx.getAverageFeatureSize()));
                blockCount++;
            }
        }
    }

    // purely for testing purposes
    protected final void setTS(long ts) {
        this.indexedFileTS = ts;
    }
}

