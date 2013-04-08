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

package org.broad.tribble.index.interval;

import org.broad.tribble.index.AbstractIndex;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.IndexFactory;
import org.broad.tribble.util.LittleEndianInputStream;
import org.broad.tribble.util.LittleEndianOutputStream;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * Index based on an interval tree
 * @see IntervalTree
 * @author jrobinso
 * @date Jul 9, 2010
 */
public class IntervalTreeIndex extends AbstractIndex {

    /**
     * Default constructor -- used by factory methods.  Do not remove.
     */
    public IntervalTreeIndex() {

    }

    /**
     *
     * @param featureFile File which we are indexing
     */
    public IntervalTreeIndex(String featureFile) {
        super(featureFile);
    }

    @Override
    public Class getChrIndexClass() {
        return ChrIndex.class;
    }

    @Override
    protected int getType() {
        return IndexFactory.IndexType.INTERVAL_TREE.getHeaderValue();
    }

    /**
     * Add a new interval to this index
     * @param chr  Chromosome
     * @param interval
     */
    public void insert(String chr, Interval interval) {
        ChrIndex chrIdx = (ChrIndex) chrIndices.get(chr);
        if (chrIdx == null) {
            chrIdx = new ChrIndex(chr);
            chrIndices.put(chr, chrIdx);
        }
        chrIdx.insert(interval);
    }

    protected void setChrIndex(List<ChrIndex> indicies) {
        for (ChrIndex index : indicies) {
            chrIndices.put(index.getName(), index);
        }
    }

    public void printTree() {

        for (String chr : chrIndices.keySet()) {
            System.out.println(chr + ":");
            ChrIndex chrIdx = (ChrIndex) chrIndices.get(chr);
            chrIdx.printTree();
            System.out.println();
        }
    }

    public static class ChrIndex implements org.broad.tribble.index.ChrIndex {

        IntervalTree tree;
        String name;

        /**
         * Default constructor needed for factory methods -- DO NOT REMOVE
         */
        public ChrIndex() {

        }

        public ChrIndex(String name) {
            this.name = name;
            tree = new IntervalTree();
        }

        public String getName() {
            return name;
        }

        public void insert(Interval iv) {
            tree.insert(iv);
        }

        public List<Block> getBlocks() {
            return null;
        }


        public List<Block> getBlocks(int start, int end) {

            // Get intervals and build blocks list
            List<Interval> intervals = tree.findOverlapping(new Interval(start, end));

            // save time (and save throwing an exception) if the blocks are empty, return now
            if (intervals == null || intervals.size() == 0) return new ArrayList<Block>();

            Block[] blocks = new Block[intervals.size()];
            int idx = 0;
            for (Interval iv : intervals) {
                blocks[idx++] = iv.getBlock();
            }

            // Sort blocks by start position
            Arrays.sort(blocks, new Comparator<Block>() {
                public int compare(Block b1, Block b2) {
                    // this is a little cryptic because the normal method (b1.getStartPosition() - b2.getStartPosition()) wraps in int space and we incorrectly sort the blocks in extreme cases
                    return b1.getStartPosition() - b2.getStartPosition() < 1 ? -1 : (b1.getStartPosition() - b2.getStartPosition() > 1 ? 1 : 0);
                }
            });

            // Consolidate blocks  that are close together
            List<Block> consolidatedBlocks = new ArrayList(blocks.length);
            Block lastBlock = blocks[0];
            consolidatedBlocks.add(lastBlock);
            for (int i = 1; i < blocks.length; i++) {
                Block block = blocks[i];
                if (block.getStartPosition() < (lastBlock.getEndPosition() + 1000)) {
                    lastBlock.setEndPosition(block.getEndPosition());
                } else {
                    lastBlock = block;
                    consolidatedBlocks.add(lastBlock);
                }
            }

            return consolidatedBlocks;
        }

        public void printTree() {
            System.out.println(tree.toString());
        }

        public void write(LittleEndianOutputStream dos) throws IOException {

            dos.writeString(name);
            List<Interval> intervals = tree.getIntervals();

            dos.writeInt(intervals.size());
            for (Interval interval : intervals) {
                dos.writeInt(interval.start);
                dos.writeInt(interval.end);
                dos.writeLong(interval.getBlock().getStartPosition());
                dos.writeInt((int)interval.getBlock().getSize());
            }

        }

        public void read(LittleEndianInputStream dis) throws IOException {

            tree = new IntervalTree();

            name = dis.readString();
            int nIntervals = dis.readInt();
            while (nIntervals-- > 0) {

                int start = dis.readInt();
                int end = dis.readInt();
                long pos = dis.readLong();
                int size = dis.readInt();

                Interval iv = new Interval(start, end, new Block(pos, size));
                tree.insert(iv);
            }


        }

    }
}
