/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package org.broad.tribble.index.linear;

import org.broad.tribble.Feature;
import org.broad.tribble.index.Block;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexCreator;

import java.io.File;
import java.util.ArrayList;
import java.util.LinkedList;

/**
 * @author jrobinso
 */
public class LinearIndexCreator implements IndexCreator {
    public static int DEFAULT_BIN_WIDTH = 8000;
    // the set bin width
    private int binWidth = DEFAULT_BIN_WIDTH;

    // the input file
    private File inputFile;

    private LinkedList<LinearIndex.ChrIndex> chrList = new LinkedList<LinearIndex.ChrIndex>();
    private int longestFeature= 0;

    private ArrayList<Block> blocks = new ArrayList<Block>();

    public void initialize(File inputFile, int binSize) {
        this.inputFile = inputFile;
        binWidth = binSize;
    }

    /**
     * add a feature to the index
     * @param feature the feature, from which we use the contig, start, and stop
     * @param filePosition the position of the file at the BEGINNING of the current feature
     */
    public void addFeature(Feature feature, long filePosition) {
        // fi we don't have a chrIndex yet, or if the last one was for the previous contig, create a new one
        if (chrList.size() == 0 || !chrList.getLast().getName().equals(feature.getChr())) {
            // if we're creating a new chrIndex (not the first), make sure to dump the blocks to the old chrIndex
            if (chrList.size() != 0)
                for (int x = 0; x < blocks.size(); x++) {
                    blocks.get(x).setEndPosition((x + 1 == blocks.size()) ? filePosition : blocks.get(x + 1).getStartPosition());
                    chrList.getLast().addBlock(blocks.get(x));
                }
            chrList.add(new LinearIndex.ChrIndex(feature.getChr(),binWidth));
            blocks.clear();

            // Add the first block
            blocks.add(new Block(filePosition, 0));
            longestFeature = 0;
        }

        // if start > current bin location, make new bins until we're at the correct location
        while (feature.getStart() > blocks.size() * binWidth) {
            blocks.add(new Block(filePosition,0));
        }
        if ((feature.getEnd()- feature.getStart())+1 > longestFeature) {
            longestFeature = (feature.getEnd()- feature.getStart())+1;
            chrList.getLast().updateLongestFeature(longestFeature);
        }
        chrList.getLast().incrementFeatureCount();
    }

    /**
     * finalize the index; producing an index object
     * @param finalFilePosition the final file position, for indexes that have to close out with the final position
     * @return an Index object
     */
    public Index finalizeIndex(long finalFilePosition) {
        if (finalFilePosition == 0)
            throw new IllegalArgumentException("finalFilePosition != 0, -> " + finalFilePosition);

        for (int x = 0; x < blocks.size(); x++) {
            blocks.get(x).setEndPosition((x + 1 == blocks.size()) ? finalFilePosition : blocks.get(x+1).getStartPosition());
            chrList.getLast().addBlock(blocks.get(x));
        }
        blocks.clear();

        LinearIndex index = new LinearIndex(chrList,inputFile);
        index.finalizeIndex();
        return index.optimize();
    }

    /**
     * the current default bin size
     * @return
     */
    public int defaultBinSize() {
        return DEFAULT_BIN_WIDTH;
    }

    public int getBinSize() { return binWidth; }
}

