/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package net.sf.picard.sam;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * This strategy was designed for TopHat output, but could be of general utility.  It picks the alignment with best MAPQ.
 * If paired-end, it is the alignment in which the sum of the MAPQs of both ends is the best.  In case of ties, one
 * is selected arbitrarily.
 */
public class BestMapqPrimaryAlignmentSelectionStrategy implements PrimaryAlignmentSelectionStrategy {
    // Give random number generator a seed so results are repeatable.  Used to pick a primary alignment from
    // multiple alignments with equal mapping quality.
    private final Random random = new Random(1);

    /**
     * Primary alignment was filtered out.  Need to select a new one.
     */
    public void pickPrimaryAlignment(final MultiHitAlignedReadIterator.HitsForInsert hits) {

        if (hits.numHits() == 0) throw new IllegalArgumentException("No alignments to pick from");

        // Find all the hits with the best MAPQ.
        final List<Integer> primaryAlignmentIndices = new ArrayList<Integer>(hits.numHits());
        int bestMapQ = -1;
        for (int i = 0; i < hits.numHits(); ++i) {
            int thisMapQ = 0;
            if (hits.getFirstOfPair(i) != null) {
                thisMapQ += hits.getFirstOfPair(i).getMappingQuality();
            }
            if (hits.getSecondOfPair(i) != null) {
                thisMapQ += hits.getSecondOfPair(i).getMappingQuality();
            }
            if (thisMapQ > bestMapQ) {
                bestMapQ = thisMapQ;
                primaryAlignmentIndices.clear();
            }
            if (thisMapQ == bestMapQ) primaryAlignmentIndices.add(i);
        }

        // Of all the hits with the best MAPQ, randomly select one to be primary.
        final int primaryAlignmentIndex;
        if (primaryAlignmentIndices.size() == 1) primaryAlignmentIndex = primaryAlignmentIndices.get(0);
        else if (primaryAlignmentIndices.size() > 1) primaryAlignmentIndex =
                primaryAlignmentIndices.get(random.nextInt(primaryAlignmentIndices.size()));
        else throw new IllegalStateException("Never found a best MAPQ -- should never happen");

        hits.setPrimaryAlignment(primaryAlignmentIndex);
    }
}
