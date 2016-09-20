/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

package picard.sam.markduplicates;

import htsjdk.samtools.util.Histogram;
import picard.sam.markduplicates.util.OpticalDuplicateFinder;

import java.util.List;

import static picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;

/**
 * Abstract class for search duplicates algorithm which is used in the EstimateLibraryComplexity
 *
 * @author Pavel_Silin@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
abstract class ElcDuplicatesFinder {

    protected final int maxReadLength;
    protected final double maxDiffRate;
    protected int minIdenticalBases;
    protected OpticalDuplicateFinder opticalDuplicateFinder;

    ElcDuplicatesFinder(double maxDiffRate, int maxReadLength, int minIdenticalBases,
                        OpticalDuplicateFinder opticalDuplicateFinder) {
        this.maxDiffRate = maxDiffRate;
        this.minIdenticalBases = minIdenticalBases;
        this.opticalDuplicateFinder = opticalDuplicateFinder;
        this.maxReadLength = (maxReadLength <= 0) ? Integer.MAX_VALUE : maxReadLength;
    }

    abstract void searchDuplicates(List<PairedReadSequence> sequences,
                                   Histogram<Integer> duplicationHisto, Histogram<Integer> opticalHisto);

    /**
     * Fill histograms based on duplicates.
     */
    protected void fillHistogram(Histogram<Integer> duplicationHisto, Histogram<Integer> opticalHisto,
                                 PairedReadSequence prs, List<PairedReadSequence> dupes) {
        if (!dupes.isEmpty()) {
            dupes.add(prs);
            final int duplicateCount = dupes.size();
            duplicationHisto.increment(duplicateCount);

            final boolean[] flags = opticalDuplicateFinder.findOpticalDuplicates(dupes, prs);
            for (final boolean b : flags) {
                if (b) opticalHisto.increment(duplicateCount);
            }
        } else {
            duplicationHisto.increment(1);
        }
    }

    protected int minLength(byte[] read1, byte[] read2) {
        return Math.min(Math.min(read1.length, read2.length), maxReadLength);
    }
}