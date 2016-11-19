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
 * Resolve algorithm is used for the current PairedReadSequence group
 *
 * * @author Pavel_Silin@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
class ElcDuplicatesFinderResolver {

    /**
     * This parameter determines the choice of the algorithm: if the group size > BOUNDARY_LIBRARY_SIZE, the modified
     * algorithm applies
     */
    private static final int BOUNDARY_LIBRARY_SIZE = 100;

    private boolean useBarcodes;
    private ElcHashBasedDuplicatesFinder hashBasedDuplicatesFinder;
    private ElcIdenticalBasesDuplicatesFinder identicalBasesDuplicateFinder;

    ElcDuplicatesFinderResolver(double maxDiffRate, int maxReadLength, int minIdenticalBases, boolean useBarcodes,
                                OpticalDuplicateFinder opticalDuplicateFinder) {
        this.useBarcodes = useBarcodes;

        this.hashBasedDuplicatesFinder = new ElcHashBasedDuplicatesFinder(
                maxDiffRate,
                maxReadLength,
                minIdenticalBases,
                opticalDuplicateFinder
        );

        this.identicalBasesDuplicateFinder = new ElcIdenticalBasesDuplicatesFinder(
                maxDiffRate,
                maxReadLength,
                minIdenticalBases,
                useBarcodes,
                opticalDuplicateFinder
        );
    }

    void resolveAndSearch(List<PairedReadSequence> sequences, Histogram<Integer> duplicationHisto,
                          Histogram<Integer> opticalHisto) {
        if (useBarcodes || sequences.size() < BOUNDARY_LIBRARY_SIZE) {
            identicalBasesDuplicateFinder.searchDuplicates(sequences, duplicationHisto, opticalHisto);
        } else {
            hashBasedDuplicatesFinder.searchDuplicates(sequences, duplicationHisto, opticalHisto);
        }
    }

}
