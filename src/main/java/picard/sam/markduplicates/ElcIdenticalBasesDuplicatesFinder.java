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

import java.util.ArrayList;
import java.util.List;

import static picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;
import static picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequenceWithBarcodes;

/**
 * Algorithm for search duplicates is used in EstimateLibraryComplexity only.
 * It finds duplicates comparing bases in pairs.
 *
 * * @author Pavel_Silin@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
class ElcIdenticalBasesDuplicatesFinder extends ElcDuplicatesFinder {

    private boolean useBarcodes;

    ElcIdenticalBasesDuplicatesFinder(double maxDiffRate, int maxReadLength, int minIdenticalBases,
                                      boolean useBarcodes, OpticalDuplicateFinder opticalDuplicateFinder) {
        super(maxDiffRate, maxReadLength, minIdenticalBases, opticalDuplicateFinder);
        this.useBarcodes = useBarcodes;
    }

    /**
     * Search duplicates in library size < BOUNDARY_LIBRARY_SIZE or reads with barcodes
     */
    @Override
    void searchDuplicates(List<PairedReadSequence> sequences,
                          Histogram<Integer> duplicationHisto, Histogram<Integer> opticalHisto) {
        for (int i = 0; i < sequences.size(); ++i) {
            final PairedReadSequence lhs = sequences.get(i);
            if (lhs == null) continue;
            final List<PairedReadSequence> dupes = new ArrayList<>();

            for (int j = i + 1; j < sequences.size(); ++j) {
                final PairedReadSequence rhs = sequences.get(j);
                if (rhs == null) continue;

                if (matches(lhs, rhs, maxDiffRate, useBarcodes)) {
                    dupes.add(rhs);
                    sequences.set(j, null);
                }
            }

            fillHistogram(duplicationHisto, opticalHisto, lhs, dupes);

        }
    }

    /**
     * Checks to see if two reads pairs have sequence that are the same, give or take a few
     * errors/diffs as dictated by the maxDiffRate.
     */
    private boolean matches(final PairedReadSequence lhs, final PairedReadSequence rhs, final double maxDiffRate, final boolean useBarcodes) {
        final int read1Length = minLength(lhs.read1, rhs.read1);
        final int read2Length = minLength(lhs.read2, rhs.read2);
        final int maxErrors = (int) Math.floor((read1Length + read2Length) * maxDiffRate);
        int errors = 0;

        if (useBarcodes) {
            final PairedReadSequenceWithBarcodes lhsWithBarcodes = (PairedReadSequenceWithBarcodes) lhs;
            final PairedReadSequenceWithBarcodes rhsWithBarcodes = (PairedReadSequenceWithBarcodes) rhs;
            if (lhsWithBarcodes.barcode != rhsWithBarcodes.barcode ||
                    lhsWithBarcodes.readOneBarcode != rhsWithBarcodes.readOneBarcode ||
                    lhsWithBarcodes.readTwoBarcode != rhsWithBarcodes.readTwoBarcode) {
                return false;
            }
        }

        // The loop can start from MIN_IDENTICAL_BASES because we've already confirmed that
        // at least those first few bases are identical when sorting.
        for (int i = minIdenticalBases; i < read1Length; ++i) {
            if (lhs.read1[i] != rhs.read1[i] && ++errors > maxErrors) {
                return false;
            }
        }

        for (int i = minIdenticalBases; i < read2Length; ++i) {
            if (lhs.read2[i] != rhs.read2[i] && ++errors > maxErrors) {
                return false;
            }
        }

        return true;
    }


}
