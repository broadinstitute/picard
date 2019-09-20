/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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

package picard.sam.SamErrorMetric;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.List;
import java.util.Map;

public abstract class BaseErrorCalculator implements BaseCalculator {
    long nBases;

    // TODO mgatzen debug
    /**
     * Number of bases which are skipped because they overlap with a SNP variant site
     */
    protected long nSkippedSNPs = 0;

    /**
     * Number of insertions or deletions which are skipped because they overlap with an indel variant site. Note that
     * this is not the number of bases that are skipped, i.e. each insertion or deletion is counted only once.
     */
    protected long nSkippedIndels = 0;

    /**
     * the function by which new loci are "shown" to the calculator
     **/
    @Override
    public void addBase(final SamLocusIterator.RecordAndOffset recordAndOffset, final SamLocusAndReferenceIterator.SAMLocusAndReference locusAndRef, final Map<Integer, List<VariantContext>> potentialVariants) {

        if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Match) {
            // Search for SNP variants
            if (potentialVariants.containsKey(0)) {
                for (VariantContext variantContext : potentialVariants.get(0)) {
                    if (variantContext.getType() == VariantContext.Type.SNP) {
                        // Don't consider mismatches at a SNP variant site as an error
                        nSkippedSNPs++;
                        return;
                    }
                }
            }
        }
        else {
            // Search for indel variants
            for (final Map.Entry<Integer, List<VariantContext>> entry : potentialVariants.entrySet()) {
                for (final VariantContext variantContext : entry.getValue()) {
                    if (variantContext.getType() == VariantContext.Type.INDEL) {
                        // Don't consider records at an indel variant site (or its surrounding loci) sas an error.
                        nSkippedIndels++;
                        return;
                    }
                }
            }
        }

        // TODO mgatzen should deletions be counted towards total bases?
        if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Match) {
            if (!SequenceUtil.isNoCall(recordAndOffset.getReadBase())) {
                nBases++;
            }
        } else if (recordAndOffset.getAlignmentType() == AbstractRecordAndOffset.AlignmentType.Insertion) {
            CigarElement cigarElement = ReadBaseStratification.getIndelElement(recordAndOffset);
            if (cigarElement != null) {
                nBases += cigarElement.getLength();
            }
        }
    }
}
