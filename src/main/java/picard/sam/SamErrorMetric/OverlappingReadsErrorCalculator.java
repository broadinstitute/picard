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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.SamLocusAndReferenceIterator;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * A calculator that estimates the error rate of the bases it observes, assuming that the reference is truth.
 * This calculator only includes bases that have been read twice in the same template (once from each read) and
 * thus only includes bases that arise from the overlapping part of the reads. Over those bases the Calculator
 * distinguishes between whether the two reads agree with each other but differ from the reference (indicative
 * of a difference between the template and the reference, and when one of the reads agrees with the reference
 * but the other does not which indicates that there might have been a sequencing error in that read.
 */
public class OverlappingReadsErrorCalculator extends BaseErrorCalculator {
    private long nBothDisagreeWithReference;
    private long nDisagreeWithRefAndMate;
    private long nThreeWaysDisagreement;
    private long nTotalBasesWithOverlappingReads;

    private static int currentPosition;
    private static final Map<String, Set<SamLocusIterator.RecordAndOffset>> readNameSets = new CollectionUtil.DefaultingMap<>(s -> new HashSet<>(), true);

    private static void updateReadNameSet(final SamLocusIterator.LocusInfo locusInfo) {
        if (locusInfo.getPosition() == currentPosition) {
            return;
        }
        readNameSets.clear();
        locusInfo.getRecordAndOffsets().forEach(r -> readNameSets.get(r.getReadName()).add(r));
        currentPosition = locusInfo.getPosition();
    }

    /**
     * The function by which new loci are "shown" to the calculator
     **/
    @Override
    public void addBase(final SamLocusIterator.RecordAndOffset recordAndOffset, final SamLocusAndReferenceIterator.SAMLocusAndReference locusAndRef) {
        super.addBase(recordAndOffset, locusAndRef);
        final byte readBase = recordAndOffset.getReadBase();
        final SAMRecord record = recordAndOffset.getRecord();

        // by traversing the reads and splitting into sets with the same name we convert a O(N^2) iteration
        // into a O(N) iteration
        updateReadNameSet(locusAndRef.getLocus());

        final SamLocusIterator.RecordAndOffset mate = readNameSets.get(record.getReadName())
                .stream()
                .filter(putative -> areReadsMates(record, putative.getRecord()))
                .findFirst()
                .orElse(null);

        // we are only interested in bases for which the mate read also has a base over the same locus
        if (mate == null) return;
        // both bases need to be called for this error calculation

        final byte mateBase = mate.getReadBase();
        if (SequenceUtil.isNoCall(readBase)) return;
        if (SequenceUtil.isNoCall(mateBase)) return;

        nTotalBasesWithOverlappingReads++;

        // Only bases that disagree with the reference are counted as errors.
        if (SequenceUtil.basesEqual(readBase, locusAndRef.getReferenceBase())) return;

        final boolean agreesWithMate = SequenceUtil.basesEqual(readBase, mateBase);
        final boolean mateAgreesWithRef = SequenceUtil.basesEqual(mateBase, locusAndRef.getReferenceBase());

        if (agreesWithMate) {
            nBothDisagreeWithReference++;
        } else if (mateAgreesWithRef) {
            nDisagreeWithRefAndMate++;
        } else {
            nThreeWaysDisagreement++;
        }
    }

    /**
     * The suffix that pertains to the implementation of aggregation
     **/
    @Override
    public String getSuffix() {
        return "overlapping_error";
    }

    /**
     * Returns the metric generated by the observed loci
     **/
    @Override
    public OverlappingErrorMetric getMetric() {
        return new OverlappingErrorMetric("", nBases, nTotalBasesWithOverlappingReads, nDisagreeWithRefAndMate,
                nBothDisagreeWithReference, nThreeWaysDisagreement);
    }

    private boolean areReadsMates(final SAMRecord read1, final SAMRecord read2) {
        // must have same name
        return (read1.getReadName().equals(read2.getReadName()) &&
                // must be paired
                read1.getReadPairedFlag() &&
                // one must be first while the other is not
                read1.getFirstOfPairFlag() != read2.getFirstOfPairFlag() &&
                // one must be second while the other is not
                read1.getSecondOfPairFlag() != read2.getSecondOfPairFlag() &&
                // read1 must be mapped
                !read1.getReadUnmappedFlag() &&
                // read2 must be mapped
                !read2.getReadUnmappedFlag() &&
                // read1 must be non-secondary
                !read1.isSecondaryAlignment() &&
                // read2 must be non-secondary
                !read2.isSecondaryAlignment() &&
                // read1 mate position must agree with read2's position
                read1.getMateAlignmentStart() == read2.getAlignmentStart() &&
                // read1 mate reference must agree with read2's reference
                read1.getMateReferenceIndex().equals(read2.getReferenceIndex())
        );
    }
}
