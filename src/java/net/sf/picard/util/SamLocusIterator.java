/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import net.sf.picard.PicardException;
import net.sf.picard.filter.*;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;

import java.util.*;

/**
 * Iterator that traverses a SAM File, accumulating information on a per-locus basis.
 * Optionally takes a target interval list, in which case the loci returned are the ones covered by
 * the interval list.  If no target interval list, whatever loci are covered by the input reads are returned.
 * By default duplicate reads and non-primary alignments are filtered out.  Filtering may be changed
 * via setSamFilters().
 *
 * @author alecw@broadinstitute.org
 */
public class SamLocusIterator implements Iterable<SamLocusIterator.LocusInfo>, CloseableIterator<SamLocusIterator.LocusInfo> {
    private static final Log LOG = Log.getInstance(SamLocusIterator.class);

    /**
     * A SAMRecord plus the zero-based offset in the read corresponding to the position in LocusInfo
     */
    public static class RecordAndOffset {
        private final SAMRecord record;
        private final int offset;

        public RecordAndOffset(final SAMRecord record, final int offset) {
            this.offset = offset;
            this.record = record;
        }

        /**
         * Zero-based offset into the read corresonding to the current position in LocusInfo
         */
        public int getOffset() {
            return offset;
        }

        public SAMRecord getRecord() {
            return record;
        }

        public byte getReadBase() {
            return record.getReadBases()[offset];
        }

        public byte getBaseQuality() {
            return record.getBaseQualities()[offset];
        }
    }

    /**
     * The unit of iteration.  Holds the locus, plus a ReadAndOffset for each read that overlaps the locus
     */
    public static class LocusInfo implements Locus {
        private final SAMSequenceRecord referenceSequence;
        private final int position;
        private final List<RecordAndOffset> recordAndOffsets = new ArrayList<RecordAndOffset>(100);

        LocusInfo(final SAMSequenceRecord referenceSequence, final int position) {
            this.referenceSequence = referenceSequence;
            this.position = position;
        }

        /**
         * Accumulate info for one read at the locus.
         */
        public void add(final SAMRecord read, final int position) {
            recordAndOffsets.add(new RecordAndOffset(read, position));
        }

        public int getSequenceIndex() { return referenceSequence.getSequenceIndex(); }

        /**
         * @return 1-based reference position
         */
        public int getPosition() { return position; }

        public List<RecordAndOffset> getRecordAndPositions() {
            return Collections.unmodifiableList(recordAndOffsets);
        }

        public String getSequenceName() { return referenceSequence.getSequenceName(); }
    }



    private final SAMFileReader samReader;
    private final ReferenceSequenceMask referenceSequenceMask;
    private PeekableIterator<SAMRecord> samIterator;
    private List<SamRecordFilter> samFilters = Arrays.asList(new NotPrimaryAlignmentFilter(),
                                                             new DuplicateReadFilter());
    private final List<Interval> intervals;
    private final boolean useIndex;

    // LocusInfos on this list are ready to be returned by iterator.  All reads that overlap
    // the locus have been accumulated before the LocusInfo is moved into this list.
    private final LinkedList<LocusInfo> complete = new LinkedList<LocusInfo>();

    // LocusInfos for which accumulation is in progress
    private final LinkedList<LocusInfo> accumulator = new LinkedList<LocusInfo>();

    private int qualityScoreCutoff = Integer.MIN_VALUE;

    private int mappingQualityScoreCutoff = Integer.MIN_VALUE;

    /**
     * If true, emit a LocusInfo for every locus in the target map, or if no target map,
     * emit a LocusInfo for every locus in the reference sequence.
     * If false, emit a LocusInfo only if a locus has coverage.
     */
    private boolean emitUncoveredLoci = true;

    // When there is a target mask, these members remember the last locus for which a LocusInfo has been
    // returned, so that any uncovered locus in the target mask can be covered by a 0-coverage LocusInfo
    private int lastReferenceSequence = 0;
    private int lastPosition = 0;

    // Set to true when past all aligned reads in input SAM file
    private boolean finishedAlignedReads = false;

    private final LocusComparator<Locus> locusComparator = new LocusComparator<Locus>();


    /**
     * Prepare to iterate through the given SAM records, skipping non-primary alignments.  Do not use
     * BAM index even if available.
     */
    public SamLocusIterator(final SAMFileReader samReader) {
        this(samReader, null);
    }

    /**
     * Prepare to iterate through the given SAM records, skipping non-primary alignments.  Do not use
     * BAM index even if available.
     *
     * @param intervalList Either the list of desired intervals, or null.  Note that if an intervalList is
     * passed in that is not coordinate sorted, it will eventually be coordinated sorted by this class.
     */
    public SamLocusIterator(final SAMFileReader samReader, final IntervalList intervalList) {
        this(samReader, intervalList, false);
    }

    /**
     * Prepare to iterate through the given SAM records, skipping non-primary alignments
     *
     * @param samReader must be coordinate sorted
     * @param intervalList Either the list of desired intervals, or null.  Note that if an intervalList is
     * passed in that is not coordinate sorted, it will eventually be coordinated sorted by this class.
     * @param useIndex If true, do indexed lookup to improve performance.  Not relevant if intervalList == null.
     * This can actually slow performance if the intervals are densely packed.
     */
    public SamLocusIterator(final SAMFileReader samReader, final IntervalList intervalList, final boolean useIndex) {
        if (samReader.getFileHeader().getSortOrder() == null || samReader.getFileHeader().getSortOrder() == SAMFileHeader.SortOrder.unsorted) {
            LOG.warn("SamLocusIterator constructed with samReader that has SortOrder == unsorted.  ", "" +
                    "Assuming SAM is coordinate sorted, but exceptions may occur if it is not.");
        } else if (samReader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("SamLocusIterator cannot operate on a SAM file that is not coordinate sorted.");
        }
        this.samReader = samReader;
        this.useIndex = useIndex;
        if (intervalList != null) {
            intervals = intervalList.getUniqueIntervals();
            this.referenceSequenceMask = new IntervalListReferenceSequenceMask(intervalList);
        } else {
            intervals = null;
            this.referenceSequenceMask = new WholeGenomeReferenceSequenceMask(samReader.getFileHeader());
        }
    }

    public Iterator<LocusInfo> iterator() {
        if (samIterator != null) {
            throw new IllegalStateException("Cannot call iterator() more than once on SamLocusIterator");
        }
        CloseableIterator<SAMRecord> tempIterator;
        if (intervals != null) {
            tempIterator = new SamRecordIntervalIteratorFactory().makeSamRecordIntervalIterator(samReader, intervals, useIndex);
        } else {
            tempIterator = samReader.iterator();
        }
        if (samFilters != null) {
            tempIterator = new FilteringIterator(tempIterator, new AggregateFilter(samFilters));
        }
        samIterator = new PeekableIterator<SAMRecord>(tempIterator);
        return this;
    }

    public void close() {
        this.samIterator.close();
    }

    private boolean samHasMore() {
        return !finishedAlignedReads && (samIterator.peek() != null);
    }

    /**
     * @return true if there are more aligned reads in the SAM file, LocusInfos in some stage of accumulation,
     * or loci in the target mask that have yet to be covered.
     */
    public boolean hasNext() {
        while (complete.isEmpty() && ((!accumulator.isEmpty()) || samHasMore() || hasRemainingMaskBases())) {
            final LocusInfo locusInfo = next();
            if (locusInfo != null) {
                complete.addFirst(locusInfo);
            }
        }
        return !complete.isEmpty();
    }

    /**
     * @return true if there are loci in the target mask that have yet to be covered by LocusInfos
     */
    private boolean hasRemainingMaskBases() {
        // if there are more sequences in the mask, by definition some of them must have
        // marked bases otherwise if we're in the last sequence, but we're not at the last marked position,
        // there is also more in the mask
        if (!emitUncoveredLoci) {
            // If not emitting uncovered loci, this check is irrelevant
            return false;
        }
        return (lastReferenceSequence < referenceSequenceMask.getMaxSequenceIndex() ||
               (lastReferenceSequence == referenceSequenceMask.getMaxSequenceIndex() &&
                lastPosition <= referenceSequenceMask.nextPosition(lastReferenceSequence, lastPosition+1)));
    }

    /**
     * hasNext() has been fixed so that if it returns true, next() is now guaranteed not to return null.
     */
    public LocusInfo next() {

        // if we don't have any completed entries to return, try and make some!
        while(complete.isEmpty() && samHasMore()) {
            final SAMRecord rec = samIterator.peek();

            // There might be unmapped reads mixed in with the mapped ones, but when a read
            // is encountered with no reference index it means that all the mapped reads have been seen.
            if (rec.getReferenceIndex() == -1) {
                this.finishedAlignedReads = true;
                continue;

            }
            // Skip over an unaligned read that has been forced to be sorted with the aligned reads
            if (rec.getReadUnmappedFlag()) {
                samIterator.next();
                continue;
            }

            final Locus alignmentStart = new LocusImpl(rec.getReferenceIndex(), rec.getAlignmentStart());

            // emit everything that is before the start of the current read, because we know no more
            // coverage will be accumulated for those loci.
            while (!accumulator.isEmpty() && locusComparator.compare(accumulator.getFirst(), alignmentStart) < 0) {
                final LocusInfo first = accumulator.getFirst();
                populateCompleteQueue(alignmentStart);
                if (!complete.isEmpty()) {
                    return complete.removeFirst();
                }
                if (!accumulator.isEmpty() && first == accumulator.getFirst()) {
                    throw new PicardException("Stuck in infinite loop");
                }
            }

            // at this point, either the accumulator list is empty or the head should
            // be the same position as the first base of the read
            if (!accumulator.isEmpty()) {
                if (accumulator.getFirst().getSequenceIndex() != rec.getReferenceIndex() ||
                        accumulator.getFirst().position != rec.getAlignmentStart()) {
                    throw new IllegalStateException("accumulator should be empty or aligned with current SAMRecord");
                }
            }

            // Store the loci for the read in the accumulator
            accumulateSamRecord(rec);

            samIterator.next();
        }

        final Locus endLocus = new LocusImpl(Integer.MAX_VALUE, Integer.MAX_VALUE);
        // if we have nothing to return to the user, and we're at the end of the SAM iterator,
        // push everything into the complete queue
        if (complete.isEmpty() && !samHasMore()) {
            while(!accumulator.isEmpty()) {
                populateCompleteQueue(endLocus);
                if (!complete.isEmpty()) {
                    return complete.removeFirst();
                }
            }
        }

        // if there are completed entries, return those
        if (!complete.isEmpty()) {
            return complete.removeFirst();
        } else if (emitUncoveredLoci){
            final Locus afterLastMaskPositionLocus = new LocusImpl(referenceSequenceMask.getMaxSequenceIndex(),
                    referenceSequenceMask.getMaxPosition() + 1);
            // In this case... we're past the last read from SAM so see if we can
            // fill out any more (zero coverage) entries from the mask
            return createNextUncoveredLocusInfo(afterLastMaskPositionLocus);
        } else {
            return null;
        }
    }

    /**
     * Capture the loci covered by the given SAMRecord in the LocusInfos in the accumulator,
     * creating new LocusInfos as needed.
     */
    private void accumulateSamRecord(final SAMRecord rec) {
        // interpret the CIGAR string and add the base info
        for(final AlignmentBlock alignmentBlock : rec.getAlignmentBlocks()) {
            for (int i = 0; i < alignmentBlock.getLength(); ++i) {
                // 0-based offset into the read of the current base
                final int readOffset = alignmentBlock.getReadStart() + i - 1;
                // 1-based reference position that the current base aligns to
                final int refPos = alignmentBlock.getReferenceStart() + i;

                // 0-based offset from the aligned position of the first base in the read to the aligned position
                // of the current base.
                final int refOffset =  refPos - rec.getAlignmentStart();

                // Ensure there are LocusInfos up to and including this position
                for (int j = accumulator.size(); j <= refOffset; ++j) {
                    accumulator.add(new LocusInfo(getReferenceSequence(rec.getReferenceIndex()),
                            rec.getAlignmentStart() + j));
                }
                // if the quality score cutoff is met, accumulate the base info
                if (rec.getBaseQualities()[readOffset] >= getQualityScoreCutoff() &&
                        rec.getMappingQuality() >= getMappingQualityScoreCutoff()) {
                    accumulator.get(refOffset).add(rec, readOffset);
                }
            }
        }
    }

    /**
     * Create the next relevant zero-coverage LocusInfo
     * @param stopBeforeLocus don't go up to this sequence and position
     * @return a zero-coverage LocusInfo, or null if there is none before the stopBefore locus
     */
    private LocusInfo createNextUncoveredLocusInfo(final Locus stopBeforeLocus) {
        while (lastReferenceSequence <= stopBeforeLocus.getSequenceIndex() &&
               lastReferenceSequence <= referenceSequenceMask.getMaxSequenceIndex()) {

            if (lastReferenceSequence == stopBeforeLocus.getSequenceIndex() &&
                lastPosition +1 >= stopBeforeLocus.getPosition()) {
                return null;
            }



            final int nextbit = referenceSequenceMask.nextPosition(lastReferenceSequence, lastPosition+1);

            // try the next reference sequence
            if (nextbit == -1) {
                // No more in this reference sequence
                if (lastReferenceSequence == stopBeforeLocus.getSequenceIndex()) {
                    lastPosition = stopBeforeLocus.getPosition();
                    return null;
                }
                lastReferenceSequence++;
                lastPosition = 0;
            } else if (lastReferenceSequence < stopBeforeLocus.getSequenceIndex() ||
                    nextbit < stopBeforeLocus.getPosition()) {
                lastPosition = nextbit;
                return new LocusInfo(getReferenceSequence(lastReferenceSequence), lastPosition);
            } else if (nextbit >= stopBeforeLocus.getPosition()) {
                return null;
            }

        }

        return null;
    }

    /**
     * Pop the first entry from the LocusInfo accumulator into the complete queue.  In addition,
     * check the ReferenceSequenceMask and if there are intervening mask positions between the last popped base and the one
     * about to be popped, put those on the complete queue as well.
     * Note that a single call to this method may not empty the accumulator completely, or even
     * empty it at all, because it may just put a zero-coverage LocusInfo into the complete queue.
     */
    private void populateCompleteQueue(final Locus stopBeforeLocus) {
        // Because of gapped alignments, it is possible to create LocusInfo's with no reads associated with them.
        // Skip over these.
        while (!accumulator.isEmpty() && accumulator.getFirst().getRecordAndPositions().isEmpty() &&
               locusComparator.compare(accumulator.getFirst(), stopBeforeLocus) < 0) {
            accumulator.removeFirst();
        }
        if (accumulator.isEmpty()) {
            return;
        }
        final LocusInfo locusInfo = accumulator.getFirst();
        if (locusComparator.compare(stopBeforeLocus, locusInfo) <= 0) {
            return;
        }

        // If necessary, emit a zero-coverage LocusInfo
        if (emitUncoveredLoci) {
            final LocusInfo zeroCoverage = createNextUncoveredLocusInfo(locusInfo);
            if (zeroCoverage != null) {
                complete.addLast(zeroCoverage);
                return;
            }
        }

        // At this point we know we're going to process the LocusInfo, so remove it from the accumulator.
        accumulator.removeFirst();

        // fill in any gaps based on our genome mask
        final int sequenceIndex = locusInfo.getSequenceIndex();


        // only add to the complete queue if it's in the mask (or we have no mask!)
        if (referenceSequenceMask.get(locusInfo.getSequenceIndex(), locusInfo.getPosition())) {
            complete.addLast(locusInfo);
        }

        lastReferenceSequence = sequenceIndex;
        lastPosition = locusInfo.getPosition();
    }

    private SAMSequenceRecord getReferenceSequence(final int referenceSequenceIndex) {
        return samReader.getFileHeader().getSequence(referenceSequenceIndex);
    }

    public void remove() {
        throw new UnsupportedOperationException("Can not remove records from a SAM file via an iterator!");
    }

    // --------------------------------------------------------------------------------------------
    // Helper methods below this point...
    // --------------------------------------------------------------------------------------------

    /**
     * Controls which, if any, SAMRecords are filtered.  By default duplicate reads and non-primary alignments
     * are filtered out.  The list of filters passed here replaces any existing filters.
     * @param samFilters list of filters, or null if no filtering is desired.
     */
    public void setSamFilters(final List<SamRecordFilter> samFilters) {
        this.samFilters = samFilters;
    }

    public int getQualityScoreCutoff() { return qualityScoreCutoff; }
    public void setQualityScoreCutoff(final int qualityScoreCutoff) { this.qualityScoreCutoff = qualityScoreCutoff; }

    public int getMappingQualityScoreCutoff() {
        return mappingQualityScoreCutoff;
    }

    public void setMappingQualityScoreCutoff(final int mappingQualityScoreCutoff) {
        this.mappingQualityScoreCutoff = mappingQualityScoreCutoff;
    }

    public boolean isEmitUncoveredLoci() {
        return emitUncoveredLoci;
    }

    public void setEmitUncoveredLoci(final boolean emitUncoveredLoci) {
        this.emitUncoveredLoci = emitUncoveredLoci;
    }
}
