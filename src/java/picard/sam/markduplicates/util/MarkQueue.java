/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

package picard.sam.markduplicates.util;

import htsjdk.samtools.DuplicateScoringStrategy;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SamRecordTrackingBuffer;
import picard.PicardException;
import picard.sam.DuplicationMetrics;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;

import java.util.Comparator;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * This is the mark queue.
 *
 * This stores a current set of read ends that need to be duplicate marked.  It only stores internally the "best" read end for a given
 * possible duplicate location, preferring to perform duplicate marking as read ends come in, rather than wait for all "comparable"
 * read ends to arrive.  This reduces the memory footprint of this data structure.
 */
public class MarkQueue {

    /** Comparator to order the mark queue set.  The set of all the read ends that have are compared to be the same should
     * be used for duplicate marking. */
    private class MarkQueueComparator implements Comparator<ReadEndsForMateCigar> {
        public int compare(final ReadEndsForMateCigar lhs, final ReadEndsForMateCigar rhs) {
            int retval = lhs.libraryId - rhs.libraryId;
            if (retval == 0) retval = lhs.read1ReferenceIndex - rhs.read1ReferenceIndex;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = rhs.orientation - lhs.orientation; // to get pairs first
            if (retval == 0) retval = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            return retval;
        }
    }

    /** Comparator for ReadEndsForMateCigar that orders by read1 position then pair orientation then read2 position.
     * This is the main comparator to choose the best representative read and end to store as the non-duplicate.
     */
    class ReadEndsMCComparator implements Comparator<ReadEndsForMateCigar> {
        private final ScoringStrategy duplicateScoringStrategy;

        public ReadEndsMCComparator(final ScoringStrategy duplicateScoringStrategy) {
            this.duplicateScoringStrategy = duplicateScoringStrategy;
        }

        public int compare(final ReadEndsForMateCigar lhs, final ReadEndsForMateCigar rhs) {
            int retval = lhs.libraryId - rhs.libraryId;
            if (retval == 0) retval = lhs.read1ReferenceIndex - rhs.read1ReferenceIndex;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = rhs.orientation - lhs.orientation; // IMPORTANT: reverse the order to get pairs first
            if (retval == 0 && lhs.isPaired() != rhs.isPaired()) return lhs.isPaired() ? -1 : 1; // unpaired goes first...
            if (retval == 0) retval = lhs.hasUnmapped - rhs.hasUnmapped;
            if (retval == 0) retval = lhs.read2ReferenceIndex - rhs.read2ReferenceIndex;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            // TODO: cache the scores?
            if (retval == 0) retval = DuplicateScoringStrategy.compare(lhs.getRecord(), rhs.getRecord(), this.duplicateScoringStrategy, true);
            if (retval == 0) retval = lhs.getRecordReadName().compareTo(rhs.getRecordReadName());
            /**
             * If both reads are paired and both ends mapped, always prefer the first end over the second end.  This is needed to
             * properly choose the first end for optical duplicate identification when both ends are mapped to the same position etc.
             */
            if (retval == 0 && lhs.isPaired() && rhs.isPaired() && null != lhs.getsamRecordIndex()) {
                if (lhs.getRecord().getFirstOfPairFlag()) retval = -1;
                else retval = 1;
            }

            return retval;
        }
    }

    /** The genomic distance needed to assure that we have considered all reads for duplicate marking. */
    private int toMarkQueueMinimumDistance = -1;

    /** The total number of duplicates detected */
    private int numDuplicates = 0;

    /** The set of all read ends sorted by 5' start unclipped position.  Some read ends in this set may eventually be duplicates. */
    private final TreeSet<ReadEndsForMateCigar> set = new TreeSet<ReadEndsForMateCigar>(new MarkQueueComparator());

    /** Reads in the main set may occasionally have mates with the same chromosome, coordinate, and orientation, causing collisions
     * We store the 'best' end of the mate pair in the main set, and the other end in this set.  We only remove from this.pairSet when
     * we remove something from this.set.
     */
    private final TreeSet<ReadEndsForMateCigar> pairSet = new TreeSet<ReadEndsForMateCigar>(new MarkQueueComparator());

    /** Physical locations used for optical duplicate tracking.  This is only stored for paired end reads where both ends are mapped,
     * and when we see the first mate.
     */
    //private final Map<ReadEndsForMateCigar, LocationSet> locations = new HashMap<ReadEndsForMateCigar, LocationSet>();

    /** If we have two items that are the same with respect to being in the "set", then we must choose one.  The "one" will
     * eventually be the end that is not marked as a duplicate in most cases (see poll() for the exceptions).
     */
    private final Comparator<ReadEndsForMateCigar> comparator;

    /** temporary so we do not need to create many objects */
    private ReadEndsForMateCigar tmpReadEnds = null;

    public MarkQueue(final ScoringStrategy duplicateScoringStrategy) {
        comparator = new ReadEndsMCComparator(duplicateScoringStrategy);
    }

    /** Returns the number of duplicates detected */
    public int getNumDuplicates() { return this.numDuplicates; }

    /** The number of records currently in this queue. **/
    public int size() {
        return this.set.size() + this.pairSet.size();
    }

    public boolean isEmpty() {
        return this.set.isEmpty();
    }

    /** Sets the minimum genomic distance such that we can be assured that all duplicates have been considered. */
    public void setToMarkQueueMinimumDistance(final int toMarkQueueMinimumDistance) {
        this.toMarkQueueMinimumDistance = toMarkQueueMinimumDistance;
    }

    /** Returns the minimum genomic distance such that we can be assured that all duplicates have been considered. */
    public int getToMarkQueueMinimumDistance() { return this.toMarkQueueMinimumDistance; }

    /** Returns true if we should track this for optical duplicate detection, false otherwise */
    public boolean shouldBeInLocations(final ReadEndsForMateCigar current) {
        return (current.isPaired() && 0 == current.hasUnmapped);
    }

    /** Returns the set of read ends that should be considered for tracking optical duplicates. */
    public Set<ReadEnds> getLocations(final ReadEndsForMateCigar current) {
        // NB: only needed for pairs!!!
        if (!shouldBeInLocations(current)) throw new NotImplementedException();
        final Set<ReadEnds> locationSet = current.getReadEndSetForOpticalDuplicates();
        if (null == locationSet) throw new PicardException("Locations was empty: unexpected error");
        return locationSet;
    }

    /** Returns the first element in this queue */
    public ReadEndsForMateCigar peek() {
        return this.set.first();
    }

    /** Updates the duplication metrics given the provided duplicate */
    private void updateDuplicationMetrics(final ReadEndsForMateCigar duplicate, final DuplicationMetrics metrics) {
        // Update the duplication metrics
        if (!duplicate.getRecord().getReadPairedFlag() || duplicate.getRecord().getMateUnmappedFlag()) {
            ++metrics.UNPAIRED_READ_DUPLICATES;
        }
        else {
            ++metrics.READ_PAIR_DUPLICATES;// will need to be divided by 2 at the end
        }
        this.numDuplicates++;
    }

    /**
     * The poll method will return the read end that is *not* the duplicate of all comparable read ends that
     * have been seen.  All comparable read ends and the returned read end will have their seen duplicate flag set.  We use
     * the minimum genomic distance to determine when all the comparable reads have been examined.
     */
    public ReadEndsForMateCigar poll(final SamRecordTrackingBuffer outputBuffer,
                           final SAMFileHeader header,
                           final OpticalDuplicateFinder opticalDuplicateFinder,
                           final LibraryIdGenerator libraryIdGenerator) {
        /**
         * NB: we must remove all fragments or unpaired reads if this is a mapped pair.
         * NB: we must remove all fragments if this is an unpaired read.
         */

        final ReadEndsForMateCigar current = this.set.pollFirst();

        // Remove this record's comparable pair, if present.
        if (this.pairSet.contains(current)) { // the pair of this end is not a duplicate, if found
            final ReadEndsForMateCigar pair = this.pairSet.subSet(current, true, current, true).first();
            outputBuffer.setResultState(pair.getsamRecordIndex(), false); // you are not a duplicate!
            this.pairSet.remove(current);
            // NB: do not need to update metrics since this record is not a duplicate
        }

        // If we are a paired read end, we need to make sure we remove unpaired (if we are not also unpaired), as
        // well as fragments from the set, as they should all be duplicates.
        if (current.isPaired()) {

            // NB: only care about read1ReferenceIndex, read1Coordinate, and orientation in the set
            if (null == this.tmpReadEnds) { // initialize
                this.tmpReadEnds = new ReadEndsForMateCigar(header, current.getsamRecordIndex(), opticalDuplicateFinder, current.libraryId);
                this.tmpReadEnds.read2ReferenceIndex = this.tmpReadEnds.read2Coordinate = -1;
                this.tmpReadEnds.samRecordWithOrdinal = null; // nullify so we do not hold onto it forever
            }
            else {
                this.tmpReadEnds.read1ReferenceIndex = current.read1ReferenceIndex;
                this.tmpReadEnds.read1Coordinate = current.read1Coordinate;
            }

            // We should search for one of F/R
            if (current.orientation == ReadEnds.FF || current.orientation == ReadEnds.FR || current.orientation == ReadEnds.F) {
                this.tmpReadEnds.orientation = ReadEnds.F;
            }
            else {
                this.tmpReadEnds.orientation = ReadEnds.R;
            }

            // remove from the set fragments and unpaired, which only have two possible orientations
            //this.tmpReadEnds.orientation = orientation;
            if (this.set.contains(this.tmpReadEnds)) { // found in the set
                // get the duplicate read end
                final SortedSet<ReadEndsForMateCigar> sortedSet = this.set.subSet(this.tmpReadEnds, true, this.tmpReadEnds, true);
                if (1 != sortedSet.size()) throw new PicardException("SortedSet should have size one (has size " + sortedSet.size() + " )");
                final ReadEndsForMateCigar duplicate = sortedSet.first();

                /** mark as duplicate and set that it has been through duplicate marking
                 * duplicate.getRecord().setDuplicateReadFlag(true); HANDLED BY THE METHOD CALL BELOW*/
                outputBuffer.setResultState(duplicate.getsamRecordIndex(), true);

                // remove from the set
                this.set.remove(this.tmpReadEnds);

                // update the metrics
                updateDuplicationMetrics(duplicate, libraryIdGenerator.getMetricsByLibrary(libraryIdGenerator.getLibraryName(header, duplicate.getRecord())));
            }
        }

        // This read end is now ok to be emitted.  Track that it has been through duplicate marking.
        outputBuffer.setResultState(current.getsamRecordIndex(), false);

        return current;
    }

    /**
     * Add a record to the mark queue.
     */
    private boolean debug = false; // TODO: remove me
    public void add(final ReadEndsForMateCigar other,
                    final SamRecordTrackingBuffer outputBuffer,
                    final DuplicationMetrics metrics) {
        /**
         * OK this is the most complicated function in this class.  Please pay attention.
         */
        PhysicalLocationForMateCigarSet locationSet = null;
        boolean addToLocationSet = true; // only false if we have mates mapped to the same position
        ReadEndsForMateCigar duplicate = null;

        if (debug) System.err.print("TMQ add: " + other.getRecord().getSAMString());

        // 1. check the queue to see if there exists a comparable record at the location, if so compare and keep the best.
        // 2. add physical location info if paired

        /**
         * Check if we have a comparable record in our set.
         */
        if (this.set.contains(other)) { // a comparable record to "other" record already in the set
            if (debug) System.err.print("TMQ add: other contains " + other.getRecord().getSAMString());

            /**
             * Get the subset of records that are comparable, which should be of size one
             *
             * Sometimes, the ends that are comparable are in fact from the same pair.  In this case, we need to choose the best end
             * from the pair, and track the sub-optimal end.
             */
            final SortedSet<ReadEndsForMateCigar> sortedSet = this.set.subSet(other, true, other, true);
            if (1 != sortedSet.size()) throw new PicardException("SortedSet should have size one (has size " + sortedSet.size() + " )");
            final ReadEndsForMateCigar current = sortedSet.first();
            final String otherName = SAMUtils.getCanonicalRecordName(other.getRecord());
            final String currentName = SAMUtils.getCanonicalRecordName(current.getRecord());

            final int comparison = this.comparator.compare(current, other);

            /** Check if the ends to compare are in fact from the same pair */
            if (currentName.equals(otherName)) {
                /**
                 * "other" is paired-end mate of "current". We need to choose the best end to store in the main set.  We also have
                 * to maintain the pair set, as well as update the location set.
                 * */
                if (debug) System.err.print("TMQ add: other contains : name equals " + other.getRecord().getSAMString());
                if (0 < comparison) { // "other" is the best end. Swap for "current".
                    if (debug) System.err.print("TMQ add: other contains : name equals : swap : " + other.getRecord().getSAMString());
                    // Swap them in the set
                    this.set.remove(current);
                    this.set.add(other);
                    this.pairSet.add(current); // add "current" to the pairset
                    // Swap "current" and "other" in the locations
                    if (shouldBeInLocations(other)) {
                        locationSet = current.removeLocationSet();
                        locationSet.replace(current, other); // swap "current" and "other"
                        other.setLocationSet(locationSet);
                        addToLocationSet = false;
                    }
                } else { // other is less desirable. Store it in the pair set.
                    if (debug) System.err.print("TMQ add: other contains : name equals : no swap : " + other.getRecord().getSAMString());
                    this.pairSet.add(other); // add "other" to the pairset
                    if (shouldBeInLocations(current)) {
                        locationSet = current.getLocationSet();
                        addToLocationSet = false;
                    }
                }
                if (debug) System.err.print("TMQ add: other contains : name equals : DONE : " + other.getRecord().getSAMString());
            } else {
                /**
                 * "other" is a not the pair of "current" at the same location and must be compared against "current".
                 */
                if (debug) System.err.print("TMQ add: other contains : name not equals " + other.getRecord().getSAMString());

                if (0 < comparison) { // remove the current, and add the other in its place
                    if (debug) System.err.print("TMQ add: other contains : name not equals : 0<comparison " + other.getRecord().getSAMString());
                    if (shouldBeInLocations(current)) { // was this in the location set?
                        // NB we could also just check if locationSet == null after remove?
                        locationSet = current.removeLocationSet();
                    }
                    else { // make a new one
                        locationSet = new PhysicalLocationForMateCigarSet();
                    }
                    other.setLocationSet(locationSet); // update locations to use "other" as the identifier for the location set
                    // remove current and add the other
                    this.set.remove(current);
                    this.set.add(other);

                    // update the pair set in case current's pair is in that set
                    if (this.pairSet.contains(current)) {
                        final ReadEndsForMateCigar pair = this.pairSet.subSet(current, true, current, true).first();
                        this.pairSet.remove(current);
                        outputBuffer.setResultState(pair.getsamRecordIndex(), true); // track that this samRecordWithOrdinal has been through duplicate marking
                        updateDuplicationMetrics(pair, metrics);
                    }

                    // current is now a duplicate
                    duplicate = current;
                }
                else { // keep the current record, and the "other" is now a duplicate
                    if (shouldBeInLocations(current)) { // Get the location set
                        locationSet = current.getLocationSet();
                    }
                    // NB: else is technically not needed, since if this was not paired and the other one was, we would enter here and add it later

                    // other is a duplicate :/
                    duplicate = other;
                }
            }
        } else { // 'other' ReadEndMC is not in the main set, thus the first record at this location. Store it for now.
            if (debug) System.err.print("TMQ add: other not contains " + other.getRecord().getSAMString());

            if (shouldBeInLocations(other)) {
                locationSet = new PhysicalLocationForMateCigarSet();
                other.setLocationSet(locationSet);
            }
            this.set.add(other);
        }
        //System.err.println("\tSET SIZE=" + this.set.size());

        // add to the physical locations for optical duplicate tracking
        final SAMRecord record = other.getRecord();
        if (debug && record.getReadPairedFlag()) System.err.print("TESTING FOR LOCATION SET (" + record.getFirstOfPairFlag() + ") : " + record.getSAMString());
        if (record.getReadPairedFlag()
                && !record.getReadUnmappedFlag()
                && !record.getMateUnmappedFlag()
                && addToLocationSet) {
            if (null == locationSet) throw new PicardException("location set was null: " + record.getSAMString());
            locationSet.add(other);
            if (debug) System.err.print("TESTING FOR LOCATION SET OK (" + locationSet.size() + " " +
            other.getTile() + ":" + other.getX() + ":" + other.getY() + ") : " + record.getSAMString());
            if (debug) {
                for (ReadEnds loc : locationSet.getReadEnds()) {
                    System.err.println("LOCATION SET loc: + " + loc.getX() + ":" + loc.getY());
                }
            }
        }

        // if we have a duplicate, update it for duplicate tracking and update the metrics
        if (null != duplicate) {
            outputBuffer.setResultState(duplicate.getsamRecordIndex(), true);
            // count the duplicate metrics
            updateDuplicationMetrics(duplicate, metrics);
        }
    }
}
