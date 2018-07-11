/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

/**
 * This acts as an iterator over duplicate sets.  If a particular duplicate
 * set consists of records that contain UMIs this iterator breaks up a single
 * duplicate set into multiple duplicate based on the content of the UMIs.
 * Since there may be sequencing errors in the UMIs, this class allows for
 * simple error correction based on edit distances between the UMIs.
 *
 * @author fleharty
 */

package picard.sam.markduplicates;

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import picard.PicardException;

import java.util.*;

import static htsjdk.samtools.util.StringUtil.hammingDistance;

/**
 * UmiAwareDuplicateSetIterator is an iterator that wraps a duplicate set iterator
 * in such a way that each duplicate set may be broken up into subsets according
 * to UMIs in the records.  Some tolerance for errors in the UMIs is allowed, and
 * the degree of this is controlled by the maxEditDistanceToJoin parameter.
 */
class UmiAwareDuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    private final DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;
    private final int maxEditDistanceToJoin;
    private final String umiTag;
    private final String inferredUmiTag;
    private final boolean allowMissingUmis;
    private boolean isOpen = false;
    private Map<String, UmiMetrics> umiMetricsMap;
    private boolean haveWeSeenFirstRead = false;

    private long observedUmiBases = 0;

    /**
     * Creates a UMI aware duplicate set iterator
     *
     * @param wrappedIterator       Iterator of DuplicatesSets to use and break-up by UMI.
     * @param maxEditDistanceToJoin The edit distance between UMIs that will be used to union UMIs into groups
     * @param umiTag                The tag used in the bam file that designates the UMI
     * @param assignedUmiTag        The tag in the bam file that designates the assigned UMI
     * @param allowMissingUmis      Allow for SAM Records that do not have UMIs
     * @param umiMetricsMap         Map of UMI Metrics indexed by library name
     */
    UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int maxEditDistanceToJoin,
                                 final String umiTag, final String assignedUmiTag, final boolean allowMissingUmis,
                                 final Map<String, UmiMetrics> umiMetricsMap) {
        this.wrappedIterator = wrappedIterator;
        this.maxEditDistanceToJoin = maxEditDistanceToJoin;
        this.umiTag = umiTag;
        this.inferredUmiTag = assignedUmiTag;
        this.allowMissingUmis = allowMissingUmis;
        this.umiMetricsMap = umiMetricsMap;
        isOpen = true;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        isOpen = false;
        wrappedIterator.close();

        // Calculate derived fields for UMI metrics over each library
        for (final UmiMetrics metrics : umiMetricsMap.values()) {
            metrics.calculateDerivedFields();
        }
    }

    @Override
    public boolean hasNext() {
        if (!isOpen) {
            return false;
        } else {
            if (nextSetsIterator.hasNext() || wrappedIterator.hasNext()) {
                return true;
            } else {
                isOpen = false;
                return false;
            }
        }
    }

    @Override
    public DuplicateSet next() {
        if (!nextSetsIterator.hasNext()) {
            process(wrappedIterator.next());
        }
        return nextSetsIterator.next();
    }

    /**
     * Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
     * and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
     *
     * @param set Duplicate set that may be broken up into subsets according the UMIs
     */
    private void process(final DuplicateSet set) {

        // Ensure that the nextSetsIterator isn't already occupied
        if (nextSetsIterator.hasNext()) {
            throw new PicardException("nextSetsIterator is expected to be empty, but already contains data.");
        }

        final UmiGraph umiGraph = new UmiGraph(set, umiTag, inferredUmiTag, allowMissingUmis);

        // Get the UMI metrics for the library of this duplicate set, creating a new one if necessary.
        final String library = set.getRepresentative().getReadGroup().getLibrary();
        UmiMetrics metrics = umiMetricsMap.computeIfAbsent(library, UmiMetrics::new);

        final List<DuplicateSet> duplicateSets = umiGraph.joinUmisIntoDuplicateSets(maxEditDistanceToJoin);

        // Collect statistics on numbers of observed and inferred UMIs
        // and total numbers of observed and inferred UMIs
        for (final DuplicateSet ds : duplicateSets) {
            final List<SAMRecord> records = ds.getRecords();
            final SAMRecord representativeRead = ds.getRepresentative();
            final String inferredUmi = representativeRead.getStringAttribute(inferredUmiTag);

            for (final SAMRecord rec : records) {
                final String currentUmi = UmiUtil.getSanitizedUMI(rec, umiTag);

                if (currentUmi != null) {
                    // All UMIs should be the same length, the code presently does not support variable length UMIs.
                    // If the UMI contains a N, we don't want to include it in our other metrics but we still want
                    // to keep track of it.
                    if (currentUmi.contains("N")) {
                        metrics.addUmiObservationN();
                    } else {
                        if (!haveWeSeenFirstRead) {
                            metrics.MEAN_UMI_LENGTH = currentUmi.length();
                            haveWeSeenFirstRead = true;
                        } else {
                            if (metrics.MEAN_UMI_LENGTH != currentUmi.length()) {
                                throw new PicardException("UMIs of differing lengths were found.");
                            }
                        }

                        // Update UMI metrics associated with each record
                        // The hammingDistance between N and a base is a distance of 1. Comparing N to N is 0 distance.
                        metrics.OBSERVED_BASE_ERRORS += hammingDistance(currentUmi, inferredUmi);
                        observedUmiBases += currentUmi.length();
                        metrics.addUmiObservation(currentUmi, inferredUmi);
                    }
                }
            }
        }

        // Update UMI metrics associated with each duplicate set
        metrics.DUPLICATE_SETS_WITH_UMI += duplicateSets.size();
        metrics.DUPLICATE_SETS_IGNORING_UMI++;

        nextSetsIterator = duplicateSets.iterator();
    }
}
