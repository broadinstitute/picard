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

import com.google.common.math.LongMath;
import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Histogram;
import picard.PicardException;
import picard.sam.UmiMetrics;

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
    private UmiMetrics metrics;
    private boolean haveWeSeenFirstRead = false;
    private long duplicateSetsWithUmi = 0;
    private long duplicateSetsWithoutUmi = 0;
    private double expectedCollisions = 0;
    private int observedUmiBases = 0;

    private Histogram<String> observedUmis = new Histogram<>();
    private Histogram<String> inferredUmis = new Histogram<>();

    /**
     * Creates a UMI aware duplicate set iterator
     *
     * @param wrappedIterator       UMI aware duplicate set iterator is a wrapper
     * @param maxEditDistanceToJoin The edit distance between UMIs that will be used to union UMIs into groups
     * @param umiTag                The tag used in the bam file that designates the UMI
     * @param assignedUmiTag        The tag in the bam file that designates the assigned UMI
     */
    UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int maxEditDistanceToJoin,
                                 final String umiTag, final String assignedUmiTag, final boolean allowMissingUmis,
                                 final UmiMetrics metrics) {
        this.wrappedIterator = wrappedIterator;
        this.maxEditDistanceToJoin = maxEditDistanceToJoin;
        this.umiTag = umiTag;
        this.inferredUmiTag = assignedUmiTag;
        this.allowMissingUmis = allowMissingUmis;
        this.metrics = metrics;
        isOpen = true;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        isOpen = false;
        wrappedIterator.close();

        if (metrics.UMI_LENGTH > 0) {
            collectMetrics();
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

        List<DuplicateSet> duplicateSets = umiGraph.joinUmisIntoDuplicateSets(maxEditDistanceToJoin);

        // Collect statistics on numbers of observed and inferred UMIs
        // and total numbers of observed and inferred UMIs
        for (DuplicateSet ds : duplicateSets) {
            List<SAMRecord> records = ds.getRecords();
            SAMRecord representativeRead = ds.getRepresentative();

            String inferredUmi = representativeRead.getStringAttribute(inferredUmiTag);

            for (SAMRecord rec : records) {
                String currentUmi = rec.getStringAttribute(umiTag);

                if (currentUmi != null) {
                    // All UMIs should be the same length
                    if (!haveWeSeenFirstRead) {
                        metrics.UMI_LENGTH = currentUmi.length();
                        haveWeSeenFirstRead = true;
                    } else {
                        if (metrics.UMI_LENGTH != currentUmi.length()) {
                            throw new PicardException("UMIs of differing lengths were found.");
                        }
                    }

                    metrics.OBSERVED_BASE_ERRORS += hammingDistance(currentUmi, inferredUmi);
                    observedUmiBases += metrics.UMI_LENGTH;

                    observedUmis.increment(currentUmi);
                    inferredUmis.increment(inferredUmi);
                }
            }
        }
        duplicateSetsWithUmi += duplicateSets.size();
        duplicateSetsWithoutUmi++;

        // For each duplicate set estimate the number of expected UMI collisions
        double nWaysUmisCanDiffer = 0; // Number of ways two UMIs may contain errors, but be considered the same
        for (int k = 0; k <= maxEditDistanceToJoin; k++) {
            if (metrics.UMI_LENGTH > 0) {
                nWaysUmisCanDiffer = nWaysUmisCanDiffer + LongMath.binomial(metrics.UMI_LENGTH, k) * Math.pow(3.0, k);
            }
        }

        // The probability of two non-duplicate UMIs are drawn from a uniform distribution being correctly labeled as
        // not belonging to the same UMI family is given as, pCorrectlyLabeled = nWaysUmisCanDiffer / 4^metrics.UMI_LENGTH.
        // Estimate of probability all members in the duplicate set are correctly labeled is given by
        // pAllMembersCorrectlyLabeled = (1 - pCorrectlyLabeled)^duplicateSets.size().
        // The expected number of reads incorrectly labeled as belonging to a duplicate set is given as,
        // (1 - pAllMembersCorrectlyLabeled) * duplicateSets.size().
        expectedCollisions = expectedCollisions +
                (1 - Math.pow(1 - nWaysUmisCanDiffer / Math.pow(4, metrics.UMI_LENGTH), duplicateSets.size() - 1)) * duplicateSets.size();

        nextSetsIterator = duplicateSets.iterator();
    }

    private void collectMetrics() {
        metrics.OBSERVED_UNIQUE_UMIS = observedUmis.size();
        metrics.INFERRED_UNIQUE_UMIS = inferredUmis.size();

        metrics.OBSERVED_UMI_ENTROPY = effectiveNumberOfBases(observedUmis);
        metrics.INFERRED_UMI_ENTROPY = effectiveNumberOfBases(inferredUmis);

        metrics.DUPLICATE_SETS_WITH_UMI = duplicateSetsWithUmi;
        metrics.DUPLICATE_SETS_WITHOUT_UMI = duplicateSetsWithoutUmi;
        metrics.UMI_COLLISION_EST = expectedCollisions;

        metrics.UMI_COLLISION_Q = -10 * Math.log10(expectedCollisions / inferredUmis.size());
        metrics.estimateBaseQualities(observedUmiBases);
    }

    private double effectiveNumberOfBases(Histogram<?> observations) {
        double entropyBase4 = 0.0;

        double totalObservations = observations.getSumOfValues();
        for (Histogram.Bin observation : observations.values()) {
            double pObservation = observation.getValue() / totalObservations;
            entropyBase4 = entropyBase4 - pObservation * Math.log(pObservation);
        }

        // Convert to log base 4 so that the entropy is now a measure
        // of the effective number of DNA bases.  If we used log(2.0)
        // our result would be in bits.
        return entropyBase4 / Math.log(4.0);
    }
}
