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

import java.util.*;

import static picard.sam.markduplicates.EstimateLibraryComplexity.PairedReadSequence;

/**
 * Algorithm for search duplicates is used in EstimateLibraryComplexity only. It splits each PairedReadSequence
 * on maximum error number + 1 and hash them. To find a duplicate we need to compare hashes and if they are similar
 * we compare bases corresponding to the hashes.
 *
 * * @author Pavel_Silin@epam.com, EPAM Systems, Inc. <www.epam.com>
 */
class ElcHashBasedDuplicatesFinder extends ElcDuplicatesFinder {

    // We split each read on several parts, each part represent set of nucleotides,
    // which behind each other on numberOfHashesInGroup.
    private int numberOfHashesInGroup = -1; //initial value
    private int minReadLenInGroup = Integer.MAX_VALUE; //initial value

    private Map<Integer, List<PairedReadSequence>> readsByHashInGroup;

    ElcHashBasedDuplicatesFinder(double maxDiffRate, int maxReadLength, int minIdenticalBases,
                                 OpticalDuplicateFinder opticalDuplicateFinder) {
        super(maxDiffRate, maxReadLength, minIdenticalBases, opticalDuplicateFinder);
        readsByHashInGroup = new HashMap<>();
    }

    @Override
    void searchDuplicates(List<PairedReadSequence> sequences, Histogram<Integer> duplicationHisto,
                          Histogram<Integer> opticalHisto) {

        initHashLength(sequences);
        fillHashValues(sequences);
        populateDupCandidates(sequences);

        //Duplicate set to filter the treated PairedReadSequences out
        final Set<PairedReadSequence> dupSet = new HashSet<>();

        for (PairedReadSequence lhs : sequences) {
            if (dupSet.contains(lhs)) continue;
            final List<PairedReadSequence> dupes = new ArrayList<>();
            for (PairedReadSequence rhs : getSimilarReads(lhs)) {
                if (dupSet.contains(rhs)) continue;
                if (isDuplicate(lhs, rhs)) {
                    dupes.add(rhs);
                }
            }
            dupSet.addAll(dupes);
            fillHistogram(duplicationHisto, opticalHisto, lhs, dupes);
        }
    }

    private Set<PairedReadSequence> getSimilarReads(PairedReadSequence lhs) {
        final Set<PairedReadSequence> toCheck = new HashSet<>();
        for (int[] hashesForRead: new int[][]{lhs.hashes1, lhs.hashes2}){
            for (int hash : hashesForRead) {
                List<PairedReadSequence> readsWithSameHash = readsByHashInGroup.get(hash);
                // if size == 1, then this list contains only lhs read
                if (readsWithSameHash.size() > 1) {
                    toCheck.addAll(readsWithSameHash);
                }
            }
        }
        return toCheck;
    }

    /**
     * Populate PairedReadSequence.dupCandidates with duplicate candidates
     * based on the same hash on the same position.
     */
    private void populateDupCandidates(List<PairedReadSequence> seqs) {
        // Contains hash value as a key and PairedReadSequences match this key as a value
        readsByHashInGroup.clear();

        // Iterate over all PairedReadSequence and split in sets according to the hash value
        // in a certain position
        for (PairedReadSequence prs : seqs) {
            int[][] readHashValues = {prs.hashes1, prs.hashes2};
            for (int[] readHashValue : readHashValues) {
                for (int key : readHashValue) {
                    List<PairedReadSequence> dupCandidates
                            = readsByHashInGroup.computeIfAbsent(key, k -> new ArrayList<>());
                    dupCandidates.add(prs);
                }
            }
        }
    }

    /**
     * Populate hashes for each PairedReadSequence.
     */
    private void fillHashValues(List<PairedReadSequence> sequences) {
        for (PairedReadSequence prs : sequences) {
            prs.initHashes(numberOfHashesInGroup, minIdenticalBases, minReadLenInGroup);
        }
    }

    /**
     * Calculate hash length based on minReadLength and numberOfHashesInGroup
     */
    private void initHashLength(List<PairedReadSequence> sequences) {
        for (PairedReadSequence prs : sequences) {
            int minReadLength = Math.min(Math.min(prs.read1.length, prs.read2.length), maxReadLength);

            //search maximum number of hashes (if length of reads is similar number of hashes will be similar too)
            int numberOfHashes = (int) ((minReadLength - minIdenticalBases) * maxDiffRate) + 1;
            if (numberOfHashes > numberOfHashesInGroup) {
                numberOfHashesInGroup = numberOfHashes;
            }

            //search minimum length of read for calculating hashes value
            if (minReadLenInGroup > minReadLength) {
                minReadLenInGroup = minReadLength;
            }
        }
    }

    private boolean isDuplicate(final PairedReadSequence lhs, final PairedReadSequence rhs) {
        // self is not duplicate
        if (lhs == rhs) {
            return false;
        }

        final int read1Length = minLength(lhs.read1, rhs.read1);
        final int read2Length = minLength(lhs.read2, rhs.read2);
        final int maxErrors = (int) Math.floor((read1Length + read2Length) * maxDiffRate);

        int errors = compareReadToRead(
                lhs.read1, lhs.hashes1,
                rhs.read1, rhs.hashes1,
                maxErrors
        );

        if (errors > maxErrors) {
            return false;
        }

        errors += compareReadToRead(
                lhs.read2, lhs.hashes2,
                rhs.read2, rhs.hashes2,
                maxErrors
        );

        return errors <= maxErrors;
    }
    /**
     * Compare hashes and if they are similar we compare bases corresponding to the hashes.
     */
    private int compareReadToRead(byte[] read1, int[] hashes1, byte[] read2, int[] hashes2, int maxErrors) {
        int errors = 0;
        final int minReadLength = minLength(read1, read2);

        for (int k = 0; k < numberOfHashesInGroup; ++k) {
            if (hashes1[k] != hashes2[k]) {
                errors += compareHashes(read1, read2, k);
                if (errors > maxErrors) {
                    return errors;
                }
            }
        }

        if (minReadLength > minReadLenInGroup) {
            errors += compareTails(read1, read2, minReadLenInGroup, minReadLength);
        }

        return errors;
    }

    /**
     * Compare bases corresponding to the hashes.
     *
     * @return errors number for current part of the read
     */
    private int compareHashes(byte[] read1, byte[] read2, int k) {
        int errors = 0;
        int position = minIdenticalBases + k;
        while (position < minReadLenInGroup) {
            if (read1[position] != read2[position]) {
                errors++;
            }
            position += numberOfHashesInGroup;
        }
        return errors;
    }

    /**
     * Compare bases corresponding to the hashes.
     *
     * @return errors number for current part of the read
     */
    private int compareTails(byte[] read1, byte[] read2, int start, int stop) {
        int errors = 0;
        for (int i = start; i < stop; ++i) {
            if (read1[i] != read2[i]) {
                errors++;
            }
        }
        return errors;
    }
}
