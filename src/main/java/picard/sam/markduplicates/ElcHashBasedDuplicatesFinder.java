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

    // We split each read on several parts and each part has length hashLength.
    // hashLength = (min read length) / (number of errors + 1)
    private int hashLength = Integer.MAX_VALUE;

    ElcHashBasedDuplicatesFinder(double maxDiffRate, int maxReadLength, int minIdenticalBases,
                                 OpticalDuplicateFinder opticalDuplicateFinder) {
        super(maxDiffRate, maxReadLength, minIdenticalBases, opticalDuplicateFinder);
    }

    @Override
    void searchDuplicates(List<PairedReadSequence> sequences, Histogram<Integer> duplicationHisto,
                          Histogram<Integer> opticalHisto) {

        initHashLength(sequences);
        fillHashValues(sequences);
        populateDupCandidates(sequences);

        //Duplicate set to filter the treated PairedReadSequences out
        Set<PairedReadSequence> dupSet = new HashSet<>();
        for (PairedReadSequence lhs : sequences) {
            if (dupSet.contains(lhs)) continue;

            final List<PairedReadSequence> dupes = new ArrayList<>();

            for (PairedReadSequence rhs : lhs.dupCandidates) {
                if (dupSet.contains(rhs)) continue;

                if (isDuplicate(lhs, rhs)) {
                    dupes.add(rhs);
                }
            }

            dupSet.addAll(dupes);
            fillHistogram(duplicationHisto, opticalHisto, lhs, dupes);
        }
    }

    /**
     * Populate PairedReadSequence.dupCandidates with duplicate candidates
     * based on the same hash on the same position.
     */
    private void populateDupCandidates(List<PairedReadSequence> seqs) {
        // Contains hash value as a key and PairedReadSequences match this key as a value
        Map<Integer, List<PairedReadSequence>> readsByHash = new HashMap<>();

        // Iterate over all PairedReadSequence and split in sets according to the hash value 
        // in a certain position
        for (PairedReadSequence prs : seqs) {
            int[][] readHashValues = {prs.hashes1, prs.hashes2};
            for (int[] readHashValue : readHashValues) {
                for (int key : readHashValue) {
                    List<PairedReadSequence> dupCandidates = readsByHash.get(key);

                    if (dupCandidates == null) {
                        dupCandidates = new ArrayList<>();
                        readsByHash.put(key, dupCandidates);
                    }
                    dupCandidates.add(prs);
                }
            }
        }

        // Fill for each PairedReadSequence a set of possible duplicates
        for (List<PairedReadSequence> dupCandidates : readsByHash.values()) {
            for (PairedReadSequence prs : dupCandidates) {
                prs.dupCandidates.addAll(dupCandidates);
            }
        }
    }

    /**
     * Populate hashes for each PairedReadSequence.
     */
    private void fillHashValues(List<PairedReadSequence> sequences) {
        for (PairedReadSequence prs : sequences) {
            prs.initHashes(maxReadLength, hashLength, minIdenticalBases);
        }
    }

    /**
     * Calculate hash length based on minReadLength and hashNumber
     */
    private void initHashLength(List<PairedReadSequence> sequences) {
        if (sequences == null || sequences.isEmpty()) {
            return;
        }
        for (PairedReadSequence prs : sequences) {
            int minReadLength = Math.min(
                    Math.min(prs.read1.length, prs.read2.length) - minIdenticalBases,
                    maxReadLength - minIdenticalBases
            );

            //if read.length % shingle.length != 0, we have tail of the read which is not included in the hashValues,
            // and we need extra hash value to get prs.hashValues.length = (number of errors) + 1
            int hashNumber = (int) (minReadLength * maxDiffRate) + 1;
            int currentHashLength = minReadLength / hashNumber;
            if (hashLength > currentHashLength) {
                hashLength = currentHashLength;
            }
        }
    }

    private int getReadOffset(int k) {
        return k * hashLength + minIdenticalBases;
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
        int minHashesLength = Math.min(hashes1.length, hashes2.length);

        for (int k = 0; k < minHashesLength; ++k) {
            if (hashes1[k] != hashes2[k]) {
                errors += compareHashes(read1, read2, getReadOffset(k), getReadOffset((k + 1)));
                if (errors > maxErrors) {
                    return errors;
                }
            }
        }

        if (minReadLength > getReadOffset(minHashesLength)) {
            errors += compareHashes(read1, read2, getReadOffset(minHashesLength), minReadLength);
        }

        return errors;
    }

    /**
     * Compare bases corresponding to the hashes.
     *
     * @return errors number for current part of the read
     */
    private int compareHashes(byte[] read1, byte[] read2, int start, int stop) {
        int errors = 0;
        for (int i = start; i < stop; ++i) {
            if (read1[i] != read2[i]) {
                errors++;
            }
        }
        return errors;
    }
}
