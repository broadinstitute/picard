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

import htsjdk.samtools.DuplicateSet;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.StringUtils;
import picard.PicardException;
import picard.util.GraphUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.counting;
/**
 * UmiGraph is used to identify UMIs that come from the same original source molecule.  The assumption
 * is that UMIs with small edit distances are likely to be read errors on the sequencer rather than
 * distinct molecules.
 *
 * The algorithm used here is to join all pairs of UMIs that are within maxEditDistanceToJoin.  It is possible
 * for a set of UMIs A, B and C to all be considered as part of the same source molecule even if two of the UMIs
 * have a Hamming distance larger than maxEditDistanceToJoin.  Suppose A = "ATCC", B = "AACC", and C = "AACG"
 * and maxEditDistanceToJoin = 1.  In this case, A and B are 1 Hamming distance so they are joined, and B and C
 * are 1 Hamming distance so they are joined.  Because A and B are joined and because B and C are joined, this results
 * in A and C being joined even though they have a distance of 2.
 *
 * @author fleharty
 */
public class UmiGraph {
    private final List<SAMRecord> records;      // SAMRecords from the original duplicate set considered to break up by UMI
    private final Map<String, Long> umiCounts;  // Map of UMI sequences and how many times they have been observed
    private final int[] duplicateSetID;         // ID of the duplicate set that the UMI belongs to, the index is the UMI ID
    private final String[] umi;                 // Sequence of actual UMI, the index is the UMI ID
    private final int numUmis;                  // Number of observed UMIs
    private final String umiTag;                // UMI tag used in the SAM/BAM/CRAM file ie. RX
    private final String assignedUmiTag;        // Assigned UMI tag used in the SAM/BAM/CRAM file ie. MI
    private final boolean allowMissingUmis;     // Allow for missing UMIs

    public UmiGraph(DuplicateSet set, String umiTag, String assignedUmiTag, boolean allowMissingUmis) {
        this.umiTag = umiTag;
        this.assignedUmiTag = assignedUmiTag;
        this.allowMissingUmis = allowMissingUmis;
        records = set.getRecords();

        // First ensure that all the reads have a UMI, if any reads are missing a UMI throw an exception unless allowMissingUmis is true
        for (SAMRecord rec : records) {
            if (UmiUtil.getSanitizedUMI(rec, umiTag) == null) {
                if (allowMissingUmis) {
                    rec.setAttribute(umiTag, "");
                } else {
                    throw new PicardException("Read " + rec.getReadName() + " does not contain a UMI with the " + umiTag + " attribute.");
                }
            }
        }

        // Count the number of times each UMI occurs
        umiCounts = records.stream().collect(Collectors.groupingBy(p -> UmiUtil.getSanitizedUMI(p, umiTag), counting()));

        // At first we consider every UMI as if it were its own duplicate set
        numUmis = umiCounts.size();
        umi = new String[numUmis];
        duplicateSetID = IntStream.rangeClosed(0, numUmis-1).toArray();

        int i = 0;
        for (String key : umiCounts.keySet()) {
            umi[i] = key;
            i++;
        }
    }

    List<DuplicateSet> joinUmisIntoDuplicateSets(final int maxEditDistanceToJoin) {
        // Compare all UMIs to each other.  If they are within maxEditDistanceToJoin
        // join them to the same duplicate set using the union-find algorithm.

        GraphUtils.Graph<Integer> umiGraph = new GraphUtils.Graph<>();
        for (int i = 0; i < numUmis; i++) {
            umiGraph.addNode(i);
            for (int j = i + 1; j < numUmis; j++) {
                if (StringUtil.isWithinHammingDistance(umi[i], umi[j], maxEditDistanceToJoin)) {
                    umiGraph.addEdge(i, j);
                }
            }
        }

        final Map<Integer, Integer> umiClusterMap = umiGraph.cluster();
        // This ensures that all duplicate sets have unique IDs.  During Union-Find a tree is constructed
        // where each UMI points to parent UMI.  This ensures that all UMIs that belong to the same duplicate
        // set point to the same parent UMI.  Note that the parent UMI is only used as a representative UMI and
        // is not at all related to the assigned UMI.
        for (int i = 0; i < numUmis; i++) {
            duplicateSetID[i] = umiClusterMap.get(i);
        }

        final Map<Integer, List<SAMRecord>> duplicateSets = new HashMap<>();

        // Assign UMIs to duplicateSets
        final Map<String, Integer> duplicateSetsFromUmis = getDuplicateSetsFromUmis();
        for (SAMRecord rec : records) {
            final String umi = UmiUtil.getSanitizedUMI(rec, umiTag);
            final Integer duplicateSetIndex = duplicateSetsFromUmis.get(umi);

            if (duplicateSets.containsKey(duplicateSetIndex)) {
                duplicateSets.get(duplicateSetIndex).add(rec);
            } else {
                final List<SAMRecord> n = new ArrayList<>();
                n.add(rec);
                duplicateSets.put(duplicateSetIndex, n);
            }
        }

        final List<DuplicateSet> duplicateSetList = new ArrayList<>();
        for (final Map.Entry<Integer, List<SAMRecord>> entry : duplicateSets.entrySet()) {
            final DuplicateSet ds = new DuplicateSet();
            final List<SAMRecord> recordList = entry.getValue();

            // Add records to the DuplicateSet
            recordList.forEach(ds::add);

            // For a particular duplicate set, identify the most common UMI
            // and use this as an assigned UMI.
            long maxCount = 0;
            String assignedUmi = null;
            String fewestNUmi = null;
            long nCount = 0;

            for (SAMRecord rec : recordList) {
                final String umi = UmiUtil.getSanitizedUMI(rec, umiTag);

                // If there is another choice, we don't want to choose the UMI with a N
                // as the assignedUmi
                if (umi.contains("N")) {

                    int count = StringUtils.countMatches(umi, "N");
                    // If an UMI containing a N hasn't been seen before, the current UMI is now the fewestNUmi
                    if (nCount == 0) {
                        nCount = count;
                        fewestNUmi = umi;
                    } else if (count < nCount) { // If N containing UMI already exists, reset it if we find a lower one
                        nCount = count;
                        fewestNUmi = umi;
                    }
                } else if (umiCounts.get(umi) > maxCount) {
                    maxCount = umiCounts.get(umi);
                    assignedUmi = umi;
                }
            }

            // If we didn't find a Umi w/o a N to represent
            // then choose the one with the fewest Ns
            if (assignedUmi == null) { assignedUmi = fewestNUmi; }

            // Set the records to contain the assigned UMI
            for (final SAMRecord rec : recordList) {
                if (allowMissingUmis && rec.getStringAttribute(umiTag).isEmpty()) {
                    // The SAM spec doesn't support empty tags, so we set it to null if it is empty.
                    rec.setAttribute(umiTag, null);
                } else {
                    rec.setAttribute(assignedUmiTag, assignedUmi);
                }
            }

            duplicateSetList.add(ds);
        }

        return duplicateSetList;
    }

    // Create a map that maps a umi to the duplicateSetID
    private Map<String, Integer> getDuplicateSetsFromUmis() {
        final Map<String, Integer> duplicateSetsFromUmis = new HashMap<>();
        for (int i = 0; i < duplicateSetID.length; i++) {
            duplicateSetsFromUmis.put(umi[i], duplicateSetID[i]);
        }
        return duplicateSetsFromUmis;
    }
}
