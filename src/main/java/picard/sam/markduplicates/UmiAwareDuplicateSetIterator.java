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
import htsjdk.samtools.DuplicateSetIterator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import picard.PicardException;

import java.util.*;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.counting;

/**
 * Created by fleharty on 5/23/16.
 */
public class UmiAwareDuplicateSetIterator implements CloseableIterator<DuplicateSet> {
    private final DuplicateSetIterator wrappedIterator;
    private Iterator<DuplicateSet> nextSetsIterator;
    private final int editDistanceToJoin;
    private final boolean addInferredUmi;
    private final String umiTag;
    private final String inferredUmiTag;

    public UmiAwareDuplicateSetIterator(final DuplicateSetIterator wrappedIterator, final int editDistanceToJoin, final boolean addInferredUmi, final String umiTag, final String inferredUmiTag) {
        this.wrappedIterator = wrappedIterator;
        this.editDistanceToJoin = editDistanceToJoin;
        this.addInferredUmi = addInferredUmi;
        this.umiTag = umiTag;
        this.inferredUmiTag = inferredUmiTag;
        nextSetsIterator = Collections.emptyIterator();
    }

    @Override
    public void close() {
        wrappedIterator.close();
    }

    @Override
    public boolean hasNext() {
        return nextSetsIterator.hasNext() || wrappedIterator.hasNext();
    }

    @Override
    public DuplicateSet next() {
        if(!nextSetsIterator.hasNext())
            process(wrappedIterator.next());

        return nextSetsIterator.next();
    }

    // Takes a duplicate set and breaks it up into possible smaller sets according to the UMI,
    // and updates nextSetsIterator to be an iterator on that set of DuplicateSets.
    private void process(final DuplicateSet set) {

        List<SAMRecord> records = set.getRecords();

        // If any records are missing the UMI_TAG proceed as if there were no UMIs
        // and return nextSetsIterator without breaking it up into smaller sets.
        for(SAMRecord rec : records) {
            if(rec.getStringAttribute(umiTag) == null) {
                nextSetsIterator = Collections.singleton(set).iterator();
                return;
            }
        }

        // Count all the uniquely observed UMI sequences
        Map<String, Long> umiFrequencies = records.stream()
                        .collect(Collectors.groupingBy(p -> p.getStringAttribute(umiTag),
                                counting()));

        UmiGraph umiGraph = new UmiGraph();

        // Add each of the uniquely observed UMI sequences to the graph
        int firstGroup = 0;
        for (Map.Entry<String, Long> umiFrequency : umiFrequencies.entrySet()) {
            umiGraph.addNode(new UmiSequence(umiFrequency.getKey(), umiFrequency.getValue()));
        }
        umiGraph.processNodes();

        Map<Integer, List<UmiSequence>> umisByDuplicateSet = umiGraph.getUmisByDuplicateSet();

        // Construct DuplicateSetList
        List<DuplicateSet> duplicateSetList = new ArrayList<>();
        for(int i = 0;i < umiGraph.duplicateSetCounter;i++) {
            DuplicateSet e = new DuplicateSet();
            duplicateSetList.add(e);
        }

        for(SAMRecord rec : records) {
            String umi = records.get(records.indexOf(rec)).getStringAttribute(umiTag);

            // Figure out which group this belongs to
            int recordGroup = umiGraph.node.get(umi).duplicateSet;
            duplicateSetList.get(recordGroup).add(records.get(records.indexOf(rec)));
        }

        // Optionally add the inferred Umi
        if(addInferredUmi) {
            for(DuplicateSet ds : duplicateSetList) {
                // For each duplicate set identify the most common UMI
                List<SAMRecord> recordsFromDuplicateSet = ds.getRecords();

                // Assign the inferred UMI to all reads in the current group
                for(SAMRecord record : recordsFromDuplicateSet) {
                    record.setAttribute(inferredUmiTag, umiGraph.getInferredUmiFromOriginalSequence(record.getStringAttribute(umiTag)));
                }
            }
        }

        nextSetsIterator = duplicateSetList.iterator();
    }

    private int getEditDistance(final String s1, final String s2) {
        if(s1 == null || s2 == null) {
            throw new PicardException("Attempt to compare two incomparable UMIs.  At least one of the UMIs was null.");
        }
        if(s1.length() != s2.length()) {
            throw new PicardException("Barcode " + s1 + " and " + s2 + " do not have matching lengths.");
        }
        int count = 0;
        for(int i = 0;i < s1.length();i++) {
            if(s1.charAt(i) != s2.charAt(i)) {
                count++;
            }
        }
        return count;
    }

    private class UmiSequence {
        private final String sequence;
        private String inferredSequence;
        private final Long occupancy;
        private Integer duplicateSet;

        UmiSequence(final String sequence, final Long occupancy) {
            this.sequence = sequence;
            this.occupancy = occupancy;
        }
    }

    private class UmiGraph {
        int duplicateSetCounter = 0;
        private final Map<String, UmiSequence> node = new HashMap<>();

        void addNode(final UmiSequence umi) {
            node.put(umi.sequence, umi);
        }

        void processNodes() {
            boolean isFirstEntry = true;
            for(Map.Entry<String, UmiSequence> umi1 : node.entrySet()) {

                if(isFirstEntry) {
                    isFirstEntry = false;
                    umi1.getValue().duplicateSet = duplicateSetCounter;
                    duplicateSetCounter++;
                }
                // Test to see if this is similar to something I've seen before (that has been assigned to a duplicate set)
                for (Map.Entry<String, UmiSequence> umi2 : node.entrySet()) {
                    // If it is similar (within edit distance editDistanceToJoin) and not the same umi
                    // join it to the appropriate existing duplicate set
                    if (getEditDistance(umi1.getKey(), umi2.getKey()) <= editDistanceToJoin) {
                        // Since it is similar to something we have seen before we
                        // assign it to it's appropriate duplicateSet
                        if(umi1.getValue().duplicateSet != null) {
                            umi2.getValue().duplicateSet = umi1.getValue().duplicateSet;
                        } else if(umi2.getValue().duplicateSet != null) {
                            umi1.getValue().duplicateSet = umi2.getValue().duplicateSet;
                        }
                    }
                }
                // If we haven't added the UMI sequence to a duplicate set yet, then we need to put it in a
                // in its own duplicate set.
                if (umi1.getValue().duplicateSet == null) {
                    // If we haven't seen a similar UMI before, we add it, and create a
                    // new duplicate set for it.
                    umi1.getValue().duplicateSet = duplicateSetCounter;
                    //node.put(umi1.getValue().sequence, umi1.getValue());
                    duplicateSetCounter++;
                }
            }
        }

        private String getInferredUmiFromOriginalSequence(String seq) {
            return node.get(seq).inferredSequence;
        }

        Map<Integer, List<UmiSequence>> getUmisByDuplicateSet() {
            Map<Integer, List<UmiSequence>> umisByDuplicateSet = new HashMap<>();

            // Create UmiSequence Lists for each group that exists
            for(int i = 0;i < duplicateSetCounter;i++) {
                umisByDuplicateSet.put(i, new ArrayList<>());
            }

            // Assign each UmiSequence to a Duplicate set
            for(Map.Entry<String, UmiSequence> s : node.entrySet()) {
                umisByDuplicateSet.get(s.getValue().duplicateSet).add(s.getValue());
            }

            // Infer the Umi from each duplicate set
            for(Map.Entry<Integer, List<UmiSequence>> s : umisByDuplicateSet.entrySet())  {

                // Find most common UMI sequence in this duplicate set
                Long maxOccupancy = new Long(0);
                String mostCommonUmiSequence = null;
                for(UmiSequence seq : s.getValue()) {
                    if(seq.occupancy > maxOccupancy) {
                        maxOccupancy = seq.occupancy;
                        mostCommonUmiSequence = seq.sequence;
                    }
                }

                // Set the most common UMI sequence to the inferred UMI for all UmiSequences
                // in this duplicate set
                for(UmiSequence seq : s.getValue()) {
                    seq.inferredSequence = mostCommonUmiSequence;
                }
            }
            return umisByDuplicateSet;
        }
    }
}

