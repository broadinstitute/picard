/*
 * The MIT License
 *
 * Copyright (c) 2013 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

/**
 * Holds all the hits (alignments) for a read or read pair.  For single-end reads, all the alignments are
 * stored in firstOfPairOrFragment list.  For paired-end, alignments are stored in firstOfPairOrFragment list and
 * secondOfPair list.
 *
 * If there is more than one alignment, the selected PrimaryAlignmentSelectionStrategy is used to decide which
 * alignment should be marked as primary.  The rest are marked as secondary.
 *
 * When AbstractAlignmentMerger emits these reads, for paired end reads it assumes that the ith first end alignment
 * and the ith second end alignment are correlated for the purpose of setting mate information in the SAMRecord.
 * If it is not appropriate for the ends to be linked like that, then the alignments should be staggered in
 * the lists so that there is a null in the other end list for each alignment.  E.g. for the firstOfPair(5),
 * secondOfPair(5) should be null in order not to set the mate info.  In that case the mate info will indicate that
 * the other end is unmapped.
 */
class HitsForInsert {

    public enum NumPrimaryAlignmentState {
        NONE, ONE, MORE_THAN_ONE
    }

    // These are package-visible to make life easier for the PrimaryAlignmentSelectionStrategies.
    final List<SAMRecord> firstOfPairOrFragment = new ArrayList<>();
    final List<SAMRecord> secondOfPair = new ArrayList<>();

    private final List<SAMRecord> supplementalFirstOfPairOrFragment = new ArrayList<>();
    private final List<SAMRecord> supplementalSecondOfPair = new ArrayList<>();

    /**
     * @throws IllegalStateException if numHits() == 0
     */
    public String getReadName() {
        return getRepresentativeRead().getReadName();
    }

    /**
     * @throws IllegalStateException if numHits() == 0
     */
    public boolean isPaired() {
        return getRepresentativeRead().getReadPairedFlag();
    }

    public SAMRecord getRepresentativeRead() {
        for (final SAMRecord rec : firstOfPairOrFragment) {
            if (rec != null) return rec;
        }
        for (final SAMRecord rec : secondOfPair) {
            if (rec != null) return rec;
        }
        throw new IllegalStateException("Should not be called if numHits == 0");
    }

    /**
     * Note that a single alignment for each end of a read pair is counted as a single hit.
     */
    public int numHits() {
        return Math.max(firstOfPairOrFragment.size(), secondOfPair.size());
    }

    /** True if either the first or second of pair has supplementary alignments, otherwise false. */
    public boolean hasSupplementalHits() {
        return !(this.supplementalFirstOfPairOrFragment.isEmpty() && this.supplementalSecondOfPair.isEmpty());
    }

    /**
     * @return Returns the ith hit for the first end, or null if the first end is not aligned.
     */
    public SAMRecord getFirstOfPair(final int i) {
        if (i >= firstOfPairOrFragment.size()) {
            return null;
        } else {
            return firstOfPairOrFragment.get(i);
        }
    }

    public void addFirstOfPairOrFragment(final SAMRecord rec) {
        firstOfPairOrFragment.add(rec);
    }

    public void addSecondOfPair(final SAMRecord rec) {
        secondOfPair.add(rec);
    }

    public void addSupplementalFirstOfPairOrFragment(final SAMRecord rec) {
        supplementalFirstOfPairOrFragment.add(rec);
    }

    public void addSupplementalSecondOfPair(final SAMRecord rec) {
        supplementalSecondOfPair.add(rec);
    }

    /**
     * @return The ith hit for a un-paired read.  Never returns null.
     * Do not call if paired read.
     */
    public SAMRecord getFragment(final int i) {
        final SAMRecord samRecord = firstOfPairOrFragment.get(i);
        if (samRecord.getReadPairedFlag()) {
            throw new UnsupportedOperationException("getFragment called for paired read: " + samRecord.toString());
        }
        return samRecord;
    }

    /**
     * @return Returns the ith hit for the second end, or null if the second end is not aligned.
     */
    public SAMRecord getSecondOfPair(final int i) {
        if (i >= secondOfPair.size()) {
            return null;
        } else {
            return secondOfPair.get(i);
        }
    }

    /**
     * Get the index of the first primary we see in the list of hits (either read1 or read2).
     * NOTE: if the PrimaryAlignmentSelectionStrategy has not been run, the returned value may not represent the ONLY primary.
     *
     * @return the index, or -1 if no primary was found.
     */
    public int getIndexOfEarliestPrimary() {
        for (int i = 0; i < numHits(); ++i) {
            final SAMRecord firstAligned = getFirstOfPair(i);
            final SAMRecord secondAligned = getSecondOfPair(i);
            final boolean isPrimaryAlignment = (firstAligned != null && !firstAligned.isSecondaryOrSupplementary()) ||
                    (secondAligned != null && !secondAligned.isSecondaryOrSupplementary());
            if (isPrimaryAlignment) return i;
        }
        return -1;
    }

    /**
     * Used by PrimaryAlignmentSelectionStrategy to set all alignments to not primary, except for the one specified by the argument.
     * If paired, and set the alignment for both ends if there is an alignment for both ends, otherwise just for the end for which
     * there is an alignment at the given index.
     * @param primaryAlignmentIndex
     */
    public void setPrimaryAlignment(final int primaryAlignmentIndex) {
        if (primaryAlignmentIndex < 0 || primaryAlignmentIndex >= this.numHits()) {
            throw new IllegalArgumentException("primaryAlignmentIndex(" + primaryAlignmentIndex +
                    ") out of range for numHits(" + numHits() + ")");
        }
        // Set all alignment to be not primary except the selected one.
        for (int i = 0; i < this.numHits(); ++i) {
            final boolean notPrimary = (i != primaryAlignmentIndex);
            if (this.getFirstOfPair(i) != null) {
                this.getFirstOfPair(i).setNotPrimaryAlignmentFlag(notPrimary);
            }
            if (this.getSecondOfPair(i) != null) {
                this.getSecondOfPair(i).setNotPrimaryAlignmentFlag(notPrimary);
            }
        }
    }

    /**
     * Some alignment strategies expect to receive alignments for ends that are coordinated by
     * hit index (HI) tag.  This method lines up alignments for each end by HI tag value, and if there is
     * no corresponding alignment for an alignment, there is a null in the array at that slot.
     *
     * This method then renumbers the HI values so that they start at zero and have no gaps, because some
     * reads may have been filtered out.
     */
    public void coordinateByHitIndex() {
        // Sort by HI value, with reads with no HI going at the end.
        final Comparator<SAMRecord> comparator = Comparator.comparing(
                record -> record.getIntegerAttribute(SAMTag.HI.name()),
                Comparator.nullsLast(Comparator.naturalOrder())
        );
        firstOfPairOrFragment.sort(comparator);
        secondOfPair.sort(comparator);

        // Insert nulls as necessary in the two lists so that correlated alignments have the same index
        // and uncorrelated alignments have null in the other list at the corresponding index.
        for (int i = 0; i < Math.min(firstOfPairOrFragment.size(), secondOfPair.size()); ++i) {
            final Integer leftHi = firstOfPairOrFragment.get(i).getIntegerAttribute(SAMTag.HI.name());
            final Integer rightHi = secondOfPair.get(i).getIntegerAttribute(SAMTag.HI.name());
            if (leftHi != null) {
                if (rightHi != null) {
                    if (leftHi < rightHi) secondOfPair.add(i, null);
                    else if (rightHi < leftHi) firstOfPairOrFragment.add(i, null);
                    // else the are correlated
                }
            } else if (rightHi != null) {
                firstOfPairOrFragment.add(i, null);
            } else {
                // Both alignments do not have HI, so push down the ones on the right.
                // Right is arbitrarily picked to push down.
                secondOfPair.add(i, null);
            }
        }
    }

    /**
     * Renumber any correlated alignments, and remove hit index if no correlated read.
     */
    private void renumberHitIndex() {
        for (int i = 0; i < numHits(); ++i) {
            final SAMRecord first = getFirstOfPair(i);
            final SAMRecord second = getSecondOfPair(i);
            if (first != null && second != null) {
                first.setAttribute(SAMTag.HI.name(), i);
                second.setAttribute(SAMTag.HI.name(), i);
            } else if (first != null) {
                first.setAttribute(SAMTag.HI.name(), null);
            } else {
                second.setAttribute(SAMTag.HI.name(), null);
            }
        }
    }

    /**
     * This method lines up alignments in firstOfPairOrFragment and secondOfPair lists so that the first and the second
     * records of each pair have the same index.
     */
    void coordinateByMate() {
        final List<SAMRecord> newFirstOfPairOrFragment = new ArrayList<>();
        final List<SAMRecord> newSecondOfPair = new ArrayList<>();

        for (int i = 0; i < firstOfPairOrFragment.size(); ++i) {
            final SAMRecord first = firstOfPairOrFragment.get(i);
            newFirstOfPairOrFragment.add(i, first);
            newSecondOfPair.add(null);

            for (int k = 0; k < secondOfPair.size(); ++k) {
                final SAMRecord second = secondOfPair.get(k);
                if (arePair(first, second)) {
                    newSecondOfPair.set(i, second);
                    secondOfPair.set(k, null);
                    break;
                }
            }
        }

        for (int i = 0; i < secondOfPair.size(); ++i) {
            if (secondOfPair.get(i) != null) {
                newFirstOfPairOrFragment.add(i, null);
                newSecondOfPair.add(i, secondOfPair.get(i));
            }
        }

        firstOfPairOrFragment.clear();
        firstOfPairOrFragment.addAll(newFirstOfPairOrFragment);
        secondOfPair.clear();
        secondOfPair.addAll(newSecondOfPair);

        renumberHitIndex();
    }

    /**
     * Identifies whether records present pairwise alignment or not.
     *
     * It is unnecessary to check QNAME here cause an object of {@code HitsForInsert}
     * presents only records with the same QNAME.
     *
     * @return {@code true} if records belong to the same pairwise alignment
     */
    private static boolean arePair(final SAMRecord first, final SAMRecord second) {
        return first != null && second != null
                && first.getMateAlignmentStart() == second.getAlignmentStart()
                && first.getAlignmentStart() == second.getMateAlignmentStart()
                && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(first.getMateReferenceName())
                && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(second.getMateReferenceName())
                && Objects.equals(first.getReferenceName(), second.getMateReferenceName())
                && Objects.equals(first.getMateReferenceName(), second.getReferenceName());
    }

    /**
     * Determine if there is a single primary alignment in a list of alignments.
     * @param records
     * @return NONE, ONE or MORE_THAN_ONE.
     */
    private NumPrimaryAlignmentState tallyPrimaryAlignments(final List<SAMRecord> records) {
        boolean seenPrimary = false;
        for (SAMRecord record : records) {
            if (record != null && !record.isSecondaryOrSupplementary()) {
                if (seenPrimary) return NumPrimaryAlignmentState.MORE_THAN_ONE;
                else seenPrimary = true;
            }
        }
        return seenPrimary ? NumPrimaryAlignmentState.ONE : NumPrimaryAlignmentState.NONE;
    }

    public NumPrimaryAlignmentState tallyPrimaryAlignments(final boolean firstEnd) {
        final List<SAMRecord> records = firstEnd ? firstOfPairOrFragment : secondOfPair;
        return tallyPrimaryAlignments(records);
    }

    List<SAMRecord> getSupplementalFirstOfPairOrFragment() {
        return supplementalFirstOfPairOrFragment;
    }

    List<SAMRecord> getSupplementalSecondOfPair() {
        return supplementalSecondOfPair;
    }
}
