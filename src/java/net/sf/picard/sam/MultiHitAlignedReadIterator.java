/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.filter.FilteringIterator;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.util.Log;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordQueryNameComparator;
import net.sf.samtools.SAMTag;
import net.sf.samtools.util.CloseableIterator;

import java.util.*;


/**
 * Iterate over queryname-sorted SAM, and return each group of reads with the same queryname.  Unmapped reads
 * are filtered out.  If there are multiple hits for the same read, and the first and second ends need to be correlated,
 * then they are sorted by hit index.
 * A set of hits for a single query may then be filtered with a caller-supplied filter, which will remove any
 * alignments that do not pass the filter.  If the primary alignment is removed, the best-mapping secondary alignment
 * or alignment pair will be marked as primary.
 *
 *
 * @throws IllegalStateException if the input is not queryname-sorted.
 */
class MultiHitAlignedReadIterator implements CloseableIterator<MultiHitAlignedReadIterator.HitsForInsert> {
    private static final Log log = Log.getInstance(MultiHitAlignedReadIterator.class);

    private final PeekableIterator<SAMRecord> peekIterator;
    private final HitIndexComparator comparator = new HitIndexComparator();
    private final SAMRecordQueryNameComparator queryNameComparator = new SAMRecordQueryNameComparator();
    private final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy;

    // Give random number generator a seed so results are repeatable.  Used to pick a primary alignment from
    // multiple alignments with equal mapping quality.
    private final Random random = new Random(1);

    /**
     *
     * @param querynameOrderIterator
     * @param primaryAlignmentSelectionStrategy Algorithm for selecting primary alignment when it is not clear from
     *                                          the input what should be primary.
     */
    MultiHitAlignedReadIterator(final CloseableIterator<SAMRecord> querynameOrderIterator,
                                final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy) {
        this.primaryAlignmentSelectionStrategy = primaryAlignmentSelectionStrategy;
        peekIterator = new PeekableIterator<SAMRecord>(new FilteringIterator(querynameOrderIterator,
                new SamRecordFilter() {
                    public boolean filterOut(final SAMRecord record) {
                        return record.getReadUnmappedFlag();
                    }
                    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
                        return (first.getReadUnmappedFlag() && second.getReadUnmappedFlag());
                    }
                }));
    }

    public void close() {
        peekIterator.close();
    }

    public boolean hasNext() {
        return peekIterator.hasNext();
    }

    /**
     * @throws IllegalStateException if the input is not queryname-sorted.
     */
    public HitsForInsert next() {
        final String readName = peekIterator.peek().getReadName();
        final HitsForInsert hits = new HitsForInsert();

        Boolean isPaired = null;

        // Accumulate the alignments matching readName.
        while (peekIterator.hasNext() && peekIterator.peek().getReadName().equals(readName)) {
            final SAMRecord rec = peekIterator.next();
            // It is critical to do this here, because SamAlignmentMerger uses this exception to determine
            // if the aligned input needs to be sorted.
            if (peekIterator.hasNext() && queryNameComparator.fileOrderCompare(rec, peekIterator.peek()) > 0) {
                throw new IllegalStateException("Underlying iterator is not queryname sorted: " +
                rec + " > " + peekIterator.peek());
            }
            if (isPaired == null) {
                isPaired = rec.getReadPairedFlag();
            } else if (isPaired != rec.getReadPairedFlag()) {
                throw new PicardException("Got a mix of paired and unpaired alignments for read " + readName);
            }
            if (!rec.getReadPairedFlag() || rec.getFirstOfPairFlag()) {
                hits.firstOfPairOrFragment.add(rec);
            } else if (rec.getSecondOfPairFlag()) {
                hits.secondOfPair.add(rec);
            } else throw new PicardException("Read is marked as pair but neither first or second: " + readName);
        }

        // For paired-end alignments that need to be correlated, sort by hit index.
        if (!hits.firstOfPairOrFragment.isEmpty() && !hits.secondOfPair.isEmpty()) {
            if (hits.firstOfPairOrFragment.size() != hits.secondOfPair.size()) {
                throw new PicardException("Unequal number of first and second of pair hits for read: " + readName);
            }
            if (hits.firstOfPairOrFragment.size() > 1) {
                Collections.sort(hits.firstOfPairOrFragment, comparator);
                Collections.sort(hits.secondOfPair, comparator);
                for (int i = 0; i < hits.firstOfPairOrFragment.size(); ++i) {
                    if (hits.firstOfPairOrFragment.get(i).getIntegerAttribute(SAMTag.HI.name()) == null) {
                        throw new PicardException("Missing HI tag on multi-hit alignment for read: " + readName);
                    }
                    if (!hits.firstOfPairOrFragment.get(i).getIntegerAttribute(SAMTag.HI.name()).
                            equals(hits.secondOfPair.get(i).getIntegerAttribute(SAMTag.HI.name()))) {
                        throw new PicardException("Hit index mismatch for multi-hit alignment for read: " + readName);
                    }
                }
            }
        }
        // See if primary alignment is not already unambiguously determined.
        final NumPrimaryAlignmentState firstEndAlignmentState = hits.tallyPrimaryAlignments(hits.firstOfPairOrFragment);
        final NumPrimaryAlignmentState secondEndAlignmentState = hits.tallyPrimaryAlignments(hits.secondOfPair);
        if ((firstEndAlignmentState == NumPrimaryAlignmentState.NONE && secondEndAlignmentState == NumPrimaryAlignmentState.NONE) ||
                firstEndAlignmentState == NumPrimaryAlignmentState.MORE_THAN_ONE ||
                secondEndAlignmentState == NumPrimaryAlignmentState.MORE_THAN_ONE) {
            // Need to use selected strategy for picking primary.
            primaryAlignmentSelectionStrategy.pickPrimaryAlignment(hits);
        } else {
            // Primary is already determined.  Just check that it makes sense, i.e. if specified for both ends,
            // that the specification is the same for both ends.
            final int firstEndPrimaryAlignmentIndex = hits.findPrimaryAlignment(hits.firstOfPairOrFragment);
            final int secondEndPrimaryAlignmentIndex = hits.findPrimaryAlignment(hits.secondOfPair);
            if (firstEndPrimaryAlignmentIndex != -1 && secondEndPrimaryAlignmentIndex != -1 &&
                    firstEndPrimaryAlignmentIndex != secondEndPrimaryAlignmentIndex) {
                throw new PicardException("Mismatched primary alignments for paired read " + hits.getReadName());
            }
        }
        return hits;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /**
     * Holds all the hits (alignments) for a read or read pair.  For a read pair, either there are N hits
     * for one end and 0 hits for the other end, or there are N hits for both ends.  If there
     * are N hits for both ends, then the nth hit for the first end is correlated to the nth of
     * the second end.
     */
    public class HitsForInsert {
        private final List<SAMRecord> firstOfPairOrFragment = new ArrayList<SAMRecord>();
        private final List<SAMRecord> secondOfPair = new ArrayList<SAMRecord>();

        /**
         * @throws if numHits() == 0
         */
        public String getReadName() {
            return getRepresentativeRead().getReadName();
        }

        /**
         * @throws if numHits() == 0
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

        /**
         * @return Returns the ith hit for the first end, or null if the first end is not aligned.
         */
        public SAMRecord getFirstOfPair(final int i) {
            if (firstOfPairOrFragment.isEmpty()) {
                return null;
            } else {
                return firstOfPairOrFragment.get(i);
            }
        }

        /**
         * @return The ith hit for a un-paired read.  Never returns null.
         * Do not call if paired read.
         */
        public SAMRecord getFragment(final int i) {
            final SAMRecord samRecord = firstOfPairOrFragment.get(i);
            if (samRecord.getReadPairedFlag()) throw new UnsupportedOperationException("getFragment called for paired read");
            return samRecord;
        }

        /**
         * @return Returns the ith hit for the second end, or null if the second end is not aligned.
         */
        public SAMRecord getSecondOfPair(final int i) {
            if (secondOfPair.isEmpty()) {
                return null;
            } else {
                return secondOfPair.get(i);
            }
        }

        /**
         * Set all alignments to not primary, except for the one specified by the argument.  If paired, and set the
         * alignment for both ends if there is an alignment for both ends, otherwise just for the end for which
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
         * Remove reads that filter rejects, then squeeze down lists.  If the primary alignment is filtered out,
         * another alignment will be selected as primary, by looking at best MAPQ.
         * Note that after filtering, numHits may == 0.
         * Hit indices may be revised if alignments get filtered out.
         */
        public void filterReads(final SamRecordFilter filter, final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy) {
            // Set to null all reads that filter rejects.
            filterReads(filter, firstOfPairOrFragment);
            filterReads(filter, secondOfPair);

            if (!firstOfPairOrFragment.isEmpty() && !secondOfPair.isEmpty()) {
                if (firstOfPairOrFragment.size() != secondOfPair.size()) {
                    throw new IllegalStateException("Both ends of pair have hits, but not the same number of hits for each end");
                }
                // If both ends of an alignment have been set to null, squeeze them out of the list.
                // If one has been set to null, clear HI tag on the other end.
                for (int i = firstOfPairOrFragment.size() - 1; i >= 0; --i) {
                    if (firstOfPairOrFragment.get(i) == null) {
                        if (secondOfPair.get(i) == null) {
                            // Both ends null, squeeze the list
                            firstOfPairOrFragment.remove(i);
                            secondOfPair.remove(i);
                        } else {
                            secondOfPair.get(i).setAttribute(SAMTag.HI.name(), null);
                        }
                    } else if (secondOfPair.get(i) == null) {
                        firstOfPairOrFragment.get(i).setAttribute(SAMTag.HI.name(), null);
                    }
                }
                if (areAllElementsNull(firstOfPairOrFragment)) firstOfPairOrFragment.clear();
                if (areAllElementsNull(secondOfPair)) secondOfPair.clear();
            } else  {
                // One of these is already empty, but code is simpler this way.
                squeezeList(firstOfPairOrFragment);
                squeezeList(secondOfPair);
            }

            // If all hits have been squeezed out, nothing else to be done.
            if (numHits() == 0) return;

            // reassign HI tag if necessary
            if (numHits() == 1) {
                // Clear HI if only one alignment
                if (getFirstOfPair(0) != null) getFirstOfPair(0).setAttribute(SAMTag.HI.name(), null);
                if (getSecondOfPair(0) != null) getSecondOfPair(0).setAttribute(SAMTag.HI.name(), null);
            } else {
                for (int i = 0; i < numHits(); ++i) {
                    if (getFirstOfPair(i) != null) getFirstOfPair(i).setAttribute(SAMTag.HI.name(), i);
                    if (getSecondOfPair(i) != null) getSecondOfPair(i).setAttribute(SAMTag.HI.name(), i);
                }
            }

            // See if the primary alignment has been squeezed out.
            final int firstEndPrimaryAlignmentIndex = findPrimaryAlignment(firstOfPairOrFragment);
            final int secondEndPrimaryAlignmentIndex = findPrimaryAlignment(secondOfPair);

            // Sanity check
            if (firstEndPrimaryAlignmentIndex != -1) {
                if (secondEndPrimaryAlignmentIndex == -1) {
                    if (getSecondOfPair(firstEndPrimaryAlignmentIndex) != null) {
                        throw new IllegalStateException("Mismatched primary alignments for read " + getReadName());
                    }
                } else if (firstEndPrimaryAlignmentIndex != secondEndPrimaryAlignmentIndex) {
                    throw new IllegalStateException("Mismatched primary alignments for read " + getReadName());
                }
            } else if (secondEndPrimaryAlignmentIndex != -1) {
                if (getFirstOfPair(secondEndPrimaryAlignmentIndex) != null) {
                    throw new IllegalStateException("Mismatched primary alignments for read " + getReadName());
                }
            }
            if (firstEndPrimaryAlignmentIndex == -1 && secondEndPrimaryAlignmentIndex == -1) {
                primaryAlignmentSelectionStrategy.pickPrimaryAlignment(this);
            }
        }


        private void filterReads(final SamRecordFilter filter, final List<SAMRecord> records) {
            for (int i = 0; i < records.size(); ++i) {
                if (filter.filterOut(records.get(i))) {
                    records.set(i, null);
                }
            }
        }

        private void squeezeList(final List<SAMRecord> records) {
            for (int i = records.size() - 1; i >= 0; --i) {
                if (records.get(i) == null) records.remove(i);
            }
        }

        private boolean areAllElementsNull(final List<SAMRecord> records) {
            for (final SAMRecord rec: records) {
                if (rec != null) return false;
            }
            return true;
        }

        /**
         * Determine if there is a single primary alignment in a list of alignments.
         * @param records
         * @return NONE, ONE or MORE_THAN_ONE.
         */
        private NumPrimaryAlignmentState tallyPrimaryAlignments(final List<SAMRecord> records) {
            boolean seenPrimary = false;
            for (int i = 0; i < records.size(); ++i) {
                if (records.get(i) != null && !records.get(i).getNotPrimaryAlignmentFlag()) {
                    if (seenPrimary) return NumPrimaryAlignmentState.MORE_THAN_ONE;
                    else seenPrimary = true;
                }
            }
            if (seenPrimary) return NumPrimaryAlignmentState.ONE;
            else return NumPrimaryAlignmentState.NONE;
        }

        private int findPrimaryAlignment(final List<SAMRecord> records) {
            int indexOfPrimaryAlignment = -1;
            for (int i = 0; i < records.size(); ++i) {
                if (records.get(i) != null && !records.get(i).getNotPrimaryAlignmentFlag()) {
                    if (indexOfPrimaryAlignment != -1) {
                        throw new IllegalStateException("Multiple primary alignments found for read " + getReadName());
                    }
                    indexOfPrimaryAlignment = i;
                }
            }
            return indexOfPrimaryAlignment;
        }
    }

    private enum NumPrimaryAlignmentState {
        NONE, ONE, MORE_THAN_ONE;
    }

    private static class HitIndexComparator implements Comparator<SAMRecord> {
        public int compare(final SAMRecord rec1, final SAMRecord rec2) {
            final Integer hi1 = rec1.getIntegerAttribute(SAMTag.HI.name());
            final Integer hi2 = rec2.getIntegerAttribute(SAMTag.HI.name());
            if (hi1 == null || hi2 == null) throw new PicardException("HI tag missing for multi-hit paired alignments for read: " +
            rec1.getReadName());
            return hi1.compareTo(hi2);
        }
    }
}
