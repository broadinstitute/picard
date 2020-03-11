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
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.CoordMath;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * For a paired-end aligner that aligns each end independently, select the pair of alignments that result
 * in the largest insert size.  If such a pair of alignments cannot be found, either because one end is not aligned,
 * or because all alignment pairs are chimeric, then select the best MAPQ for each end independently.
 *
 * The primary alignments are then correlated so that their mate info points to each
 * other, but all non-primary alignments are uncorrelated.
 */
public class MostDistantPrimaryAlignmentSelectionStrategy implements PrimaryAlignmentSelectionStrategy {
    // Give random number generator a seed so results are repeatable.  Used to pick a primary alignment from
    // multiple alignments with equal mapping quality.
    private final Random random = new Random(1);

    @Override
    public void pickPrimaryAlignment(final HitsForInsert hitsForInsert) {
        final BestEndAlignmentsAccumulator firstEndBest = new BestEndAlignmentsAccumulator();
        final BestEndAlignmentsAccumulator secondEndBest = new BestEndAlignmentsAccumulator();
        final CollectionUtil.MultiMap<Integer, SAMRecord> firstEndBySequence =
                new CollectionUtil.MultiMap<Integer, SAMRecord>();
        final BestPairAlignmentsAccumulator pairBest = new BestPairAlignmentsAccumulator();

        for (final SAMRecord rec : hitsForInsert.firstOfPairOrFragment) {
            if (rec.getReadUnmappedFlag()) throw new IllegalStateException();
            firstEndBest.considerBest(rec);
            firstEndBySequence.append(rec.getReferenceIndex(), rec);
        }

        for (final SAMRecord secondEnd: hitsForInsert.secondOfPair) {
            if (secondEnd.getReadUnmappedFlag()) throw new IllegalStateException();
            secondEndBest.considerBest(secondEnd);
            final Collection<SAMRecord> firstEnds = firstEndBySequence.get(secondEnd.getReferenceIndex());
            if (firstEnds != null) {
                for (final SAMRecord firstEnd : firstEnds) {
                    pairBest.considerBest(firstEnd, secondEnd);
                }
            }
        }

        final SAMRecord bestFirstEnd;
        final SAMRecord bestSecondEnd;
        if (pairBest.hasBest()) {
            final Map.Entry<SAMRecord, SAMRecord> pairEntry = pickRandomlyFromList(pairBest.bestAlignmentPairs);
            bestFirstEnd = pairEntry.getKey();
            bestSecondEnd = pairEntry.getValue();
        } else {
            if (firstEndBest.hasBest()) {
                bestFirstEnd = pickRandomlyFromList(firstEndBest.bestAlignments);
            } else {
                bestFirstEnd = null;
            }
            if (secondEndBest.hasBest()) {
                bestSecondEnd = pickRandomlyFromList(secondEndBest.bestAlignments);
            } else {
                bestSecondEnd = null;
            }
        }

        if (hitsForInsert.firstOfPairOrFragment.isEmpty() != (bestFirstEnd == null)) {
            throw new IllegalStateException("Should not happen");
        }
        if (hitsForInsert.secondOfPair.isEmpty() != (bestSecondEnd == null)) {
            throw new IllegalStateException("Should not happen");
        }
        if (bestFirstEnd != null) {
            moveToHead(hitsForInsert.firstOfPairOrFragment, bestFirstEnd);
        }
        if (bestSecondEnd != null) {
            moveToHead(hitsForInsert.secondOfPair, bestSecondEnd);
        }
        hitsForInsert.setPrimaryAlignment(0);

        // For non-primary alignments, de-correlate them so that the mate fields don't point at some
        // arbitrary alignment for the other end.

        // No non-primary alignments for one of the ends so nothing to do.
        if (hitsForInsert.firstOfPairOrFragment.size() <= 1 || hitsForInsert.secondOfPair.size() <= 1) return;
        final int amountToSlide = hitsForInsert.firstOfPairOrFragment.size() - 1;
        for (int i = 0; i < amountToSlide; ++i) {
            hitsForInsert.secondOfPair.add(1, null);
        }


    }

    private <T> T pickRandomlyFromList(final List<T> list) {
        return list.get(random.nextInt(list.size()));
    }

    // Uses reference equality, not .equals()
    private void moveToHead(final List<SAMRecord> list, final SAMRecord rec) {
        if (list.get(0) == rec) return;
        for (int i = 1; i < list.size(); ++i) {
            if (list.get(i) == rec) {
                list.remove(i);
                list.add(0, rec);
                return;
            }
        }
        throw new IllegalStateException("Should not be reached");
    }

    private static class BestEndAlignmentsAccumulator {
        public int bestMapq = -1;
        public List<SAMRecord> bestAlignments = new ArrayList<SAMRecord>();

        public void considerBest(final SAMRecord rec) {
            if (bestMapq == -1) {
                bestMapq = rec.getMappingQuality();
                bestAlignments.add(rec);
            } else {
                final int cmp = SAMUtils.compareMapqs(bestMapq, rec.getMappingQuality());
                if (cmp < 0) {
                    bestMapq = rec.getMappingQuality();
                    bestAlignments.clear();
                    bestAlignments.add(rec);
                } else if (cmp == 0) {
                    bestAlignments.add(rec);
                }
            }
        }

        public boolean hasBest() {
            return bestMapq != -1;
        }
    }

    private static class BestPairAlignmentsAccumulator {
        public int bestDistance = -1;
        public int bestPairMapq = -1;
        public List<Map.Entry<SAMRecord, SAMRecord>> bestAlignmentPairs =
                new ArrayList<Map.Entry<SAMRecord, SAMRecord>>();

        public void considerBest(final SAMRecord firstEnd, final SAMRecord secondEnd) {
            final int thisPairMapq = SAMUtils.combineMapqs(firstEnd.getMappingQuality(), secondEnd.getMappingQuality());
            final int thisDistance = CoordMath.getLength(Math.min(firstEnd.getAlignmentStart(), secondEnd.getAlignmentStart()),
                    Math.max(firstEnd.getAlignmentEnd(), secondEnd.getAlignmentEnd()));
            if (thisDistance > bestDistance || (thisDistance == bestDistance && thisPairMapq > bestPairMapq)) {
                bestDistance = thisDistance;
                bestPairMapq = thisPairMapq;
                bestAlignmentPairs.clear();
                bestAlignmentPairs.add(new AbstractMap.SimpleEntry<SAMRecord, SAMRecord>(firstEnd, secondEnd));
            } else if (thisDistance == bestDistance && thisPairMapq == bestPairMapq) {
                bestAlignmentPairs.add(new AbstractMap.SimpleEntry<SAMRecord, SAMRecord>(firstEnd, secondEnd));
            }
        }

        public boolean hasBest() {
            return bestDistance != -1;
        }
    }
}
