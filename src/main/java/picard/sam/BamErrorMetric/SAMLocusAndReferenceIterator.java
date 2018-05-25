/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
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
package picard.sam.BamErrorMetric;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.IterableOnceIterator;
import htsjdk.samtools.util.SamLocusIterator;

import java.util.List;

import static htsjdk.samtools.util.SamLocusIterator.*;

/**
 * Iterator that traverses a SAM File and a ReferenceFile, accumulating information on a per-locus basis.
 * Only, whatever loci are covered by the input reads are returned.
 * Duplicate reads and non-primary alignments are filtered out.
 * Iterator element holds both pileup (in the form of a LocusInfo object) and the reference base
 *
 * @author Yossi Farjoun
 */
public class SAMLocusAndReferenceIterator extends IterableOnceIterator<SAMLocusAndReferenceIterator.SAMLocusAndReference> {

    private final ReferenceSequenceFileWalker referenceSequenceFileWalker;
    private final SamLocusIterator locusIterator;

    public SAMLocusAndReferenceIterator(final ReferenceSequenceFileWalker referenceFile, final SamLocusIterator locusIterator) {
        this.referenceSequenceFileWalker = referenceFile;
        this.locusIterator = locusIterator;
    }

    @Override
    public boolean hasNext() {
        return locusIterator.hasNext();
    }

    @Override
    public SAMLocusAndReference next() {
        final LocusInfo locus = locusIterator.next();
        final ReferenceSequence referenceSequence = referenceSequenceFileWalker.get(locus.getSequenceIndex(), locus.getSequenceName(),
                locus.getSequenceLength());

        //position is 1-based...arrays are 0-based!
        return new SAMLocusAndReference(locus, referenceSequence.getBases()[locus.getPosition() - 1]);
    }

    public static class SAMLocusAndReference {
        public final LocusInfo locus;
        public final byte referenceBase;

        public SAMLocusAndReference(final LocusInfo locus, final byte referenceBase) {
            this.locus = locus;
            this.referenceBase = referenceBase;
        }

        public List<RecordAndOffset> getRecordAndOffsets() {
            return locus.getRecordAndOffsets();
        }
    }
}
