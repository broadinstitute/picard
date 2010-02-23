/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMSequenceRecord;

/**
 * Implementation of ReferenceSequenceMask that indicates that all the loci in the sequence dictionary are of interest.
 * @author alecw at broadinstitute dot org
 */
public class WholeGenomeReferenceSequenceMask implements ReferenceSequenceMask {

    SAMFileHeader header;

    public WholeGenomeReferenceSequenceMask(final SAMFileHeader header) {
        this.header = header;
    }

    /**
     * @return true if the mask is set for the given sequence and position
     */
    public boolean get(final int sequenceIndex, final int position) {
        if (sequenceIndex < 0) {
            throw new IllegalArgumentException("Negative sequence index " + sequenceIndex);
        }
        if (sequenceIndex >= header.getSequenceDictionary().size()) {
            return false;
        }
        final SAMSequenceRecord sequenceRecord = header.getSequence(sequenceIndex);
        return position <= sequenceRecord.getSequenceLength();
    }

    /**
     * @return the next pos on the given sequence >= position that is set, or -1 if there are no more set positions
     */
    public int nextPosition(final int sequenceIndex, final int position) {
        if (get(sequenceIndex, position + 1)) {
            return position + 1;
        } else {
            return -1;
        }
    }

    /**
     * @return Largest sequence index for which there are set bits.
     */
    public int getMaxSequenceIndex() {
        return header.getSequenceDictionary().size() - 1;
    }

    /**
     * @return the largest position on the last sequence index
     */
    public int getMaxPosition() {
        SAMSequenceRecord lastSequenceRecord = header.getSequence(getMaxSequenceIndex());
        return lastSequenceRecord.getSequenceLength();
    }
}
