/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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
package net.sf.picard.reference;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMSequenceDictionary;

import java.io.File;

/**
 * Manages a ReferenceSequenceFile.  Loads the requested sequence, ensuring that
 * access is in order, and confirming that sequence name and index agree.
 *
 * @author alecw@broadinstitute.org
 */
public class ReferenceSequenceFileWalker {
    private final ReferenceSequenceFile referenceSequenceFile;
    private ReferenceSequence referenceSequence = null;

    public ReferenceSequenceFileWalker(final ReferenceSequenceFile referenceSequenceFile) {
        this.referenceSequenceFile = referenceSequenceFile;
    }

    public ReferenceSequenceFileWalker(final File file) {
        this(ReferenceSequenceFileFactory.getReferenceSequenceFile(file));
    }

    /**
     * Ensure that the requested sequence is loaded.  Throws an exception if out-of-order
     * request is made, or if there is a mismatch between the requested name and the name
     * found in the ReferenceSequenceFile
     */
    public ReferenceSequence get(final int sequenceIndex, final String sequenceName, final int length) {
        // Has the side-effect of setting referenceSequence member
        get(sequenceIndex);
        if (!referenceSequence.getName().equals(sequenceName)) {
            // Sanity check the sequence names against the sequence dictionary while scanning through.
            throw new PicardException("Sequence name mismatch at sequence index (" + referenceSequence.getContigIndex() +
                    ", " + referenceSequence.getName() + ") != " + sequenceName);
        }
        if (referenceSequence.getBases().length != length) {
            throw new PicardException("Sequence length mismatch for (" + sequenceIndex + ", " + sequenceName +
            ").  expected " + length + " but found " + referenceSequence.getBases().length);
        }
        return referenceSequence;
    }

    /**
     * Get reference sequence without validating name or length.  This is OK if the entire sequence
     * dictionary was validated before reading sequences.
     */
    public ReferenceSequence get(final int sequenceIndex) {
        if (referenceSequence != null && referenceSequence.getContigIndex() == sequenceIndex) {
            return referenceSequence;
        }
        if (referenceSequence != null && referenceSequence.getContigIndex() > sequenceIndex) {
            throw new PicardException("Requesting earlier reference sequence: " + sequenceIndex + " < " +
            referenceSequence.getContigIndex());
        }
        referenceSequence = null;
        for(referenceSequence = referenceSequenceFile.nextSequence();
                referenceSequence != null && referenceSequence.getContigIndex() < sequenceIndex;
                referenceSequence = referenceSequenceFile.nextSequence()) {
        }
        if (referenceSequence == null || referenceSequence.getContigIndex() != sequenceIndex) {
            throw new PicardException("Reference sequence (" + sequenceIndex +
            ") not found in " + referenceSequenceFile.toString());
        }
        return referenceSequence;
    }

    public SAMSequenceDictionary getSequenceDictionary() {
        return referenceSequenceFile.getSequenceDictionary();
    }
}
