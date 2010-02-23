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

/**
 * Interface for specifying loci of interest for genotype calling and other operations.
 * It is a requirement that the sequences be probed in ascending order.
 *
 * @author alecw at broadinstitute dot oh are gee
 */
public interface ReferenceSequenceMask {

    /**
     * It is required that sequenceIndex is >= any previous sequenceIndex passed to this class.
     * @return true if the mask is set for the given sequence and position
     */
    boolean get(int sequenceIndex, int position);

    /**
     * It is required that sequenceIndex is >= any previous sequenceIndex passed to this class.
     * @return the next pos on the given sequence >= position that is set, or -1 if there are no more set positions
     */
    int nextPosition(int sequenceIndex, int position);

    /**
     * @return Largest sequence index for which there are set bits.
     */
    int getMaxSequenceIndex();

    /**
     * @return the largest position on the last sequence index
     */
    int getMaxPosition();
}
