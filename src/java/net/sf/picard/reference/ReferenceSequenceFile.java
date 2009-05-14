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

import net.sf.samtools.SAMSequenceDictionary;

/**
 * An interface for working with files of reference sequences regardless of the file format
 * being used.
 *
 * @author Tim Fennell
 */
public interface ReferenceSequenceFile {

    /**
     * Must return a sequence dictionary with at least the following fields completed
     * for each sequence: name, length.
     *
     * @return a list of sequence records representing the sequences in this reference file
     */
    public SAMSequenceDictionary getSequenceDictionary();

    /**
     * Retrieves the next whole sequences from the file.
     * @return a ReferenceSequence or null if at the end of the file
     */
    public ReferenceSequence nextSequence();

    /**
     * @return Reference name, file name, or something other human-readable representation.
     */
    public String toString();
}
