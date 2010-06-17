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
package net.sf.samtools;

import java.io.File;

/**
 * Class for writing BAI (BAM index files) as text or binary.
 * Requires an existing index file. Used for testing only.
 */
public class BAMIndexTextWriter extends CachingBAMFileIndex { // *not* DiskBasedBAMFileIndex 

    /**
     * The number of references (chromosomes) in the BAI file
     * The name n_ref corresponds to the name in the SAM Format Specification Document
     */
    public final int n_ref;

    private final File OUTPUT;  // todo could use mFile in super, though currently private

    /**
     * Constructor
     * @param INPUT     A BAM Index File, .bai
     * @param OUTPUT    Reformatted BAM Index File, .bai.txt, or (sorted) binary .bai file
     */
    public BAMIndexTextWriter(final File INPUT, final File OUTPUT) {
        super(INPUT);
        this.OUTPUT = OUTPUT;
        if (OUTPUT.exists()) {
            OUTPUT.delete();
        }
        n_ref = getNumberOfReferences();
    }

    public void writeText(final boolean sortBins) throws Exception {
        writeText(n_ref, OUTPUT, sortBins);
    }

    public void writeBinary(final boolean sortBins, final long bamFileSize) throws Exception {
        writeBinary(n_ref, OUTPUT, sortBins, bamFileSize);
    }
}