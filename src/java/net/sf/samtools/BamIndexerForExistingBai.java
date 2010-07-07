/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sub-license, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools;

import java.io.File;

/**
 * Class used only for testing. Constructs BAM index content from an existing bai file
 * and writes it out (as binary bai file or textual bai.txt file)
 */
public class BamIndexerForExistingBai {

    // input either built from bam file, or (for debugging) from existing bai file
    private final File inputFile;

    // whether to sort the bins in the index.
    private final boolean sortBins;

    /**
     * Constructor
     *
     * @param input       BAM Index (.bai) file
     * @param sortBins    Whether to sort the bins in the output
     */
    public BamIndexerForExistingBai(final File input, final boolean sortBins) {

        this.inputFile = input;
        this.sortBins= sortBins;
    }

    /**
     * Generates a BAM index file, either textual or binary, sorted or not, from an input BAI file
     *
     *  @param output      BAM Index (.bai) file (or bai.txt file when text)
     *  @param textOutput  Whether to create text output or binary
     */
    public void createIndex(final File output, final boolean textOutput) {

        // content is (for debugging) from an existing bai file.
        final CachingBAMFileIndex existingIndex = new CachingBAMFileIndex(inputFile);
        final int n_ref = existingIndex.getNumberOfReferences();
        final BAMIndexWriter outputWriter;
        if (textOutput){
            outputWriter = new TextualBAMIndexWriter(n_ref, output, sortBins);
        } else {
            outputWriter = BAMIndexWriterFactory.makeBAMIndexWriter(n_ref, output, sortBins, inputFile.length());
        }
        outputWriter.writeHeader();

        // write the content one reference at a time
        try {
            for (int i = 0; i < n_ref; i++) {
                outputWriter.writeReference(existingIndex.getQueryResults(i), i);
            }
            outputWriter.close();

        } catch (Exception e) {
            outputWriter.deleteIndexFile();
            throw new SAMException("Exception creating BAM index", e);
        }
    }
}