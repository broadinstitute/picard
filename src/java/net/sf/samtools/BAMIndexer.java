/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import net.sf.samtools.util.CloseableIterator;

import java.io.*;

/**
 * Class for both constructing BAM index content and writing it out.
 */
public class BAMIndexer {

    // The number of references (chromosomes) in the BAM file
    private final int numReferences;

    // input bam file
    private final File inputFile;

    // content is built up from the input bam file using this
    private final BAMIndexBuilder indexBuilder;

    // output written as binary, or (for debugging) as text
    private final BAMIndexWriter outputWriter;

    private int currentReference = 0;

    /**
     * @param input       BAM (.bam) file
     * @param output      binary BAM Index (.bai) file
     * @param nReferences Number of references in the input BAM file
     */
    public BAMIndexer(final File input, final File output, final int nReferences, final boolean sortBins) {

        numReferences = nReferences;
        inputFile = input;
        indexBuilder = new BAMIndexBuilder();
        outputWriter = new BinaryBAMIndexWriter(nReferences, output, sortBins);
        outputWriter.writeHeader();
    }

    /**
     * Generates a BAM index file from an input BAM file
     */
    public void createIndex() {

        SAMFileReader reader = new SAMFileReader(inputFile);
        reader.enableFileSource(true);
        CloseableIterator<SAMRecord> alignmentIterator = reader.iterator();
        int totalRecords = 0;
        SAMRecord rec = null;

        // create and write the content
        try {
            while (alignmentIterator.hasNext()) {
                if (++totalRecords % 1000000 == 0) {
                    verbose(totalRecords + " reads processed ...");
                }
                rec = alignmentIterator.next();
                processAlignment(rec);
            }
            alignmentIterator.close();
            finish();

        } catch (Exception e) {
            outputWriter.deleteIndexFile();
            throw new SAMException("Exception creating BAM index on record " + rec, e);
        }
    }

    /**
     * Record any index information for a given BAM record.
     * If this alignment starts a new reference, write out the old reference
     *
     * @param rec The BAM record
     */
    public void processAlignment(final SAMRecord rec) {
        final int reference = rec.getReferenceIndex();
        if (reference != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && reference != currentReference){
            // process any completed references
            while (currentReference < reference) {
                BAMIndexContent content = indexBuilder.processReference(currentReference);
                outputWriter.writeReference(content, currentReference);
                indexBuilder.startNewReference();
                currentReference++;
            }
            currentReference = reference;
        }
        indexBuilder.processAlignment(rec);
    }

    /**
     * After all the alignment records have been processed, finish is called.
     * Note, we can do this processing per-reference instead of per-file if desired
     */
    public void finish() {
        // process any remaining references
        while (currentReference < numReferences) {
            BAMIndexContent content = indexBuilder.processReference(currentReference);
            outputWriter.writeReference(content, currentReference);
            indexBuilder.startNewReference();
            currentReference++;
        }
        long noCoordinateRecords = indexBuilder.finish();
        outputWriter.close(noCoordinateRecords); // writes count of noCoordinateRecords
    }

    /**
     * Deletes old or partial index file
     * Called whenever exceptions occur.
     */
    public void deleteIndex(){
        if (outputWriter != null){
            outputWriter.deleteIndexFile();
        }
    }

    private void verbose(String message) {
        boolean verbose = true;
        if (verbose) {
            System.out.println("BAMIndexer: " + message);
        }
    }
}