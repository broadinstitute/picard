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

import net.sf.samtools.util.CloseableIterator;

import java.io.*;

/**
 * Class for constructing BAM index content from a BAM file
 * and writing it out (as binary bai file or textual bai.txt file for debugging)
 */
public class BAMIndexer {

    /**
     * The number of references (chromosomes) in the BAM file
     */
    private final int numReferences;

    // input bam file
    private final File inputFile;

    // content is built up from the input bam file using this
    private final BAMIndexBuilder indexBuilder;

    // output written as binary, or (for debugging) as text
    private final BAMIndexWriter outputWriter;

    private int currentReference = -1;

    /**
     * Constructor
     *
     * @param input       BAM (.bam) file
     * @param output      binary BAM Index (.bai) file
     * @param nReferences Number of references in the input BAM file, or -1 if unknown
     */
    public BAMIndexer(final File input, final File output, final int nReferences, final boolean sortBins) {

        numReferences = nReferences;
        inputFile = input;
        indexBuilder = new BAMIndexBuilder();
        outputWriter = BAMIndexWriterFactory.makeBAMIndexWriter(nReferences, output, sortBins, input.length());
        outputWriter.writeHeader();
    }

    /**
     * Constructor that allows specifying sorted text output
     * * Todo - hide this.  Only called by BAMIndexWriterText and BuildBamIndex(for shorthand to get text output)
     *
     * @param input       BAM (.bam) file
     * @param output      BAM Index (.bai) file (or bai.txt file when text)
     * @param nReferences Number of references in the input BAM file, or -1 if unknown
     * @param sortBins    Whether to sort the bins in the output
     * @param textOutput  Whether to create text output or binary
     */
    public BAMIndexer(final File input, final File output, final int nReferences, final boolean sortBins, final boolean textOutput) {

        numReferences = nReferences;
        inputFile = input;
        indexBuilder = new BAMIndexBuilder();
        outputWriter = textOutput ? new TextualBAMIndexWriter(nReferences, output, sortBins)
                : BAMIndexWriterFactory.makeBAMIndexWriter(nReferences, output, sortBins, input.length());
        outputWriter.writeHeader();
    }

    /**
     * Generates a BAM index file from an input BAM file
     */
    public void createIndex() {

        int noCoordinateRecords = 0;
        SAMFileReader reader = new SAMFileReader(inputFile);
        reader.enableFileSource(true);
        CloseableIterator<SAMRecord> alignmentIterator = reader.iterator();
        int totalRecords = 0;

        // create and write the content
        try {
            while (alignmentIterator.hasNext()) {
                if (++totalRecords % 1000000 == 0) {
                    verbose(totalRecords + " reads processed ...");
                }
                final SAMRecord rec = alignmentIterator.next();
                if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                    noCoordinateRecords++;
                    continue;  // do nothing for un-aligned records, except count
                    // this count should match           noCoordinateRecord in BAMIndexBuilder
                }

                processAlignment(rec);
            }
            alignmentIterator.close();
            finish();
            verbose("Records without coordinates = " + noCoordinateRecords);

        } catch (Exception e) {
            outputWriter.deleteIndexFile();
            throw new SAMException("Exception creating BAM index", e);
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
            // process any skipped references
            for (int i = currentReference + 1; i < reference; i++) {
                // System.err.println("### skipping reference " + i);
                BAMIndexContent skippedContent = indexBuilder.processReference(i);
                outputWriter.writeReference(skippedContent, i);
                indexBuilder.startNewReference();
            }
            currentReference = reference;
        }
        final BAMIndexContent finishedContent = indexBuilder.processAlignment(rec);
        if (finishedContent != null) {
            outputWriter.writeReference(finishedContent, finishedContent.getReferenceSequence());
        }
    }

    /**
     * After all the alignment records have been processed, finish is called.
     * Note, we can do this processing per-reference instead of per-file if desired
     */
    public void finish() {
        // process any skipped references
        for (int i = currentReference; i < numReferences; i++) {
            // System.err.println("### ending reference " + i);
            BAMIndexContent skippedContent = indexBuilder.processReference(i);
            outputWriter.writeReference(skippedContent, i);
            indexBuilder.startNewReference();
        }
        // todo add meta information of noCoordinateRecordCount
        outputWriter.close();
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