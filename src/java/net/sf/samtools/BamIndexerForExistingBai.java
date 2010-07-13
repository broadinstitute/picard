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
import java.io.IOException;
import java.util.Collections;
import java.util.List;

/**
 * Class used only for testing. Constructs BAM index content from an existing bai file
 * and writes it out (as binary bai file or textual bai.txt file)
 */
public class BamIndexerForExistingBai {

    // input either built from bam file, or (for debugging) from existing bai file
    private final File inputFile;

    /**
     * Constructor
     *
     * @param input       BAM Index (.bai) file
     */
    public BamIndexerForExistingBai(final File input) {

        this.inputFile = input;
    }

    /**
     * Generates a BAM index file, either textual or binary, sorted or not, from an input BAI file
     *
     * @param output      BAM Index (.bai) file (or bai.txt file when text)
     * @param textOutput  Whether to create text output or binary
     * @param sortBins     Whether to sort the bins in the output
     */
    public void createIndex(final File output, final boolean textOutput, final boolean sortBins) {

        // content is from an existing bai file.

        // final DiskBasedBAMFileIndex existingIndex = new DiskBasedBAMFileIndex(inputFile); // doesn't work todo
        final CachingBAMFileIndex existingIndex = new CachingBAMFileIndex(inputFile); // also works
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
            outputWriter.close(existingIndex.getNoCoordinateCount());
            existingIndex.close();

        } catch (Exception e) {
            outputWriter.close(null); // let's keep the partial file for now  todo
 //            outputWriter.deleteIndexFile();
            throw new SAMException("Exception creating BAM index", e);
        }
    }

    /**
     * Prints meta-data statistics from BAM index (.bai) file
     * Statistics include count of aligned and unaligned reads for each reference sequence
     * and a count of all records with no start coordinate
     */
    public void indexStats() {
        try {
            final BAMFileReader bam = new BAMFileReader(inputFile, null, false, SAMFileReader.ValidationStringency.SILENT);
            // bam.enableIndexCaching(true);  // to test CachingBAMFileIndex

            if (!bam.hasIndex()) {
                throw new SAMException("No index for bam file " + inputFile);
            }

            AbstractBAMFileIndex index = (AbstractBAMFileIndex) bam.getIndex();
            index.open();
            // read through all the bins of every reference.
            int nRefs = index.getNumberOfReferences();
            for (int i = 0; i < nRefs; i++) {
                BAMIndexContent content = index.query(i, 0, -1); // todo: it would be faster just to skip to the last bin

                final SAMSequenceRecord seq = bam.getFileHeader().getSequence(i);
                if (seq == null) continue;
                final String sequenceName = seq.getSequenceName();
                final int sequenceLength = seq.getSequenceLength();
                System.out.print(sequenceName + ' ' + "length =" + sequenceLength);

                if (content == null || content.getBins() == null || content.getBins().isEmpty()) {
                    System.out.println();
                    continue;
                }
                Bin bin = Collections.max(content.getBins()); // might be expensive?

                boolean firstChunk = true;
                if (bin.equals(new Bin(i, BAMIndex.MAX_BINS))) {
                    List<Chunk> chunkList = bin.getChunkList();
                    for (Chunk c : chunkList) {
                        long start = c.getChunkStart();
                        long end = c.getChunkEnd();
                        if (firstChunk){
                            // samtools idxstats doesn't print this, so we won't either
                            // System.out.print(sequenceName + ' ' + "Start=" + start + "    End=" + end);
                            firstChunk = false;
                        } else {
                            firstChunk = true;
                            System.out.println("    Aligned= " + start + "   Unaligned= " + end);
                        }
                    }
                } else {
                    // no meta data for this index
                }
            }

            System.out.println("NoCoordinateCount= " + index.getNoCoordinateCount());

        } catch (IOException e) {
            throw new SAMException("Exception in getting index statistics", e);
        }

    }
}