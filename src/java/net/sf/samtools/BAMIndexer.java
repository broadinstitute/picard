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

import net.sf.samtools.util.BlockCompressedFilePointerUtil;

import java.io.*;
import java.util.Arrays;
import java.util.List;

import static net.sf.samtools.AbstractBAMFileIndex.MAX_BINS;

/**
 * Class for both constructing BAM index content and writing it out.
 * There are two usage patterns:
 * 1) Building a bam index from an existing bam file
 * 2) Building a bam index while building the bam file
 * In both cases, processAlignment is called for each alignment record and
 * finish() is called at the end.
 */
public class BAMIndexer {

    // The number of references (chromosomes) in the BAM file
    private final int numReferences;

    // output written as binary, or (for debugging) as text
    private final BAMIndexWriter outputWriter;

    private int currentReference = 0;

    // content is built up from the input bam file using this
    private final BAMIndexBuilder indexBuilder;

    /**
     * @param output     binary BAM Index (.bai) file
     * @param fileHeader header for the corresponding bam file
     */
    public BAMIndexer(final File output, SAMFileHeader fileHeader) {

        numReferences = fileHeader.getSequenceDictionary().size();
        indexBuilder = new BAMIndexBuilder(fileHeader);
        outputWriter = new BinaryBAMIndexWriter(numReferences, output);
    }

    /**
     * Prepare to index a BAM.
     * @param output Index will be written here.  output will be closed when finish() method is called.
     * @param fileHeader header for the corresponding bam file.
     */
    public BAMIndexer(final OutputStream output, SAMFileHeader fileHeader) {

        numReferences = fileHeader.getSequenceDictionary().size();
        indexBuilder = new BAMIndexBuilder(fileHeader);
        outputWriter = new BinaryBAMIndexWriter(numReferences, output);
    }

    /**
     * Record any index information for a given BAM record.
     * If this alignment starts a new reference, write out the old reference.
     * Requires a non-null value for rec.getFileSource().
     *
     * @param rec The BAM record
     */
    public void processAlignment(final SAMRecord rec) {
        try {
            final int reference = rec.getReferenceIndex();
            if (reference != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX && reference != currentReference) {
                // process any completed references
                advanceToReference(reference);
            }
            indexBuilder.processAlignment(rec);
        } catch (Exception e) {
            throw new SAMException("Exception creating BAM index for record " + rec, e);
        }
    }

    /**
     * After all the alignment records have been processed, finish is called.
     * Writes any final information and closes the output file.
     */
    public void finish() {
        // process any remaining references
        advanceToReference(numReferences);
        outputWriter.writeNoCoordinateRecordCount(indexBuilder.getNoCoordinateRecordCount());
        outputWriter.close();
    }

    /** write out any references between the currentReference and the nextReference */
    private void advanceToReference(int nextReference) {
        while (currentReference < nextReference) {
            BAMIndexContent content = indexBuilder.processReference(currentReference);
            outputWriter.writeReference(content);
            currentReference++;
            indexBuilder.startNewReference();
        }
    }

    /**
     * Generates a BAM index file, either textual or binary, from an input BAI file.
     * Only used for testing, but located here for visibility into CachingBAMFileIndex.
     *
     * @param output     BAM Index (.bai) file (or bai.txt file when text)
     * @param textOutput Whether to create text output or binary
     */
    static public void createAndWriteIndex(final File input, final File output, final boolean textOutput) {

        // content is from an existing bai file.

        final CachingBAMFileIndex existingIndex = new CachingBAMFileIndex(input, null);
        final int n_ref = existingIndex.getNumberOfReferences();
        final BAMIndexWriter outputWriter;
        if (textOutput) {
            outputWriter = new TextualBAMIndexWriter(n_ref, output);
        } else {
            outputWriter = new BinaryBAMIndexWriter(n_ref, output);
        }

        // write the content one reference at a time
        try {
            for (int i = 0; i < n_ref; i++) {
                outputWriter.writeReference(existingIndex.getQueryResults(i));
            }
            outputWriter.writeNoCoordinateRecordCount(existingIndex.getNoCoordinateCount());
            outputWriter.close();

        } catch (Exception e) {
            throw new SAMException("Exception creating BAM index", e);
        }
    }

    /**
     * Class for constructing BAM index files.
     * One instance is used to construct an entire index.
     * processAlignment is called for each alignment until a new reference is encountered, then
     * processReference is called when all records for the reference have been processed.
     */
    private class BAMIndexBuilder {

        private final SAMFileHeader bamHeader;

        // the bins for the current reference
        private Bin[] bins; // made only as big as needed for each reference
        private int binsSeen = 0;

        // linear index for the current reference
        private final long[] index = new long[LinearIndex.MAX_LINEAR_INDEX_SIZE];
        private int largestIndexSeen = -1;

        // information in meta data
        private BAMIndexMetaData indexStats = new BAMIndexMetaData();

        /**
         * @param header SAMFileheader used for reference name (in index stats) and for max bin number
         */
        BAMIndexBuilder(SAMFileHeader header) {
            this.bamHeader = header;
        }

        /**
         * Record any index information for a given BAM record
         *
         * @param rec The BAM record. Requires rec.getFileSource() is non-null.
         */
        public void processAlignment(final SAMRecord rec) {

            // metadata
            indexStats.recordMetaData(rec);

            final int alignmentStart = rec.getAlignmentStart();
            if (alignmentStart == SAMRecord.NO_ALIGNMENT_START) {
                return; // do nothing for records without coordinates, but count them
            }

            // various checks
            final int reference = rec.getReferenceIndex();
            if (reference != currentReference) {
                throw new SAMException("Unexpected reference " + reference +
                        " when constructing index for " + currentReference + " for record " + rec);
            }

            // process bins

            final Integer binNumber = rec.getIndexingBin();
            final int binNum = binNumber == null ? rec.computeIndexingBin() : binNumber;

            // has the bins array been allocated? If not, do so
            if (bins == null) {
                final SAMSequenceRecord seq = bamHeader.getSequence(reference);
                if (seq == null) {
                    bins = new Bin[MAX_BINS + 1];
                } else {
                    bins = new Bin[AbstractBAMFileIndex.getMaxBinNumberForSequenceLength(seq.getSequenceLength()) + 1];
                }
            }

            // is there a bin already represented for this index?  if not, add one
            final Bin bin;
            if (bins[binNum] != null) {
                bin = bins[binNum];
            } else {
                bin = new Bin(reference, binNum);
                bins[binNum] = bin;
                binsSeen++;
            }

            // process chunks

            final SAMFileSource source = rec.getFileSource();
            if (source == null) {
                throw new SAMException("No source (virtual file offsets); needed for indexing on BAM Record " + rec);
            }
            final Chunk newChunk = ((BAMFileSpan) source.getFilePointer()).getSingleChunk();
            final long chunkStart = newChunk.getChunkStart();
            final long chunkEnd = newChunk.getChunkEnd();

            final List<Chunk> oldChunks = bin.getChunkList();
            if (!bin.containsChunks()) {
                bin.addInitialChunk(newChunk);

            } else {
                final Chunk lastChunk = bin.getLastChunk();

                // Coalesce chunks that are in the same or adjacent file blocks.
                // Similar to AbstractBAMFileIndex.optimizeChunkList,
                // but no need to copy the list, no minimumOffset, and maintain bin.lastChunk
                if (BlockCompressedFilePointerUtil.areInSameOrAdjacentBlocks(lastChunk.getChunkEnd(),chunkStart)) {
                    lastChunk.setChunkEnd(chunkEnd);  // coalesced
                } else {
                    oldChunks.add(newChunk);
                    bin.setLastChunk(newChunk);
                }
            }

            // process linear index

            // the smallest file offset that appears in the 16k window for this bin
            final int alignmentEnd = rec.getAlignmentEnd();
            int startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart); // the 16k window
            final int endWindow;

            if (alignmentEnd == SAMRecord.NO_ALIGNMENT_START) {   // assume alignment uses one position
                // Next line for C (samtools index) compatibility. Differs only when on a window boundary
                startWindow = LinearIndex.convertToLinearIndexOffset(alignmentStart - 1);
                endWindow = startWindow;
            } else {
                endWindow = LinearIndex.convertToLinearIndexOffset(alignmentEnd);
            }

            if (endWindow > largestIndexSeen) {
                largestIndexSeen = endWindow;
            }

            // set linear index at every 16K window that this alignment overlaps
            for (int win = startWindow; win <= endWindow; win++) {
                if (index[win] == 0 || chunkStart < index[win]) {
                    index[win] = chunkStart;
                }
            }
        }

        /**
         * Creates the BAMIndexContent for this reference.
         * Requires all alignments of the reference have already been processed.
         */
        public BAMIndexContent processReference(int reference) {

            if (reference != currentReference) {
                throw new SAMException("Unexpected reference " + reference + " when constructing index for " + currentReference);
            }

            // process bins
            if (binsSeen == 0) return null;  // no bins for this reference

            // process chunks
            // nothing needed

            // process linear index
            // linear index will only be as long as the largest index seen
            final long[] newIndex = new long[largestIndexSeen + 1]; // in java1.6 Arrays.copyOf(index, largestIndexSeen + 1);

            // C (samtools index) also fills in intermediate 0's with values.  This seems unnecessary, but safe
            long lastNonZeroOffset = 0;
            for (int i = 0; i <= largestIndexSeen; i++) {
                if (index[i] == 0) {
                    index[i] = lastNonZeroOffset; // not necessary, but C (samtools index) does this
                    // note, if you remove the above line BAMIndexWriterTest.compareTextual and compareBinary will have to change
                } else {
                    lastNonZeroOffset = index[i];
                }
                newIndex[i] = index[i];
            }

            final LinearIndex linearIndex = new LinearIndex(reference, 0, newIndex);

            return new BAMIndexContent(reference, bins, binsSeen, indexStats, linearIndex);
        }

        /**
         * @return the count of records with no coordinate positions
         */
        public long getNoCoordinateRecordCount() {
            return indexStats.getNoCoordinateRecordCount();
        }

        /**
         * reinitialize all data structures when the reference changes
         */
        void startNewReference() {
            bins = null;
            if (binsSeen > 0) {
                Arrays.fill(index, 0);
            }
            binsSeen = 0;
            largestIndexSeen = -1;
            indexStats.newReference();
        }
    }
}