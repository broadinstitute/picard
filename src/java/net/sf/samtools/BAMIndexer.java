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

import net.sf.samtools.util.BlockCompressedInputStream;

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
     * Record any index information for a given BAM record.
     * If this alignment starts a new reference, write out the old reference
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
     * Prints meta-data statistics from BAM index (.bai) file
     * Statistics include count of aligned and unaligned reads for each reference sequence
     * and a count of all records with no start coordinate
     */
    static public void printIndexStats(final File inputBamFile) {
        try {
            final BAMFileReader bam = new BAMFileReader(inputBamFile, null, false, SAMFileReader.ValidationStringency.SILENT);

            if (!bam.hasIndex()) {
                throw new SAMException("No index for bam file " + inputBamFile);
            }

            AbstractBAMFileIndex index = (AbstractBAMFileIndex) bam.getIndex();
            index.open();
            // read through all the bins of every reference.
            int nRefs = index.getNumberOfReferences();
            for (int i = 0; i < nRefs; i++) {

                final SAMSequenceRecord seq = bam.getFileHeader().getSequence(i);
                if (seq == null) continue;
                final String sequenceName = seq.getSequenceName();
                final int sequenceLength = seq.getSequenceLength();
                System.out.print(sequenceName + ' ' + "length=\t" + sequenceLength);

                BAMIndexContent content = index.query(i, 0, -1); // todo: it would be faster just to skip to the last bin

                if (content == null || content.getBins() == null) {
                    System.out.println();
                    continue;
                }
                List<Chunk> chunkList = content.getMetaDataChunks();
                if (chunkList == null || chunkList.size() == 0) {
                    // System.out.println("No metadata chunks");
                } else if (chunkList.size() != 2) {
                    throw new SAMException("Unexpected number of metadata chunks " + (chunkList.size()));
                }
                boolean firstChunk = true;
                if (chunkList != null) {
                    for (Chunk c : chunkList) {
                        long start = c.getChunkStart();
                        long end = c.getChunkEnd();
                        if (firstChunk) {
                            // samtools idxstats doesn't print this, so we won't either
                            // System.out.print(sequenceName + ' ' + "Start=" + start + "    End=" + end);
                            firstChunk = false;
                        } else {
                            firstChunk = true;
                            System.out.println("\tAligned= " + start + "\tUnaligned= " + end);
                        }
                    }
                }
            }

            System.out.println("NoCoordinateCount= " + index.getNoCoordinateCount());

        } catch (IOException e) {
            throw new SAMException("Exception in getting index statistics", e);
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
            existingIndex.close();

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
        private BAMIndexStats indexStats = new BAMIndexStats();

        /**
         * @param header SAMFileheader used for reference name (in index stats) and for max bin number
         */
        BAMIndexBuilder(SAMFileHeader header) {
            this.bamHeader = header;
        }

        /**
         * Record any index information for a given BAM record
         *
         * @param rec The BAM record
         */
        public void processAlignment(final SAMRecord rec) {

            // various checks
            final int reference = rec.getReferenceIndex() == -1 ? currentReference : rec.getReferenceIndex();

            if (reference != currentReference) {
                throw new SAMException("Unexpected reference " + reference +
                        " when constructing index for " + currentReference + " for record " + rec);
            }

            // meta data processing
            final int alignmentStart = rec.getAlignmentStart();
            if (alignmentStart == SAMRecord.NO_ALIGNMENT_START) {
                indexStats.incrementNoCoordinateRecordCount();
                return; // do nothing for records without coordinates, but count them
            }

            if (rec.getReadUnmappedFlag()) {
                indexStats.incrementUnAlignedRecordCount();
            } else {
                indexStats.incrementAlignedRecordCount();
            }

            // process bins

            final Integer binNumber = rec.getIndexingBin();
            final int binNum = binNumber == null ? rec.computeIndexingBin() : binNumber;

            // has the bins array been allocated? If not, do so
            if (bins == null) {
                final SAMSequenceRecord seq = bamHeader == null ? null : bamHeader.getSequence(reference);
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
            final BAMFileSpan newSpan = (BAMFileSpan) source.getFilePointer();
            Chunk newChunk = newSpan.getSingleChunk();

            final long chunkStart = newChunk.getChunkStart();
            indexStats.setFirstOffsetIfSmaller(chunkStart);

            final long chunkEnd = newChunk.getChunkEnd();
            indexStats.setLastOffsetIfLarger(chunkEnd);

            List<Chunk> oldChunks = bin.getChunkList();
            if (oldChunks == null) {
                bin.addInitialChunk(newChunk);

            } else {
                final Chunk lastChunk = bin.getLastChunk();

                // Coalesce chunks that are in adjacent file blocks.
                // Similar to AbstractBAMFileIndex.optimizeChunkList,
                // but no need to copy the list, no minimumOffset, and maintain bin.lastChunk
                final long lastFileBlock = BlockCompressedInputStream.getFileBlock(lastChunk.getChunkEnd());
                final long chunkFileBlock = BlockCompressedInputStream.getFileBlock(chunkStart);
                if (chunkFileBlock - lastFileBlock <= 1) { // todo - possibility of overflow?
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
            final Bin metaBin = new Bin(reference, MAX_BINS);
            final List<Chunk> metaData = getMetaDataChunks();
            metaBin.setChunkList(metaData);

            // process linear index
            // linear index will only be as long as the largest index seen
            final long[] newIndex = new long[largestIndexSeen + 1]; // in java1.6 Arrays.copyOf(index, largestIndexSeen + 1);

            // C (samtools index) also fills in intermediate 0's with values.  This seems unnecessary, but safe
            long lastNonZeroOffset = 0;
            for (int i = 0; i <= largestIndexSeen; i++) {
                if (index[i] == 0) {
                    index[i] = lastNonZeroOffset; // not necessary, but C (samtools index) does this
                } else {
                    lastNonZeroOffset = index[i];
                }
                newIndex[i] = index[i];
            }

            final LinearIndex linearIndex = new LinearIndex(reference, 0, newIndex);

            return new BAMIndexContent(reference, bins, binsSeen + 1, metaData, linearIndex);
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

        /**
         * @return the per-reference list of meta data in the form of Chunks
         */
        private List<Chunk> getMetaDataChunks() {

            return indexStats.getMetaDataChunks();
        }

    }
}