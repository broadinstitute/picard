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

import java.io.File;
import java.io.OutputStream;

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
    public BAMIndexer(final File output, final SAMFileHeader fileHeader) {

        numReferences = fileHeader.getSequenceDictionary().size();
        indexBuilder = new BAMIndexBuilder(fileHeader.getSequenceDictionary());
        outputWriter = new BinaryBAMIndexWriter(numReferences, output);
    }

    /**
     * Prepare to index a BAM.
     * @param output Index will be written here.  output will be closed when finish() method is called.
     * @param fileHeader header for the corresponding bam file.
     */
    public BAMIndexer(final OutputStream output, final SAMFileHeader fileHeader) {

        numReferences = fileHeader.getSequenceDictionary().size();
        indexBuilder = new BAMIndexBuilder(fileHeader.getSequenceDictionary());
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
        } catch (final Exception e) {
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
    private void advanceToReference(final int nextReference) {
        while (currentReference < nextReference) {
            final BAMIndexContent content = indexBuilder.processReference(currentReference);
            outputWriter.writeReference(content);
            currentReference++;
            if (currentReference < numReferences) {
                indexBuilder.startNewReference();
            }
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

        } catch (final Exception e) {
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

        private final SAMSequenceDictionary sequenceDictionary;

        private BinningIndexBuilder binningIndexBuilder;

        private int currentReference = -1;

        // information in meta data
        private final BAMIndexMetaData indexStats = new BAMIndexMetaData();

        BAMIndexBuilder(final SAMSequenceDictionary sequenceDictionary) {
            this.sequenceDictionary = sequenceDictionary;
            if (!sequenceDictionary.isEmpty()) startNewReference();
        }

        /**
         * Record any index information for a given BAM record
         *
         * @param rec The BAM record. Requires rec.getFileSource() is non-null.
         */
        public void processAlignment(final SAMRecord rec) {

            // metadata
            indexStats.recordMetaData(rec);

            if (rec.getAlignmentStart() == SAMRecord.NO_ALIGNMENT_START) {
                return; // do nothing for records without coordinates, but count them
            }

            // various checks
            final int reference = rec.getReferenceIndex();
            if (reference != currentReference) {
                throw new SAMException("Unexpected reference " + reference +
                        " when constructing index for " + currentReference + " for record " + rec);
            }

            binningIndexBuilder.processFeature(new BinningIndexBuilder.FeatureToBeIndexed() {
                @Override
                public int getStart() {
                    return rec.getAlignmentStart();
                }

                @Override
                public int getEnd() {
                    return rec.getAlignmentEnd();
                }

                @Override
                public Integer getIndexingBin() {
                    final Integer binNumber = rec.getIndexingBin();
                    return (binNumber == null ? rec.computeIndexingBin() : binNumber);

                }

                @Override
                public Chunk getChunk() {
                    final SAMFileSource source = rec.getFileSource();
                    if (source == null) {
                        throw new SAMException("No source (virtual file offsets); needed for indexing on BAM Record " + rec);
                    }
                    return ((BAMFileSpan) source.getFilePointer()).getSingleChunk();
                }
            });

        }

        /**
         * Creates the BAMIndexContent for this reference.
         * Requires all alignments of the reference have already been processed.
         * @return Null if there are no features for this reference.
         */
        public BAMIndexContent processReference(final int reference) {

            if (reference != currentReference) {
                throw new SAMException("Unexpected reference " + reference + " when constructing index for " + currentReference);
            }

            final BinningIndexContent indexContent = binningIndexBuilder.generateIndexContent();
            if (indexContent == null) return null;
            return new BAMIndexContent(indexContent.getReferenceSequence(), indexContent.getBins(),
                    indexStats, indexContent.getLinearIndex());

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
            ++currentReference;
            // I'm not crazy about recycling this object, but that is the way it was originally written and
            // it helps keep track of no-coordinate read count (which shouldn't be stored in this class anyway).
            indexStats.newReference();
            binningIndexBuilder = new BinningIndexBuilder(currentReference,
                    sequenceDictionary.getSequence(currentReference).getSequenceLength());
        }
    }
}