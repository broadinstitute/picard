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
package net.sf.samtools;

import net.sf.samtools.util.BlockCompressedFilePointerUtil;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Metadata about the bam index contained within the bam index.
 * One instance created per index file.
 */
public class BAMIndexMetaData {

    // information for the entire index.
    // stored at the end of the index
    private long noCoordinateRecords = 0;

    // information for each reference.
    // stored in two chunks in bin # MAX_BINS
    private long firstOffset = -1;
    private long lastOffset = 0;
    private int alignedRecords = 0;
    private int unAlignedRecords = 0;  // unmapped, but associated with this reference


    /**
     * Constructor used when writing an index
     * construct one instance for each index generated
     */
    BAMIndexMetaData() {
        noCoordinateRecords = 0;
        newReference();
    }

    /**
     * Constructor used when reading an index
     * construct one instance for each index generated
     */
    BAMIndexMetaData(List<Chunk> chunkList) {
        noCoordinateRecords = 0;

        if (chunkList == null || chunkList.size() == 0) {
            // System.out.println("No metadata chunks");
        } else if (chunkList.size() != 2) {
            throw new SAMException("Unexpected number of metadata chunks " + (chunkList.size()));
        }
        // fill in the first/lastOffset un/alignedRecords from this
        boolean firstChunk = true;
        if (chunkList != null) {
            for (Chunk c : chunkList) {
                long start = c.getChunkStart();
                long end = c.getChunkEnd();
                if (firstChunk) {
                    firstOffset = start;
                    lastOffset = end;
                    firstChunk = false;
                } else {
                    firstChunk = true;
                    alignedRecords = (int) start;
                    unAlignedRecords = (int) end;
                }
            }
        }
    }

    /**
     * @return the count of aligned records associated with this reference
     */
    public int getAlignedRecordCount() {
        return alignedRecords;
    }

    /**
     * @return the count of unaligned records associated with this reference
     */
    public int getUnalignedRecordCount() {
        return unAlignedRecords;
    }

    /**
     * Call for each new reference sequence encountered
     */
    void newReference() {
        firstOffset = -1;
        lastOffset = 0;
        alignedRecords = 0;
        unAlignedRecords = 0;
    }

    /**
     * Extract relevant metaData from the record and its filePointer
     * Call only once per record in the file being indexed
     *
     * @param rec
     */
    void recordMetaData(final SAMRecord rec) {

        final int alignmentStart = rec.getAlignmentStart();
        if (alignmentStart == SAMRecord.NO_ALIGNMENT_START) {
            incrementNoCoordinateRecordCount();
            return;
        }

        if (rec.getFileSource() == null){
            throw new SAMException("BAM cannot be indexed without setting a fileSource for record " + rec);
        }
        final Chunk newChunk = ((BAMFileSpan) rec.getFileSource().getFilePointer()).getSingleChunk();
        final long start = newChunk.getChunkStart();
        final long end = newChunk.getChunkEnd();

        if (rec.getReadUnmappedFlag()) {
            unAlignedRecords++;
        } else {
            alignedRecords++;
        }
        if (BlockCompressedFilePointerUtil.compare(start, firstOffset) < 1 || firstOffset == -1) {
            this.firstOffset = start;
        }
        if (BlockCompressedFilePointerUtil.compare(lastOffset, end) < 1) {
            this.lastOffset = end;
        }
    }

    /**
     * Call whenever a reference with no coordinate information is encountered in the bam file
     */
    void incrementNoCoordinateRecordCount() {
        noCoordinateRecords++;
    }

    /**
     * Set local variable. Normally noCoordinateRecord count accessed from AbstractBAMFileIndex when reading
     */
    private void setNoCoordinateRecordCount(long count) {
        noCoordinateRecords = count;
    }


    /**
     * @return the count of records with no coordinate information in the bam file.
     * Not public, since only used by BAMIndexer when writing bam index.
     * Readers of bam index should use AbstractBAMFileIndex.getNoCoordinateRecordCount.
     */
    long getNoCoordinateRecordCount() {
        return noCoordinateRecords;
    }

    /**
     * @return the first virtual file offset used by this reference
     */
    long getFirstOffset() {
        return firstOffset;
    }

    /**
     * @return the last virtual file offset used by this reference
     */
    long getLastOffset() {
        return lastOffset;
    }

    /**
     * Prints meta-data statistics from BAM index (.bai) file
     * Statistics include count of aligned and unaligned reads for each reference sequence
     * and a count of all records with no start coordinate
     */
    static public void printIndexStats(final File inputBamFile) {
        try {
            final BAMFileReader bam = new BAMFileReader(inputBamFile, null, false, SAMFileReader.ValidationStringency.SILENT, new DefaultSAMRecordFactory());
            if (!bam.hasIndex()) {
                throw new SAMException("No index for bam file " + inputBamFile);
            }
            BAMIndexMetaData[] data = getIndexStats(bam);
            // read through all the bins of every reference.
            int nRefs = bam.getFileHeader().getSequenceDictionary().size();
            for (int i = 0; i < nRefs; i++) {
                final SAMSequenceRecord seq = bam.getFileHeader().getSequence(i);
                if (seq == null) continue;
                final String sequenceName = seq.getSequenceName();
                final int sequenceLength = seq.getSequenceLength();
                System.out.print(sequenceName + ' ' + "length=\t" + sequenceLength);
                if (data[i] == null) {
                    System.out.println();
                    continue;
                }
                System.out.println("\tAligned= " + data[i].getAlignedRecordCount() +
                        "\tUnaligned= " + data[i].getUnalignedRecordCount());
            }
            System.out.println("NoCoordinateCount= " + data[0].getNoCoordinateRecordCount());
        } catch (IOException e) {
            throw new SAMException("Exception in getting index statistics", e);
        }
    }

    /**
     * Prints meta-data statistics from BAM index (.bai) file
     * Statistics include count of aligned and unaligned reads for each reference sequence
     * and a count of all records with no start coordinate
     */
    static public BAMIndexMetaData[] getIndexStats(final BAMFileReader bam){

        AbstractBAMFileIndex index = (AbstractBAMFileIndex) bam.getIndex();
        // read through all the bins of every reference.
        int nRefs = index.getNumberOfReferences();
        BAMIndexMetaData[] result = new BAMIndexMetaData[nRefs == 0 ? 1 : nRefs];
        for (int i = 0; i < nRefs; i++) {
            result[i] = index.getMetaData(i);
        }

        if (result[0] == null){
           result[0] = new BAMIndexMetaData();
        }
        final Long noCoordCount = index.getNoCoordinateCount();
        if (noCoordCount != null)  // null in old index files without metadata
           result[0].setNoCoordinateRecordCount(noCoordCount);

        return result;
    }
}
