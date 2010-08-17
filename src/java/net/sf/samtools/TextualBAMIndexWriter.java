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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;

/**
 * Class for writing binary BAM index files as human-readable text.
 * Used for testing only.

 */
class TextualBAMIndexWriter implements BAMIndexWriter {

    protected final int nRef;
    protected final File output;
    private final PrintWriter pw;
    private int count = 0;

    /**
     * constructor
     *
     * @param nRef    Number of reference sequences
     * @param output   BAM Index output file
     */
    public TextualBAMIndexWriter(final int nRef, final File output) {
        this.output = output;
        this.nRef = nRef;
        try {
            pw = new PrintWriter(output);
        } catch (FileNotFoundException e) {
            throw new SAMException("Can't find output file " + output, e);
        }
        writeHeader();
    }

    /**
     * Write header information at the beginning of the file
     */
    private void writeHeader() {
        pw.println("n_ref=" + nRef);
    }

    /**
     * Write this content as human-readable text
     */
    public void writeReference(final BAMIndexContent content) {

        final int reference = content.getReferenceSequence();

        if (content == null) {
            writeNullContent(reference);
            count++;
            return;
        }

        if (reference != count){
            throw new SAMException("Reference on content is " + reference + " but expecting reference " + count);
        }
        count++;

        final BAMIndexContent.BinList bins = content.getBins();
        final int size = bins == null ? 0 : content.getNumberOfNonNullBins();

        if (size == 0) {
            writeNullContent(reference);
            return;
        }

        //final List<Chunk> chunks = content.getMetaData() == null ? null
        //        : content.getMetaData().getMetaDataChunks();
        BAMIndexMetaData metaData = content.getMetaData();

        pw.println("Reference " + reference + " has n_bin= " + Integer.toString(size + (metaData != null? 1 : 0)));

        // chunks
        for (Bin bin : bins) {   // note, bins will always be sorted
            if (bin.getBinNumber() == AbstractBAMFileIndex.MAX_BINS)  break;
            if (bin.getChunkList() == null) {
                pw.println("  Ref " + reference + " bin " + bin.getBinNumber() + " has no binArray");  // remove?
                continue;
            }
            final List<Chunk> chunkList = bin.getChunkList();
            if (chunkList == null) {
                pw.println("  Ref " + reference + " bin " + bin.getBinNumber() + " has no chunkList");
                continue;
            }
            pw.print("  Ref " + reference + " bin " + bin.getBinNumber() + " has n_chunk= " + chunkList.size());
            if (chunkList.size() == 0) {
                 pw.println();
            }
            for (final Chunk c : chunkList) {
                pw.println("     Chunk: " + c.toString() +
                        " start: " + Long.toString(c.getChunkStart(), 16) +
                        " end: " + Long.toString(c.getChunkEnd(), 16));
            }
        }

        writeChunkMetaData(reference, metaData);
        
        // linear index
        final LinearIndex linearIndex = content.getLinearIndex();
        if (linearIndex == null || linearIndex.getIndexEntries() == null) {
            pw.println("Reference " + reference + " has n_intv= 0");
            return;
        }
        final long[] entries = linearIndex.getIndexEntries();
        final int indexStart = linearIndex.getIndexStart();
        // System.out.println("index start is " + indexStart);
        final int n_intv = entries.length + indexStart;
        pw.println("Reference " + reference + " has n_intv= " + n_intv);
        for (int k = 0; k < entries.length; k++) {
            if (entries[k] != 0) {
                pw.println("  Ref " + reference + " ioffset for " + (k + indexStart) + " is " + Long.toString(entries[k]));
            }
        }
        pw.flush ();  // write each reference to disk as it's being created
    }

    /**
     * Write the meta data represented by the chunkLists associated with bin MAX_BINS 37450
     *
     * @param metaData information describing numAligned records, numUnAligned, etc
     */
    private void writeChunkMetaData(int reference, BAMIndexMetaData metaData) {
        final int nChunks = metaData == null ? 0 : 2;
        pw.print("  Ref " + reference + " bin 37450 has n_chunk= " + nChunks);
        if (nChunks == 0) {
            pw.println();
        } else {
            pw.println("     Chunk: " + //  c.toString() +
                    " start: " + Long.toString(metaData.getFirstOffset(), 16) +
                    " end: " + Long.toString(metaData.getLastOffset(), 16));
            pw.println("     Chunk: " + //  c.toString() +
                    " start: " + Long.toString(metaData.getAlignedRecordCount(), 16) +
                    " end: " + Long.toString(metaData.getUnalignedRecordCount(), 16));
        }

    }
       
    private void writeNullContent(int reference) {
        pw.println("Reference " + reference + " has n_bin=0");
        pw.println("Reference " + reference + " has n_intv=0");
    }

    /**
     * Write count of records without coordinates
     *
     * @param noCoordinateCount the count of records seen with no coordinate positions in the start coordinate
     */
    public void writeNoCoordinateRecordCount(final Long noCoordinateCount) {
        pw.println("No Coordinate Count=" + noCoordinateCount);
    }

    /**
     * Any necessary processing at the end of the file
     */
    public void close() {
        pw.close();
    }
}
