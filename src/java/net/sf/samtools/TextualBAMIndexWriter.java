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
import java.util.Arrays;
import java.util.List;

/**
 * Class for writing binary BAM index files as human-readable text.
 * Used for testing only.

 */
class TextualBAMIndexWriter extends AbstractBAMIndexWriter {

    private final PrintWriter pw;
    private final boolean sortBins;

    /**
     * constructor
     *
     * @param n_ref    Number of reference sequences
     * @param output   BAM Index output file
     * @param sortBins Whether to sort the bins - useful for comparison to c-generated index
     */
    public TextualBAMIndexWriter(final int n_ref, final File output, final boolean sortBins) {
        super(output, n_ref);
        this.sortBins = sortBins;
        try {
            pw = new PrintWriter(output);
        } catch (FileNotFoundException e) {
            throw new SAMException("Can't find output file " + output, e);
        }
    }

    /**
     * Write header information at the beginning of the file
     */
    public void writeHeader() {
        pw.println("n_ref=" + n_ref);
    }

    /**
     * Write this content as human-readable text
     */
    public void writeReference(final BAMIndexContent content, int reference) {

        if (content == null) {
            writeNullContent(pw, reference);
            return;
        }
        final int ref = content.getReferenceSequence();
        if (ref != reference){
            // something is wrong. either ignore it, or add nulls till we catch up
            System.err.println("Reference on content is " + ref + " but expecting reference " + reference);
            if (ref > reference){
                for (int i = reference; i< ref; i++){
                     writeNullContent(pw, i);
                }
            }
            // if (ref < reference) // just ignore it?? todo
        }
        final List<Bin> bins = content.getBins();
        final LinearIndex linearIndex = content.getLinearIndex();

       // bins
        if (bins == null || bins.size() == 0) {
            writeNullContent(pw, ref);
            return;
        }

        final int size = bins.size();
        pw.println("Reference " + ref + " has n_bin= " + size);

        // copy into an array so that it can be sorted
        final Bin[] binArray = new Bin[size];
        if (size != 0) {
            bins.toArray(binArray);
        }
        if (sortBins) Arrays.sort(binArray);  // Sort for easy text comparisons

        // chunks
        for (int j = 0; j < size; j++) {
            if (binArray[j].getChunkList() == null) {
                 pw.println("  Ref " + reference + " bin " + binArray[j].getBinNumber() + " has no binArray");  // remove?
                continue;
            }
            final List<Chunk> chunkList = binArray[j].getChunkList();
            if (chunkList == null) {
                pw.println("  Ref " + ref + " bin " + binArray[j].getBinNumber() + " has no chunkList");
                continue;
            }
            if (chunkList.size() == 0) {
                // continue;         // todo?
                // pw.println();
            }
            pw.print("  Ref " + ref + " bin " + binArray[j].getBinNumber() + " has n_chunk= " + chunkList.size());
            if (chunkList.size() == 0) {
                // continue;         // todo?
                 pw.println();
            }
            for (final Chunk c : chunkList) {
                pw.println("     Chunk: " + c.toString() +
                        " start: " + Long.toString(c.getChunkStart(), 16) +
                        " end: " + Long.toString(c.getChunkEnd(), 16));
            }
        }
        // linear index
        if (linearIndex == null || linearIndex.getIndexEntries() == null) {
            pw.println("Reference " + ref + " has n_intv= 0");
            return;
        }
        final long[] entries = linearIndex.getIndexEntries();
        final int indexStart = linearIndex.getIndexStart();
        // System.out.println("index start is " + indexStart);
        final int n_intv = entries.length + indexStart;
        pw.println("Reference " + ref + " has n_intv= " + n_intv);
        for (int k = 0; k < entries.length; k++) {
            if (entries[k] != 0) {
                pw.println("  Ref " + ref + " ioffset for " + (k + indexStart) + " is " + Long.toString(entries[k]));
            }
        }
        pw.flush ();  // write each reference to disk as it's being created
    }

    private void writeNullContent(PrintWriter pw, int reference) {
        pw.println("Reference " + reference + " has n_bin=0");
        pw.println("Reference " + reference + " has n_intv=0");
    }

    /**
     * Any necessary processing at the end of the file
     *
     * @param noCoordinateCount the count of records seen with no coordinate positions in the start coordinate
     */
    public void close(final Long noCoordinateCount) {
        pw.println("No Coordinate Count=" + noCoordinateCount);
        pw.close();
    }
}