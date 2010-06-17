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
import java.util.*;

/**
 * Class for constructing and writing BAM index files
 */
public class BAMFileIndexWriter extends AbstractBAMFileIndex{
// note, alternatively this could extend DiskBasedBAMFileIndex or CachingBAMFileIndex both of which extend AbstractBAMFileIndex
// or it could just implement BAMIndex (with a few additions noted below)

    // One element for each reference sequence.  Not populated until after call to finish()
    private final BAMIndexContent[] content;

    /**
     * The number of references (chromosomes) in the BAM file
     * The name n_ref corresponds to the name in the SAM Format Specification Document
     */
    private final int n_ref;

    /**
     * A mapping of reference sequence index to list of bins.
     */
    private final SortedMap<Integer, List<Bin>> referenceToBins = new TreeMap<Integer, List<Bin>>();

    /**
     * A mapping from bin to the chunks contained in that bin.
     */
    private final SortedMap<Bin, List<Chunk>> binToChunks = new TreeMap<Bin, List<Chunk>>();

    /**
     * A mapping of reference sequence index to linear index in progress
     */
    private final SortedMap<Integer, List<LinearIndexItem>> referenceToIndexEntries =
            new TreeMap<Integer, List<LinearIndexItem>>();

    private final File OUTPUT;
    
    /**
     * Constructor
     * @param OUTPUT      BAM Index (.bai) file (or bai.txt file when text)
     * @param nReferences Number of references in the input BAM file
     */
    public BAMFileIndexWriter(final File OUTPUT, final int nReferences) {
        // super(OUTPUT); // don't do this; it will try to open the non-existent output file
        super(null);
        this.OUTPUT = OUTPUT;
        this.n_ref = nReferences;
        this.content = new BAMIndexContent[n_ref];
    }

    /**
     * Generates a BAM index file, either textual or binary, sorted or not, from an input BAM file
     *
     * @param INPUT      BAM file
     * @param createText Whether to create text output or binary
     * @param sortBins   Whether to sort the bins in the output
     * @return Count of aligned records processed
     */
    public int createIndex(final File INPUT, final boolean createText, final boolean sortBins) {

        final SAMFileReader reader = new SAMFileReader(INPUT);
        reader.enableFileSource(true);

        int count = 0;
        int localCount = 0;
        int totalRecords = 0;

        int reference = 0;
        int lastReference = 0;

        CloseableIterator<SAMRecord> it = reader.iterator();
        while (it.hasNext()) {
            totalRecords++;
            if (totalRecords % 1000000 == 0) {
                verbose(totalRecords + " reads processed ...");
            }
            final SAMRecord rec = it.next();
            if (rec.getAlignmentStart() == 0)
                continue;  // do nothing for un-aligned records

            reference = rec.getReferenceIndex();
            if (reference != lastReference) {
                // switching chromosome/reference
                // verbose("     " + localCount + " aligned records in reference " + lastReference);
                processReference(lastReference); // or call finish at the end
                lastReference = reference;
                localCount = 0;
            }
            count++;
            localCount++;
            processAlignment(rec);
        }
        it.close();
        processReference(lastReference);
        // verbose("     " + localCount + " aligned records in reference " + lastReference);
        // finish();  // alternative to processReference calls above
        verbose("There are " + count + " aligned records in the input BAM file " + INPUT);

        deleteIndex(); // delete old index if present. Caller has checked that this is ok and expected
        try {
            if (createText) {
                writeText(sortBins);
            } else {
                writeBinary(sortBins, INPUT.length());
            }
        } catch (Exception e) {
            deleteIndex();
            throw new SAMException("Exception creating index " + e.getMessage());
        }
        return count;
    }

    /**
     * Deletes old or partial index file
     * Called whenever exceptions occur.
     */
    public void deleteIndex(){
            OUTPUT.delete();
    }

    public void writeText(boolean sortBins) throws Exception {
        writeText(n_ref, OUTPUT, sortBins);
    }

    public void writeBinary(final boolean sortBins, final long bamFileSize) throws Exception {
        writeBinary(n_ref, OUTPUT, sortBins, bamFileSize);
    }


    /**
     * Record any index information for a given BAM record
     *
     * @param rec The BAM record
     */
    public void processAlignment(final SAMRecord rec) {

        final int binNumber;
        if (rec.getIndexingBin() == null) {
            binNumber = rec.computeIndexingBin();
        } else {
            binNumber = rec.getIndexingBin();
        }

        // is there a bin already represented for this index?  if not, add it
        final int reference = rec.getReferenceIndex();
        final Bin bin = new Bin(reference, binNumber);
        List<Bin> binList = referenceToBins.get(reference);
        if (binList == null) {
            binList = new ArrayList<Bin>();
            referenceToBins.put(reference, binList);
        }
        if (!binList.contains(bin)) {
            binList.add(bin);
        }

        // add chunk information from this record to the bin
        final BAMFileSpan newSpan;
        if (rec.getFileSource() == null) {
            // todo throw new SAMException? deleteIndex()?
            System.err.println("No source for BAM Record " + rec);
            return;
        } else {
            newSpan = (BAMFileSpan) rec.getFileSource().getFilePointer();
        }

        if (newSpan != null) {
            List<Chunk> chunks = binToChunks.get(bin);
            if (chunks == null) {
                chunks = new ArrayList<Chunk>();
                binToChunks.put(bin, chunks);
            }
            for (Chunk c : newSpan.getChunks()) {
                chunks.add(c);
            }
        } else {
            deleteIndex();
            throw new SAMException("Null coordinates on record " + rec +
                    " in reference " + reference);
        }

        // add linear index information
        // the smallest file offset that starts in the 16k window for this bin
        // if (binNumber < LEVEL_STARTS[LEVEL_STARTS.length - 1]) {  // don't index the top level 6 16k bins with bin# >= 4681
        if (binNumber < getFirstBinInLevel(getNumIndexLevels()-1)) {  // don't index the top level 6 16k bins with bin# >= 4681
            final long iOffset = ((BAMFileSpan) rec.getFileSource().getFilePointer()).toCoordinateArray()[0];
            final int alignmentStart = rec.getAlignmentStart() - 1;
            final int window = LinearIndex.convertToLinearIndexOffset(alignmentStart); // the 16k window

            /* log.debug("Reference " + reference + " bin=" + binNumber + " startPosition=" +
                    iOffset + "(" + Long.toString(iOffset, 16) + "x)" + " startAlignment=" +
                    alignmentStart + "(" + Long.toString(alignmentStart,16) + "x)" +
                    " window " + window);
            */

            List<LinearIndexItem> indexList = referenceToIndexEntries.get(reference);
            if (indexList == null) {
                referenceToIndexEntries.put(reference, new ArrayList<LinearIndexItem>());
                indexList = referenceToIndexEntries.get(reference);
            }
            // Has this window already been seen? If so, might replace existing index
            boolean replaceIndex = true;
            for (final LinearIndexItem li : indexList) {
                if (li.window == window) {
                    if (li.offset < iOffset) {
                        replaceIndex = false;  // existing window value is already minimum
                    } else {
                        // remove existingIndex in place of this one
                        indexList.remove(li);
                    }
                    break;
                }
            }
            if (replaceIndex)
                indexList.add(new LinearIndexItem(window, iOffset));
        }
    }

    /**
     * After all the alignment records have been processed, finish is called.
     * Note, we can do this processing per-reference instead of per-file if desired
     */
    public void finish() {
        // constructs a BAMIndexContent from bins, chunks, and linear index
        for (int i = 0; i < n_ref; i++) {
            processReference(i);
        }
    }

    /**  Processing after all alignments of a reference have already been processed */
    private void processReference(int reference) {

        // process bins
        final List<Bin> binList = referenceToBins.get(reference);

        if (binList == null) return;  // no bins for this reference

        // process chunks
        for (Bin bin : binList) {
            List<Chunk> chunkList = binToChunks.get(bin);
            if (chunkList != null) {
                chunkList = optimizeChunkList(chunkList, 0);
                binToChunks.put(bin, chunkList);
            }
        }

        // process linear index
        final LinearIndex linearIndex = computeLinearIndex(reference);

        content[reference] = new BAMIndexContent(reference, binList, binToChunks, linearIndex);
    }

    /** transform the sparse representation of linear index to an array representation  */
    private LinearIndex computeLinearIndex(int reference) {
        final LinearIndex linearIndex;
        final List<LinearIndexItem> indexList = referenceToIndexEntries.get(reference);
        if (indexList == null) {
            // skip linear index for this reference
            linearIndex = null;
        } else {
            // get the max window in the list; linear index will be this long
            int maxWindow = 0;
            for (final LinearIndexItem l : indexList) {
                maxWindow = Math.max(maxWindow, l.window);
            }
            // log.debug("** n_intv (maxWindow) for reference " + reference + " is " + (maxWindow + 1));
            long[] newIndex = new long[maxWindow + 1]; // newIndex is one-based
            Arrays.fill(newIndex, 0);

            // non-0 linearIndex for each entry
            for (final LinearIndexItem l : indexList) {
                newIndex[l.window] = l.offset;
            }
            linearIndex = new LinearIndex(reference, 1, newIndex); // newIndex is one-based
        }
        return linearIndex;
    }

    /** Return the content of this index for a given reference sequence */
    protected BAMIndexContent getQueryResults(int reference){
        return content[reference];
    }

    /**
     * This method is needed if this class extends AbstractBAMFileIndex;
     * not needed if it extends DiskBasedBAMFileIndex or CachingBAMFileIndex
     */
    public BAMFileSpan getSpanOverlapping(final int referenceIndex, final int startPos, final int endPos) {
        throw new UnsupportedOperationException();
    }

    /*  These needed if this class just implements BAMIndex
    public void open() {}
    public void close() {} // call finish() here?
    public long getStartOfLastLinearBin() {throw new UnsupportedOperationException();};
    private static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};
    */

    // net.sf.samtools logging is primitive ...
    private void verbose(String message) {
        boolean verbose = true;
        if (verbose) {
            System.out.println ("BAMFileIndexWriter: " + message);
        }
    }

    /**
     * This class facilitates a sparse representation of the linear index during construction,
     * when the size of the linear index is unknown.
     */
    private class LinearIndexItem {

        public final int window;   // index position in the array
        public final long offset;  // value of the linear index at that position

        public LinearIndexItem(int window, long offset) {
            this.window = window;
            this.offset = offset;
        }
    }
}
