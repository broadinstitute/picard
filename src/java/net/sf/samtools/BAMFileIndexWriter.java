/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * Class for writing BAM index files
 */
public class BAMFileIndexWriter extends DiskBasedBAMFileIndex {

    private final BAMIndexContent[] content;

    /**
     * The number of references (chromosomes) in the BA file
     * The name n_ref corresponds to the name in the SAM Format Specification Document
     */
    private final int n_ref;

    /**
     * A mapping of reference sequence index to set of bins.
     */
    private final SortedMap<Integer, Set<Bin>> binsInProgress = new TreeMap<Integer, Set<Bin>>();

     /**
     * A mapping from bin to the chunks contained in that bin.
     */
    // private final SortedMap<Bin, List<Chunk>> binToChunks = new TreeMap<Bin, List<Chunk>>(); // todo might be useful since chunk list is modifiable
    private final SortedMap<Bin, BAMFileSpan> binToSpan = new TreeMap<Bin, BAMFileSpan>();

    /**
     * A mapping of reference sequence index to linear index in progress
     */
    private final SortedMap<Integer, List<LinearIndexItem>> referenceToIndexEntries =
            new TreeMap<Integer, List<LinearIndexItem>>();

    //private final Log log = Log.getInstance(getClass());   // removed - so can be in samtools package

    private final File OUTPUT;

    /**
     * @param OUTPUT      BAM Index (.bai) file (or bai.txt file when text)
     * @param nReferences number of references in the input BAM file
     */
    public BAMFileIndexWriter(final File OUTPUT, final int nReferences) {
        //  super(OUTPUT); // don't do this; it will try to open the non-existant output file
        super(null);
        this.OUTPUT = OUTPUT;
        this.n_ref = nReferences;
        this.content = new BAMIndexContent[n_ref];
    }

    /**
     * Generates a BAM index file, either textual or binary, sorted or not from the input BAM file
     *
     * @param INPUT       BAM file
     * @param createText  Whether to create text output or binary
     * @param sortBins    Whether to sort the bins in the output
     * @return count of records processed
     * @throws Exception
     */
    public int createIndex(final File INPUT,
                           final boolean createText, final boolean sortBins) throws Exception {

        final SAMFileReader bam = new SAMFileReader(INPUT);
        bam.enableFileSource(true);

        int count = 0;
        int localCount = 0;
        int totalRecords = 0;

        int reference = 0;
        int lastReference = 0;

        for (Iterator<SAMRecord> iter = bam.iterator(); iter.hasNext();) {
            totalRecords++;
            if (totalRecords % 1000000 == 0) {
                verbose(totalRecords + " reads processed ...");
            }
            final SAMRecord rec = iter.next();
            if (rec.getAlignmentStart() == 0)
                continue;  // do nothing for un-aligned records

            reference = rec.getReferenceIndex();
            if (reference != lastReference) {
                // switching chromosome/reference
                lastReference = reference;
                // log.debug("     " + localCount + " aligned records in reference " + (reference - 1));
                // processReference(lastReference); // or call finish at the end
                localCount = 0;
            }
            count++;
            localCount++;
            processRecord(rec, localCount);
        }
        bam.close();
        // processReference(lastReference);
        finish();  // alternative to processReference calls above
        //log.debug("There are " + count + " aligned records in the input BAM file " + INPUT);

        if (OUTPUT.exists()) {
            OUTPUT.delete();
        }
        if (createText) {
            writeText(sortBins);
        } else {
            long size = INPUT.length();
            //log.info("The input file size is " + size);
            writeBinary(sortBins, size);
        }
        return count;
    }

    /**
     * Record any index information for a given BAM record
     *
     * @param rec        The SAM record
     * @param localCount The count within the reference.  Used for diagnostics only
     */
    private void processRecord(final SAMRecord rec, final int localCount) {

        final int bin = rec.getIndexingBin();   // todo: recompute using rec.computeIndexingBin?

        // is there a bin already represented for this index?  if not, add it
        final int reference = rec.getReferenceIndex();
        final Bin b = new Bin(reference, bin);
        if (binsInProgress.get(reference) == null) {
            binsInProgress.put(reference, new HashSet<Bin>());
        }
        binsInProgress.get(reference).add(b);

        // add chunk information from this record to the bin
        final BAMFileSpan newSpan = (BAMFileSpan) rec.getFileSource().getFilePointer();

        if (newSpan != null) {
            //List<Chunk> chunks = binToChunks.get(b);
            final BAMFileSpan oldSpan = binToSpan.get(b);
            final List<Chunk> oldChunks = oldSpan == null? null : oldSpan.getChunks();
            // copy the unmodifiable List from the span to a list that we can add to
            List<Chunk> chunks = new ArrayList();
            if (oldChunks != null) {
                for (Chunk c :oldChunks) {
                    chunks.add(c);
                }
            }
            for (Chunk c : newSpan.getChunks()) {
                chunks.add(c);
            }
            binToSpan.put(b, new BAMFileSpan(chunks));
        } else {
            throw new SAMException("Null coordinates on record " + rec +
                    " in reference " + reference + " count " + localCount);
        }

        // add linear index information
        // the smallest file offset that starts in the 16k window for this bin
        if (bin < LEVEL_STARTS[LEVEL_STARTS.length - 1]) {  // don't index the top level 6 16k bins with bin# >=4681
            final long ioffset = ((BAMFileSpan) rec.getFileSource().getFilePointer()).toCoordinateArray()[0];
            final int alignmentStart = rec.getAlignmentStart() - 1;
            final int window = LinearIndex.convertToLinearIndexOffset(alignmentStart); // the 16k window

            /* log.debug("Reference " + reference + " count " + localCount + " bin=" + bin + " startPosition=" +
                    ioffset + "(" + Long.toString(ioffset, 16) + "x)" + // " startAlignment=" +
                    // alignmentStart + "(" + Long.toString(alignmentStart,16) + "x)" + " endAlignment=" +
                    // alignmentEnd + "(" + Long.toString(alignmentEnd,16) + "x)" +
                    " window " + window);
            */

            List<LinearIndexItem> indexList = referenceToIndexEntries.get(reference);
            if (indexList == null) {
                referenceToIndexEntries.put(reference, (List<LinearIndexItem>) new ArrayList());
                indexList = referenceToIndexEntries.get(reference);
            }
            // Has this window already been seen? If so, might replace existing index
            boolean replaceIndex = true;
            for (final LinearIndexItem l : indexList) {
                if (l.window == window) {
                    if (l.offset < ioffset) {
                        replaceIndex = false;  // existing window value is already minimum
                        // log.debug("## Reference " + reference + " not adding new ioffset " + ioffset + " to window " + l.window);
                    } else {
                        // remove existingIndex in place of this one
                        //referenceToIndexEntries.get(reference).remove(l);
                        indexList.remove(l);
                        // log.debug("## Reference " + reference + " replacing old ioffset " + l.offset + " with new offset " + ioffset + " in window " + l.window);
                    }
                    break;
                }
            }
            if (replaceIndex)
                //referenceToIndexEntries.get(reference).add(new LinearIndexItem(window, ioffset));
                indexList.add(new LinearIndexItem(window, ioffset));
        }
    }

    // When all the records have been processed, finish is called
    // note we can do this processing per-reference instead of per-file if desired
    private void finish() {
        // constructs a BAMIndexContent from bins, chunks, and linear index
        for (int i = 0; i < n_ref; i++) {
            processReference(i);
        }
    }

    private void processReference(int reference) {

        // process bins. convert Set<Bin> to List<Bin> as required by BMAIndexContent
        final Set<Bin> myBins = binsInProgress.get(reference);

        // a 2 step conversion from Set -> array -> arrayList
        final int size = myBins == null ? 0 : myBins.size();
        final Bin[] bins = myBins == null ? null : new Bin[size];
        if (size != 0) {
            myBins.toArray(bins);
        }
        final List<Bin> binList = (bins == null || bins.length == 0) ? null :
                   (List<Bin>)new ArrayList(Arrays.asList(bins)) ;

        // process chunks.
        for (int j = 0; j < size; j++) {
            if (bins == null || bins[j] == null) continue;
            // List<Chunk> chunkList = binToChunks.get(bins[j]);
            final List<Chunk> oldChunks = binToSpan.get(bins[j]).getChunks();
            if (oldChunks != null) {
                // make it a modifiable list
                List<Chunk> chunks = (List<Chunk>)new ArrayList();
                for (Chunk c : oldChunks) {
                    chunks.add(c);
                }

                chunks = optimizeChunkList(chunks, 0);
                // binToChunks.put(bins[j], chunkList);
                binToSpan.put(bins[j], new BAMFileSpan(chunks));
            }
        }
        // process linear index
        final LinearIndex linearIndex = computeLinearIndex(reference);

        content[reference] = new BAMIndexContent(reference, binList, binToSpan, linearIndex);
    }

    private LinearIndex computeLinearIndex(int reference) {
        final LinearIndex linearIndex;
        final List<LinearIndexItem> indexList = referenceToIndexEntries.get(reference);
        if (indexList == null) {
        // skip linear index for this reference
        linearIndex = null;
    } else {
        // get the max window in the list; linear index will be this long
        int maxWindow = 0;   // todo - a shorter way to do this?
        for (final LinearIndexItem l : indexList) {
            if (l.window > maxWindow) {
                maxWindow = l.window;
            }
        }
        // log.debug("** n_intv (maxWindow) for reference " + i + " is " + (maxWindow + 1));
        long[] newIndex = new long[maxWindow + 1];   // newIndex is one-based
        for (int j = 0; j <= maxWindow; j++) {
            newIndex[j] = 0;
        }
        // linearIndex each entry
        for (final LinearIndexItem l : indexList) {
            newIndex[l.window] = l.offset;   // newIndex is one-based
        }
        linearIndex =  new LinearIndex(reference, 1, newIndex); // newIndex is one-based
    }
        return linearIndex;
    }

    public void writeText(boolean sortBins) throws Exception {
        final PrintWriter pw = new PrintWriter(OUTPUT);
        pw.println("n_ref=" + n_ref);
        for (int i = 0; i < n_ref; i++) {
            content[i].writeText(pw, sortBins);
        }
        pw.close();
    }
    
     public void writeBinary(final boolean sortBins, final long bamFileSize) throws Exception {

        final int bufferSize; //  = 1000000; // 1M  works, but doesn't need to be this big
        final int defaultBufferSize = 1000000;  // 1M
        if (bamFileSize < defaultBufferSize) {
            bufferSize = (int) bamFileSize;
        } else {
            bufferSize = defaultBufferSize;
        }
        // log.info("ByteBuffer size is " + bufferSize);

        final FileOutputStream stream = new FileOutputStream(OUTPUT, true);
        final FileChannel fileChannel = stream.getChannel();
        final ByteBuffer bb = ByteBuffer.allocateDirect(bufferSize);
        bb.order(ByteOrder.LITTLE_ENDIAN);

        // magic string
        final byte[] magic = BAMFileConstants.BAM_INDEX_MAGIC;
        bb.put(magic);
        bb.putInt(n_ref);
        for (int i = 0; i < n_ref; i++) {
            content[i].writeBinary(bb, sortBins);

            /* log.debug("Flipping region " + i + " at position " + bb.position()
                    //  + " channel.position=" + fileChannel.position()
                    + " file size =" + outputFile.length()); */

            //  write out data and reset the buffer for each reference
            bb.flip();
            fileChannel.write(bb, 0);
            // stream.flush();    // todo will flushing the stream at every reference be useful?
            bb.position(0);
            bb.limit(bufferSize);
        }
        bb.flip();
        fileChannel.write(bb, 0);
        fileChannel.close();
        stream.close();
    }

    private void verbose(String message) {
        boolean verbose = true;
        if (verbose) {
            //log.info (message);
            System.out.println(message);
        }
    }

        /**
        * This class facilitates a sparse representation of the linear index during construction,
        * when the size of the linear index is unknown.
        */
    public class LinearIndexItem {

        public final int window;
        public final long offset;

        public LinearIndexItem(int window, long offset) {
            this.window = window;
            this.offset = offset;
        }
}
}
