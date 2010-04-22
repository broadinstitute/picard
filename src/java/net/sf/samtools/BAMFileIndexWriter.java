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


import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * Class for writing BAM file indexes.
 */
public class BAMFileIndexWriter
      //  extends AbstractBAMFileIndex{  // gives MAX_BINS, BAM_LIDX_SHIFT, optimizeChunkList,
      //                                  // but requires close, getSearchBins
      //  extends CachingBAMFileIndex{   // gives the above plus LEVEL_STARTS, binToChunks, getFirstLocusInBin, etc.
      // extends BAMFileIndex{           //  gives stuff from AbstractBAMFileIndex,
{

    /** The number of references (chromosomes) in the BA file
     * The name n_ref corresponds to the name in the SAM Format Specification Document
     */
    public final int n_ref;

    /**
     * A mapping from bin to the chunks contained in that bin.
     */
    private final SortedMap<Bin, BAMFileSpan> binToChunks = new TreeMap<Bin, BAMFileSpan>();

    /**
     * What is the starting bin for each level?
     */
    private static final int[] LEVEL_STARTS = {0,1,9,73,585,4681};

    /**
     * A mapping of reference sequence index to set of bins. Overrides referencesToBins in parent class.
     * Differs from parent's binsInProgress since using Set<Bin> instead of Bin[]
     */
    private final SortedMap<Integer,Set<Bin>> binsInProgress = new TreeMap<Integer,Set<Bin>>();

    /**
     * A mapping of reference sequence index to linear indicies in progress
     */
    private final SortedMap<Integer,List<LinearIndexEntry>> referenceToIndexEntries =
            new TreeMap<Integer, List<LinearIndexEntry>>();

    /**
     * A mapping of reference sequence index to linear index
     */
    //private final List<LinearIndex> referenceToLinearIndex;
    private final LinearIndex[] referenceToLinearIndex;

    final Log log = Log.getInstance(getClass());

    public BAMFileIndexWriter(final int nReferences) {
        // super(file); // depending on what we subclass, we may need to pass the output bam index file
        this.n_ref = nReferences;
        // referenceToLinearIndex = new ArrayList(nReferences);
        referenceToLinearIndex = new LinearIndex[nReferences+1];
    }

    public static int createIndex(final File INPUT, final File OUTPUT, Log log, boolean createText) throws Exception{

        
        // final CachingBAMFileIndex bam = new CachingBAMFileIndex(IoUtil.openFileForReading(INPUT), true, SAMFileReader.ValidationStringency.STRICT);
        // final BAMFileReader bam = new BAMFileReader(IoUtil.openFileForReading(INPUT), true, SAMFileReader.ValidationStringency.STRICT);
        final SAMFileReader bam = new SAMFileReader(INPUT); //, SAMFileReader.ValidationStringency.STRICT);

        // check sort order in the header - must be coordinate
        if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
            throw new PicardException("Input BAM file must be sorted by coordinates");
        }

        final BAMFileIndexWriter instance = new BAMFileIndexWriter(
                bam.getFileHeader().getSequenceDictionary().size());

        int count = 0;
        int localCount = 0;

        int reference = 0;
        int lastReference = 0;

        for (Iterator<SAMRecord> iter = bam.iterator(); iter.hasNext();) {
            final SAMRecord rec = iter.next();
            if (rec.getAlignmentStart() == 0)
                continue;  // do nothing for un-aligned records

            reference = rec.getReferenceIndex();
            if (reference != lastReference) {
                // switching chromosome/reference
                lastReference = reference;
                log.info("     " + localCount + " aligned records in reference " + reference);
                // todo - could roll-up this reference
                localCount = 0;
            }
            count++;
            localCount++;
            instance.processRecord(rec, localCount);
        }
        bam.close();
        instance.finish();
        log.info("There are " + count + " aligned records in the input BAM file " + INPUT);
        if (createText){
           instance.writeText(OUTPUT);
        } else {
            instance.writeBinary(OUTPUT);
        }
        return count;
    }

    /** Record any index information for a given BAM record
     *
     * @param rec        The SAM record
     * @param localCount The count within the reference.  Used for diagnostics only
     */
    private void processRecord(final SAMRecord rec, final int localCount){

        final int bin = rec.getIndexingBin();   // todo: recompute using rec.computeIndexingBin?

        // is this bin already represented for this index?
        // if not, add it
        final int reference = rec.getReferenceIndex();
        final Bin b = new Bin(reference, bin);
        if (binsInProgress.get(reference) == null){
            binsInProgress.put(reference, new HashSet<Bin>());
        }
        binsInProgress.get(reference).add(b);

       // add chunk information
        BAMFileSpan span = binToChunks.get(b);
        if (span == null){
            span = new BAMFileSpan();
            binToChunks.put(b, span);
        }

        if (rec.getFilePointer() != null) {
           span.add((BAMFileSpan)rec.getFilePointer());// new Chunk(rec.getAlignmentStart(), rec.getAlignmentEnd()));
        } else {
            log.error("Null coordinates on record " + rec +
                    " in reference " + reference + " count " + localCount);
        }

        // add linear index information
        // the smallest file offset that starts in the 16k window for this bin
        if (bin < LEVEL_STARTS[LEVEL_STARTS.length -1]) {  // don't index the top level 6 16k bins with bin# >=4681
            final long ioffset = ((BAMFileSpan) rec.getFilePointer()).toCoordinateArray()[0];
            final int alignmentStart = rec.getAlignmentStart() - 1;
            final int window = linearOffset(alignmentStart) + 1; // the 16k window

            log.debug("Reference " + reference + " count " + localCount + " bin=" + bin + " startPosition=" +
                    ioffset + "(" + Long.toString(ioffset, 16) + "x)" + // " startAlignment=" +
                    // alignmentStart + "(" + Long.toString(alignmentStart,16) + "x)" + " endAlignment=" +
                    // alignmentEnd + "(" + Long.toString(alignmentEnd,16) + "x)" +
                    " window " + window);
            // LinearIndex index = linearIndex;
            List<LinearIndexEntry> indexList = referenceToIndexEntries.get(reference);
            if (indexList == null) {
                referenceToIndexEntries.put(reference, (List<LinearIndexEntry>) new ArrayList());
                indexList = referenceToIndexEntries.get(reference);
            }
            // Has this window already been seen?
            boolean useExistingIndex = false;
            for (final LinearIndexEntry l : indexList){
                if (l.window == window){
                    if (l.offset < ioffset){
                        useExistingIndex = true;
                        log.debug("## Reference " + reference + " not adding new ioffset " + ioffset + " to window " + l.window);
                    } else {
                        // remove existingIndex in place of this one
                        referenceToIndexEntries.get(reference).remove(l);
                       // todo?  indexList.remove(l);
                        log.debug("## Reference " + reference + " replacing old ioffset " + l.offset + " with new offset " + ioffset + " in window " + l.window);
                    }
                    break;
                }
            }
            if (!useExistingIndex)
                referenceToIndexEntries.get(reference).add(new LinearIndexEntry(window, ioffset));
        }
    }

    private int linearOffset(final long offset){
        // which 16k bin/window
        return ((int) offset >> BAM_LIDX_SHIFT );
    }

    // When all the records have been processed
    // todo - can do this processing per-reference instead of per-file if necessary
    private void finish(){
        // copy information to super class
        
        // iterate through this.binsInProgress and allocate super.referenceToBins
        for (int i = 0; i < n_ref; i++) {
            final Set<Bin> myBins = binsInProgress.get(i);
            final int size = myBins == null ? 0 : myBins.size();
            final Bin[] bins = myBins == null ? null : new Bin[size];
            if (size != 0) {
                myBins.toArray(bins); // todo possibly sort them
            }
            // referenceToBins.put(i, bins);

            // process linear index, i.e. setup referenceToLinearIndices
            List<LinearIndexEntry> indexList = referenceToIndexEntries.get(i);
            long[] newIndex = null;
            if (indexList == null){
                  // skip linear index for this reference?
                // referenceToLinearIndex.add(i, null);
                referenceToLinearIndex[i] = null;
            } else {
                // get the max window in the list
                int maxWindow = 0;
                for (final LinearIndexEntry l : indexList) {
                    if (l.window > maxWindow) {
                        maxWindow = l.window;
                    }
                }
                log.debug("** n_intv (maxWindow) for reference " + i + " is " + (maxWindow + 1));
                newIndex = new long[maxWindow + 1];
                for (int j = 0; j < maxWindow; j++) {
                    newIndex[j] = 0;
                }
                // linearIndex each entry
                for (final LinearIndexEntry l : indexList) {
                    newIndex[l.window] = l.offset;
                }
                final LinearIndex l = new LinearIndex(i, newIndex);
                // referenceToLinearIndex.add(i, l);
                referenceToLinearIndex[i] = l;
            }

            //  process chunks. Note, this must be done *after* processing the linear index
            for (int j = 0; j < size; j++) {
                BAMFileSpan chunkList = binToChunks.get(bins[j]);
                if (chunkList != null){ //  && chunkList.size() > 1) {     // todo?
                    /* todo  code from CachingBAMFileIndex - could be factored ? */
                    final int start = getFirstLocusInBin(bins[j]) - 1;
                    final int regionLinearBin = start >> BAM_LIDX_SHIFT;
                    long minimumOffset = 0;
                    if (newIndex != null  && regionLinearBin < newIndex.length)
                        minimumOffset = newIndex[regionLinearBin];
                    if (minimumOffset != 0)
                        log.debug(" Reference " + i + " optimizeChunkList minimumOffset=" + minimumOffset);
                                       // " chunkList.size()=" + chunkList.size());
                    /* debug - print chunkList before and after optimizeChunkList
                    if (i == 11)
                      Slog.debug(" Reference " + i + " Bin " + bins[j].binNumber +
                            " Before optimizeChunkList " + chunkList);  */

                    chunkList = optimizeChunkList(chunkList, minimumOffset);
                    /*
                    if (i == 11)
                      log.debug(" Reference " + i + " Bin " + bins[j].binNumber +
                            " After optimizeChunkList " + chunkList);  */
                    binToChunks.put(bins[j],chunkList);
                }
            }
        }        
    }

    public void writeText(final File outputFile) throws Exception {
        final PrintWriter pw = new PrintWriter(outputFile);
        // magic string
        // n_ref
        // these are only the ones with references final int n_ref = binsInProgress.size();
        pw.println("n_ref=" + n_ref);
        for (int i = 0; i < n_ref; i++) {
            int referenceSequence = i;
            final Set<Bin> myBins = binsInProgress.get(i);
            if (myBins == null){
                // no references to this sequence
                pw.println("Reference " + i + " has n_bin= " + 0);
                pw.println("Reference " + i + " has n_intv= " + 0);
                continue;
            }
            final int size = myBins.size();
            final Bin[] bins = new Bin[size];
            if (size != 0) {
                myBins.toArray(bins);
            }

            Arrays.sort(bins);  // todo - make this optional.  I'm just doing it to compare
            final int n_bin = bins.length;
            pw.println("Reference " + referenceSequence + " has n_bin= " + n_bin);
            for (int j = 0; j < n_bin; j++) {
                final List<Chunk> chunkList = binToChunks.get(bins[j]).chunks;
                // final BAMFileSpan chunkList = binToChunks.get(bins[j]);
                if (chunkList == null) continue;
                pw.print("  Ref " + referenceSequence + " bin " + bins[j].binNumber + " has n_chunk= " + chunkList.size());
                for (final Chunk c : chunkList) {
                    pw.println("     Chunk: " + c.toString() +
                            " start: " + Long.toString(c.getChunkStart(), 16) +
                            " end: " + Long.toString(c.getChunkEnd(), 16));
                }
            }

            //LinearIndex li = referenceToLinearIndex.get(i);
            LinearIndex li = referenceToLinearIndex[i];
            if (li != null && li.indexEntries != null) {
                // todo - is ii.indexEntries == null an internal error?
                final int n_intv = li.indexEntries.length;
                pw.println("Reference " + referenceSequence + " has n_intv= " + n_intv);
                for (int k = 0; k < li.indexEntries.length; k++) {
                    if (li.indexEntries[k] != 0) {
                        pw.println("ioffset for " + k + " is " + Long.toString(li.indexEntries[k]));
                    }
                }
            } else {
                log.debug("No linear index for reference " + referenceSequence);
                pw.println("Reference " + referenceSequence + " has n_intv= 0");
            }
        }
        pw.close();
    }

    public void writeBinary(final File outputFile) throws Exception {
        final int BAI_SIZE_IN_BYTES = 1000000; // 1M

        if (outputFile.exists()) {
            outputFile.delete(); // todo: do this in caller so can log warning?
        }
        final FileOutputStream stream = new FileOutputStream(outputFile);
        final FileChannel fileChannel = stream.getChannel();
        final ByteBuffer bb = ByteBuffer.allocateDirect(BAI_SIZE_IN_BYTES);
        bb.order(ByteOrder.LITTLE_ENDIAN);

        // magic string
        final byte[] magic = BAMFileConstants.BAM_INDEX_MAGIC;
        bb.put(magic);
        // n_ref
        bb.putInt(n_ref);
        for (int i = 0; i < n_ref; i++) {
            final Set<Bin> myBins = binsInProgress.get(i);
            if (myBins == null) {
                // no references to this sequence
                //pw.println("Reference " + i + " has n_bin= " + 0);
                //pw.println("Reference " + i + " has n_intv= " + 0);
                continue;
            }
            final int size = myBins.size();
            final Bin[] bins = new Bin[size];
            if (size != 0) {
                myBins.toArray(bins);
            }

            Arrays.sort(bins);  // todo - make this optional.  I'm just doing it to compare
            final int n_bin = bins.length;
            //pw.println("Reference " + i + " has n_bin= " + n_bin);
            bb.putInt(n_bin);
            for (int j = 0; j < n_bin; j++) {
                bb.putInt(bins[j].binNumber); // todo uint32_t vs int32_t in spec?
                final List<Chunk> chunkList = binToChunks.get(bins[j]).chunks;
                // final BAMFileSpan chunkList = binToChunks.get(bins[j]);
                if (chunkList == null) continue;  // todo - bb.putInt(0);
                // pw.print("  Ref " + referenceSequence + " bin " + bins[j].binNumber + " has n_chunk= " + chunkList.size());

                final int n_chunk = chunkList.size();
                bb.putInt(n_chunk);
                //pw.println("  Bin " + bins[j].binNumber + " has n_chunk= " + n_chunk);
                for (final Chunk c : chunkList) {
                    // pw.println("     Chunk: " + c.toString());
                    bb.putLong(c.getChunkStart());   // todo uint32_t vs int32_t in spec?
                    bb.putLong(c.getChunkEnd());     // todo uint32_t vs int32_t in spec?
                }
            }
            //LinearIndex li = referenceToLinearIndex.get(i);
            LinearIndex li = referenceToLinearIndex[i];
            if (li != null && li.indexEntries != null) {
                // todo - is li.indexEntries == null an internal error?
                final int n_intv = li.indexEntries.length;
                //pw.println("Reference " + i + " has n_intv= " + n_intv);
                bb.putInt(n_intv);
                for (int k = 0; k < li.indexEntries.length; k++) {
                    bb.putLong(li.indexEntries[k]); //// todo uint32_t vs int32_t in spec?
                }
            } else {
                log.debug("No linear index for reference " + i);
                 bb.putInt(0);
            }
            bb.flip();
            fileChannel.write(bb, 0);
            fileChannel.close();
            stream.close();
        }
    }


    // copied from AbstractBAMFileIndex
    protected static final int MAX_BINS = 37450; // =(8^6-1)/7+1
    protected static final int BAM_LIDX_SHIFT = 14;

    protected BAMFileSpan optimizeChunkList(final BAMFileSpan chunkList, final long minimumOffset) {
        Chunk lastChunk = null;
        Collections.sort(chunkList.chunks);
        final BAMFileSpan result = new BAMFileSpan();
        for (final Chunk chunk : chunkList.chunks) {
            if (chunk.getChunkEnd() <= minimumOffset) {
                continue;
            }
            if (result.isEmpty()) {
                result.add(chunk);
                lastChunk = chunk;
                continue;
            }
            // Coalesce chunks that are in adjacent file blocks.
            // This is a performance optimization.
            final long lastFileBlock = getFileBlock(lastChunk.getChunkEnd());
            final long chunkFileBlock = getFileBlock(chunk.getChunkStart());
            if (chunkFileBlock - lastFileBlock > 1) {
                result.add(chunk);
                lastChunk = chunk;
            } else {
                if (chunk.getChunkEnd() > lastChunk.getChunkEnd()) {
                    lastChunk.setChunkEnd(chunk.getChunkEnd());
                }
            }
        }
        return result;
    }

    protected long getFileBlock(final long bgzfOffset) {
        return ((bgzfOffset >> 16L) & 0xFFFFFFFFFFFFL);
    }

    // copied from CachingBAMFileIndex.java
    
    /**
     * Reports the total amount of genomic data that any bin can index.
     */
    private static final int BIN_SPAN = 512 * 1024 * 1024;

    /**
     * Get the number of levels employed by this index.
     * @return Number of levels in this index.
     */
    private int getNumIndexLevels() {
        return LEVEL_STARTS.length;
    }

    /**
     * Gets the level associated with the given bin number.
     * @param bin The bin  for which to determine the level.
     * @return the level associated with the given bin number.
     */
    private int getLevelForBin(final Bin bin) {
        if(bin.binNumber >= MAX_BINS)
            throw new SAMException("Tried to get level for invalid bin.");
        for(int i = getNumIndexLevels()-1; i >= 0; i--) {
            if(bin.binNumber >= LEVEL_STARTS[i])
                return i;
        }
        throw new SAMException("Unable to find correct bin for bin "+bin);
    }

    /**
     * Gets the first locus that this bin can index into.
     * @param bin The bin to test.
     * @return The first position that the given bin can represent.
     */
    private int getFirstLocusInBin(final Bin bin) {
        final int level = getLevelForBin(bin);
        final int levelStart = LEVEL_STARTS[level];
        final int levelSize = ((level==getNumIndexLevels()-1) ? MAX_BINS-1 : LEVEL_STARTS[level+1]) - levelStart;
        return (bin.binNumber - levelStart)*(BIN_SPAN/levelSize)+1;
    }
}
