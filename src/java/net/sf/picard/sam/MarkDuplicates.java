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

package net.sf.picard.sam;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.util.Log;
import net.sf.picard.PicardException;
import net.sf.samtools.*;
import net.sf.samtools.util.SortingCollection;
import net.sf.samtools.util.SortingLongCollection;

import java.io.*;
import java.util.*;

/**
 * A better duplication marking algorithm that handles all cases including clipped
 * and gapped alignments.
 *
 * @author Tim Fennell
 */
public class MarkDuplicates extends CommandLineProgram {
    private static final Log log = Log.getInstance(MarkDuplicates.class);


    /**
     * If more than this many sequences in SAM file, don't spill to disk because there will not
     * be enough file handles.
     */
    private static final int MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP = 500;

    @Usage public final String USAGE =
            "Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. " +
            "All records are then written to the output file with the duplicate records flagged.";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file to analyze.  Must be coordinate sorted.") public File INPUT;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file to right marked records to") public File OUTPUT;
    @Option(shortName="M", doc="File to write duplication metrics to") public File METRICS_FILE;

    private SortingCollection<ReadEnds> pairSort;
    private SortingCollection<ReadEnds> fragSort;
    private SortingLongCollection duplicateIndexes;
    private int numDuplicateIndices = 0;
    private int nextIndex = 0; // The next offset into duplicateIndexes to use


    /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new MarkDuplicates().instanceMain(args));
    }

    /**
     * Main work method.  Reads the BAM file once and collects sorted information about
     * the 5' ends of both ends of each read (or just one end in the case of pairs).
     * Then makes a pass through those determining duplicates before re-reading the
     * input file and writing it out with duplication flags set correctly.
     */
    protected int doWork() {
        reportMemoryStats("Start of doWork");
        log.info("Reading input file and constructing read end information.");
        buildSortedReadEndLists();
        reportMemoryStats("After buildSortedReadEndLists");
        generateDuplicateIndexes();
        reportMemoryStats("After generateDuplicateIndexes");
        log.info("Marking " + this.numDuplicateIndices + " records as duplicates.");
        final DuplicationMetrics metrics = new DuplicationMetrics();
        final SAMFileReader in  = new SAMFileReader(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(),
                                                                          true,
                                                                          OUTPUT);

        // Now copy over the file while marking all the necessary indexes as duplicates
        long recordInFileIndex = 0;
        long nextDuplicateIndex = (this.duplicateIndexes.hasNext() ? this.duplicateIndexes.next(): -1);
        int  arrayIndex = 1;

        for (final SAMRecord rec : in) {
            // First bring the simple metrics up to date
            if (rec.getReadUnmappedFlag()) {
                ++metrics.UNMAPPED_READS;
            }
            else if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                ++metrics.UNPAIRED_READS_EXAMINED;
            }
            else if (rec.getFirstOfPairFlag()){
                ++metrics.READ_PAIRS_EXAMINED;
            }


            if (recordInFileIndex++ == nextDuplicateIndex) {
                rec.setDuplicateReadFlag(true);

                // Update the duplication metrics
                if (!rec.getReadPairedFlag() || rec.getMateUnmappedFlag()) {
                    ++metrics.UNPAIRED_READ_DUPLICATES;
                }
                else if (rec.getFirstOfPairFlag()) {
                    ++metrics.READ_PAIR_DUPLICATES;
                }

                // Now try and figure out the next duplicate index
                if (this.duplicateIndexes.hasNext()) {
                    nextDuplicateIndex = this.duplicateIndexes.next();
                } else {
                    // Only happens once we've marked all the duplicates
                    nextDuplicateIndex = -1;
                    arrayIndex = -1;
                }
            }
            else {
                rec.setDuplicateReadFlag(false);
            }

            out.addAlignment(rec);
        }

        reportMemoryStats("Before output close");
        out.close();
        reportMemoryStats("After output close");


        // Write out the metrics
        metrics.calculateDerivedMetrics();
        final MetricsFile<DuplicationMetrics,Double> file = getMetricsFile();
        file.addMetric(metrics);
        file.setHistogram(metrics.calculateRoiHistogram());
        file.write(METRICS_FILE);        

        return 0;
    }

    private void reportMemoryStats(final String stage) {
        System.gc();
        final Runtime runtime = Runtime.getRuntime();
        log.info(stage + " freeMemory: " + runtime.freeMemory() + "; totalMemory: " + runtime.totalMemory() +
                "; maxMemory: " + runtime.maxMemory());
    }

    /**
     * Goes through all the records in a file and generates a set of ReadEnds objects that
     * hold the necessary information (reference sequence, 5' read coordinate) to do
     * duplication, caching to disk as necssary to sort them.
     */
    private void buildSortedReadEndLists() {
        final int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * 0.25) / ReadEnds.SIZE_OF);
        log.info("Will retain up to " + maxInMemory + " data points before spilling to disk.");

        this.pairSort = SortingCollection.newInstance(ReadEnds.class,
                                                      new ReadEndsCodec(),
                                                      new ReadEndsComparator(),
                                                      maxInMemory);

        this.fragSort = SortingCollection.newInstance(ReadEnds.class,
                                                      new ReadEndsCodec(),
                                                      new ReadEndsComparator(),
                                                      maxInMemory);

        final SAMFileReader sam = new SAMFileReader(INPUT);
        final SAMFileHeader header = sam.getFileHeader();
        if (header.getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException(INPUT + " is not coordinate sorted.");
        }
        final ReadEndsMap tmp;
        if (header.getSequenceDictionary().getSequences().size() > MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP) {
            tmp = new RAMReadEndsMap();
        } else {
            tmp = new DiskReadEndsMap();
        }
        long index = 0;

        for (final SAMRecord rec : sam) {
            if (rec.getReadUnmappedFlag()) {
                if (rec.getReferenceIndex() == -1) {
                    // When we hit the unmapped reads with no coordinate, no reason to continue.
                    break;
                }
                // If this read is unmapped but sorted with the mapped reads, just skip it.
            }
            else {
                final ReadEnds fragmentEnd = buildReadEnds(header, index, rec);
                this.fragSort.add(fragmentEnd);

                if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
                    final String key = rec.getAttribute(ReservedTagConstants.READ_GROUP_ID) + ":" + rec.getReadName();
                    ReadEnds pairedEnds = tmp.remove(rec.getReferenceIndex(), key);

                    // See if we've already seen the first end or not
                    if (pairedEnds == null) {
                        pairedEnds = buildReadEnds(header, index, rec);
                        tmp.put(pairedEnds.read2Sequence, key, pairedEnds);
                    }
                    else {
                        final int sequence = fragmentEnd.read1Sequence;
                        final int coordinate = fragmentEnd.read1Coordinate;

                        // If the second read is actually later, just add the second read data, else flip the reads
                        if (sequence > pairedEnds.read1Sequence ||
                                (sequence == pairedEnds.read1Sequence && coordinate >= pairedEnds.read1Coordinate)) {
                            pairedEnds.read2Sequence    = sequence;
                            pairedEnds.read2Coordinate  = coordinate;
                            pairedEnds.read2IndexInFile = index;
                            pairedEnds.orientation = getOrientationByte(pairedEnds.orientation == ReadEnds.R,
                                                                        rec.getReadNegativeStrandFlag());
                        }
                        else {
                            pairedEnds.read2Sequence    = pairedEnds.read1Sequence;
                            pairedEnds.read2Coordinate  = pairedEnds.read1Coordinate;
                            pairedEnds.read2IndexInFile = pairedEnds.read1IndexInFile;
                            pairedEnds.read1Sequence    = sequence;
                            pairedEnds.read1Coordinate  = coordinate;
                            pairedEnds.read1IndexInFile = index;
                            pairedEnds.orientation = getOrientationByte(rec.getReadNegativeStrandFlag(),
                                                                        pairedEnds.orientation == ReadEnds.R);
                        }

                        pairedEnds.score += getScore(rec);
                        this.pairSort.add(pairedEnds);
                    }
                }
            }

            // Print out some stats every 1m reads
            if (++index % 1000000 == 0) {
                log.info("Read " + index + " records. Tracking " + tmp.size() + " as yet unmatched pairs. " +
                tmp.sizeInRam() + " records in RAM.  Last sequence index: " + rec.getReferenceIndex());
            }
        }

        log.info("Read " + index + " records. " + tmp.size() + " pairs never matched.");
        sam.close();

        // Tell these collections to free up memory if possible.
        this.pairSort.doneAdding();
        this.fragSort.doneAdding();
    }

    /** Builds a read ends object that represents a single read. */
    private ReadEnds buildReadEnds(final SAMFileHeader header, final long index, final SAMRecord rec) {
        final ReadEnds ends = new ReadEnds();
        ends.read1Sequence    = rec.getReferenceIndex();
        ends.read1Coordinate  = rec.getReadNegativeStrandFlag() ? rec.getUnclippedEnd() : rec.getUnclippedStart();
        ends.orientation = rec.getReadNegativeStrandFlag() ? ReadEnds.R : ReadEnds.F;
        ends.read1IndexInFile = index;
        ends.score = getScore(rec);

        // Doing this lets the ends object know that it's part of a pair
        if (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag()) {
            ends.read2Sequence = rec.getMateReferenceIndex();
        }

        return ends;
    }

    /**
     * Returns a single byte that encodes the orientation of the two reads in a pair.
     */
    private byte getOrientationByte(final boolean read1NegativeStrand, final boolean read2NegativeStrand) {
        if (read1NegativeStrand) {
            if (read2NegativeStrand)  return ReadEnds.RR;
            else return ReadEnds.RF;
        }
        else {
            if (read2NegativeStrand)  return ReadEnds.FR;
            else return ReadEnds.FF;
        }
    }



    /** Calculates a score for the read which is the sum of scores over Q20. */
    private short getScore(final SAMRecord rec) {
        short score = 0;
        for (final byte b : rec.getBaseQualities()) {
            if (b >= 15) score += b;
        }

        return score;
    }

    /**
     * Goes through the accumulated ReadEnds objects and determines which of them are
     * to be marked as duplicates.
     *
     * @return an array with an ordered list of indexes into the source file
     */
    private void generateDuplicateIndexes() {
        final int maxInMemory = (int) ((Runtime.getRuntime().maxMemory() * 0.25) / SortingLongCollection.SIZEOF);
        log.info("Will retain up to " + maxInMemory + " duplicate indices before spilling to disk.");
        this.duplicateIndexes = new SortingLongCollection(maxInMemory, TMP_DIR);

        ReadEnds firstOfNextChunk = null;
        final List<ReadEnds> nextChunk  = new ArrayList<ReadEnds>(200);

        // First just do the pairs
        log.info("Traversing read pair information and detecting duplicates.");
        for (final ReadEnds next : this.pairSort) {
            if (firstOfNextChunk == null) {
                firstOfNextChunk = next;
                nextChunk.add(firstOfNextChunk);
            }
            else if (areComparableForDuplicates(firstOfNextChunk, next, true)) {
                nextChunk.add(next);
            }
            else {
                if (nextChunk.size() > 1) {
                    markDuplicatePairs(nextChunk);
                }

                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
            }
        }
        markDuplicatePairs(nextChunk);
        this.pairSort = null;

        // Now deal with the fragments
        log.info("Traversing fragment information and detecting duplicates.");
        boolean containsPairs = false;
        boolean containsFrags = false;

        for (final ReadEnds next : this.fragSort) {
            if (firstOfNextChunk != null && areComparableForDuplicates(firstOfNextChunk, next, false)) {
                nextChunk.add(next);
                containsPairs = containsPairs || next.isPaired();
                containsFrags = containsFrags || !next.isPaired();
            }
            else {
                if (nextChunk.size() > 1 && containsFrags) {
                    markDuplicateFragments(nextChunk, containsPairs);
                }

                nextChunk.clear();
                nextChunk.add(next);
                firstOfNextChunk = next;
                containsPairs = next.isPaired();
                containsFrags = !next.isPaired();
            }
        }
        markDuplicateFragments(nextChunk, containsPairs);
        this.fragSort = null;

        log.info("Sorting list of duplicate records.");
        this.duplicateIndexes.doneAddingStartIteration();
    }

    private boolean areComparableForDuplicates(final ReadEnds lhs, final ReadEnds rhs, final boolean compareRead2) {
        boolean retval =  lhs.read1Sequence   == rhs.read1Sequence &&
                          lhs.read1Coordinate == rhs.read1Coordinate &&
                          lhs.orientation     == rhs.orientation;

        if (retval && compareRead2) {
            retval = lhs.read2Sequence   == rhs.read2Sequence &&
                     lhs.read2Coordinate == rhs.read2Coordinate;
        }

        return retval;
    }

    private void addIndexAsDuplicate(final long bamIndex) {
        this.duplicateIndexes.add(bamIndex);
        ++this.numDuplicateIndices;
    }

    /**
     * Takes a list of ReadEnds objects and removes from it all objects that should
     * not be marked as duplicates.
     *
     * @param list
     */
    private void markDuplicatePairs(final List<ReadEnds> list) {
        short maxScore = 0;
        ReadEnds best = null;

        for (final ReadEnds end : list) {
            if (end.score > maxScore || best == null) {
                maxScore = end.score;
                best = end;
            }
        }

        for (final ReadEnds end : list) {
            if (end != best) {
                addIndexAsDuplicate(end.read1IndexInFile);
                addIndexAsDuplicate(end.read2IndexInFile);
            }
        }
    }

    /**
     * Takes a list of ReadEnds objects and removes from it all objects that should
     * not be marked as duplicates.
     *
     * @param list
     */
    private void markDuplicateFragments(final List<ReadEnds> list, final boolean containsPairs) {
        if (containsPairs) {
            for (final ReadEnds end : list) {
                if (!end.isPaired()) addIndexAsDuplicate(end.read1IndexInFile);
            }
        }
        else {
            short maxScore = 0;
            ReadEnds best = null;
            for (final ReadEnds end : list) {
                if (end.score > maxScore || best == null) {
                    maxScore = end.score;
                    best = end;
                }
            }

            for (final ReadEnds end : list) {
                if (end != best) {
                    addIndexAsDuplicate(end.read1IndexInFile);
                }
            }
        }
    }

    /** Comparator for ReadEnds that orders by read1 position then pair orientation then read2 position. */
    static class ReadEndsComparator implements Comparator<ReadEnds> {
        public int compare(final ReadEnds lhs, final ReadEnds rhs) {
            int retval = lhs.read1Sequence - rhs.read1Sequence;
            if (retval == 0) retval = lhs.read1Coordinate - rhs.read1Coordinate;
            if (retval == 0) retval = lhs.orientation - rhs.orientation;
            if (retval == 0) retval = lhs.read2Sequence   - rhs.read2Sequence;
            if (retval == 0) retval = lhs.read2Coordinate - rhs.read2Coordinate;
            if (retval == 0) retval = (int) (lhs.read1IndexInFile - rhs.read1IndexInFile);
            if (retval == 0) retval = (int) (lhs.read2IndexInFile - rhs.read2IndexInFile);

            return retval;
        }
    }
}
