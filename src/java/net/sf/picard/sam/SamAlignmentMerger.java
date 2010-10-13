package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class that takes in a set of alignment information in SAM format and merges it with the set
 * of all reads for which alignment was attempted, stored in an unmapped SAM file.
 *
 * @author ktibbett@broadinstitute.org
 */
public class SamAlignmentMerger extends AbstractAlignmentMerger {

    private final Log log = Log.getInstance(SamAlignmentMerger.class);
    private final File alignedSamFile;
    private final boolean pairedRun;
    private final int maxGaps;
    private SAMFileReader reader;
    private boolean forceSort = false;

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map
     * @param targetBamFile     The file to which to write the merged SAM records
     * @param referenceFasta    The reference sequence for the map files
     * @param programRecord     Program record for taget file SAMRecords created.
     * @param clipAdapters      Whether adapters marked in BAM file are clipped from the read
     * @param bisulfiteSequence Whether the reads are bisulfite sequence
     * @param pairedRun         Whether the run is a paired-end run
     * @param jumpingLibrary    Whether this is a jumping library
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param alignedSamFile    The SAM file with alignment information
     * @param maxGaps           The maximum number of insertions or deletions permitted in an
     *                          alignment.  Alignments with more than this many gaps will be ignored.
     * @param attributesToRetain  private attributes from the alignment record that should be
     *                          included when merging.  This overrides the exclusion of
     *                          attributes whose tags start with the reserved characters
     *                          of X, Y, and Z
     *
     */
    public SamAlignmentMerger (final File unmappedBamFile, final File targetBamFile, final File referenceFasta,
                 final SAMProgramRecord programRecord, final boolean clipAdapters, final boolean bisulfiteSequence,
                 final boolean pairedRun, final boolean jumpingLibrary, final boolean alignedReadsOnly,
                 final File alignedSamFile, final int maxGaps, final List<String> attributesToRetain) {

        super(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
              jumpingLibrary, alignedReadsOnly, programRecord, attributesToRetain);

        IoUtil.assertFileIsReadable(alignedSamFile);
        this.alignedSamFile = alignedSamFile;
        this.pairedRun = pairedRun;
        this.maxGaps = maxGaps;
        this.reader = new SAMFileReader(this.alignedSamFile);
        this.reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
        if (programRecord == null) {
            if (reader.getFileHeader().getProgramRecords().size() == 1) {
                setProgramRecord(reader.getFileHeader().getProgramRecords().get(0));
            }
        }
        // If not null, the program record was already added in the superclass.  DO NOT RE-ADD!

        if (getProgramRecord() != null) {
            SAMFileReader tmp = new SAMFileReader(unmappedBamFile);
            try {
                for (SAMProgramRecord pg : tmp.getFileHeader().getProgramRecords()) {
                    if (pg.getId().equals(getProgramRecord().getId())) {
                        throw new PicardException("Program Record ID already in use in unmapped BAM file.");
                    }
                }
            }
            finally {
                tmp.close();
            }
        }
        log.info("Processing SAM file: " + alignedSamFile.getAbsolutePath());
    }

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map
     * @param targetBamFile     The file to which to write the merged SAM records
     * @param referenceFasta    The reference sequence for the map files
     * @param programRecord     Program record for taget file SAMRecords created.
     * @param clipAdapters      Whether adapters marked in BAM file are clipped from the read
     * @param bisulfiteSequence Whether the reads are bisulfite sequence
     * @param pairedRun         Whether the run is a paired-end run
     * @param jumpingLibrary    Whether this is a jumping library
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param alignedSamFile    The SAM file with alignment information
     * @param maxGaps           The maximum number of insertions or deletions permitted in an
     *                          alignment.  Alignments with more than this many gaps will be ignored.
     */
    public SamAlignmentMerger (final File unmappedBamFile, final File targetBamFile, final File referenceFasta,
                 final SAMProgramRecord programRecord, final boolean clipAdapters, final boolean bisulfiteSequence,
                 final boolean pairedRun, final boolean jumpingLibrary, final boolean alignedReadsOnly,
                 final File alignedSamFile, final int maxGaps) {

        this(unmappedBamFile, targetBamFile, referenceFasta, programRecord, clipAdapters, bisulfiteSequence,
              pairedRun, jumpingLibrary, alignedReadsOnly, alignedSamFile, maxGaps, null);
    }


    /**
     * Merges the alignment from the map file with the non-aligned records from the source BAM file.
     */
    public void mergeAlignment() {
        try {
            super.mergeAlignment();
        }
        catch(Exception e) {
            forceSort = true;
            setRefSeqFileWalker();
            super.mergeAlignment();
        }
    }

    /**
     * Reads the aligned SAM records into a SortingCollection and returns an iterator over that collection
     */
    protected CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords() {
        final SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        header.setSequenceDictionary(getSequenceDictionary());
        if (getProgramRecord() != null) {
            header.addProgramRecord(getProgramRecord());
        }

        if (!forceSort) {
            return new SortedAlignmentIterator(reader.iterator());
        }
        this.reader = new SAMFileReader(this.alignedSamFile);

        SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                    new BAMRecordCodec(header), new SAMRecordQueryNameComparator(), MAX_RECORDS_IN_RAM);

        Map<String, SAMRecord> readNameToReadPending = new HashMap<String,SAMRecord>();
        int count = 0;
        for (CloseableIterator<SAMRecord> it = new CigarClippingIterator(reader.iterator()); it.hasNext();) {
            SAMRecord sam = it.next();
            sam.setReadName(cleanReadName(sam.getReadName()));
            if (pairedRun) {
                SAMRecord mate = readNameToReadPending.remove(sam.getReadName());
                if (mate != null) {
                    if ((!sam.getReadUnmappedFlag()) || (!mate.getReadUnmappedFlag())) {
                        fixMates(sam, mate);
                        alignmentSorter.add(sam);
                        alignmentSorter.add(mate);
                        count += 2;
                    }
                }
                else {
                    readNameToReadPending.put(sam.getReadName(), sam);
                }
            }
            else {
                if (!sam.getReadUnmappedFlag()) {
                    alignmentSorter.add(sam);
                    count++;
                }
            }
            if (count > 0 && count % 1000000 == 0) {
                log.info("Read " + count + " records from alignment SAM/BAM.");
            }
        }
        log.info("Finished reading " + count + " total records from alignment SAM/BAM.");

        if (readNameToReadPending.size() > 0) {
            throw new PicardException("Unmatched reads left in pending map.");
        }
        
        reader.close();
        return alignmentSorter.iterator();
    }


    private class SortedAlignmentIterator implements CloseableIterator<SAMRecord> {

        CloseableIterator<SAMRecord> wrappedIterator;
        private SAMRecord next = null;

        public SortedAlignmentIterator(SAMRecordIterator it) {
            wrappedIterator = new CigarClippingIterator(it);
            it.assertSorted(SAMFileHeader.SortOrder.queryname);
        }

        public void close() {
            wrappedIterator.close();
        }

        public boolean hasNext() {
            return wrappedIterator.hasNext() || next != null;
        }

        public SAMRecord next() {

            if (!pairedRun) {
                return wrappedIterator.next();
            }

            if (next == null) {
                SAMRecord firstOfPair = wrappedIterator.next();
                next = wrappedIterator.next();
                fixMates(firstOfPair, next);
                firstOfPair.setReadName(cleanReadName(firstOfPair.getReadName()));
                return firstOfPair;
            }
            else {
                SAMRecord secondOfPair = next;
                secondOfPair.setReadName(cleanReadName(secondOfPair.getReadName()));
                next = null;
                return secondOfPair;
            }
        }

        public void remove() {
            throw new UnsupportedOperationException("remove() not supported");
        }
    }

    private void fixMates(SAMRecord first, SAMRecord second) {
        if ((!first.getReadUnmappedFlag()) || (!second.getReadUnmappedFlag())) {
            final boolean proper = SamPairUtil.isProperPair(first, second, isJumpingLibrary());
            first.setProperPairFlag(proper);
            second.setProperPairFlag(proper);
            SamPairUtil.setMateInfo(first, second, getHeader());
        }
    }

    /**
     * For now, we only ignore those alignments that have more than <code>maxGaps</code> insertions
     * or deletions.
     */
    protected boolean ignoreAlignment(SAMRecord sam) {
        int gaps = 0;
        for (CigarElement el : sam.getCigar().getCigarElements()) {
            if (el.getOperator() == CigarOperator.I || el.getOperator() == CigarOperator.D ||
                el.getOperator() == CigarOperator.N) {
                gaps++;
            }
        }
        return gaps > maxGaps;
    }

    // Accessor for testing
    public boolean getForceSort() {return this.forceSort; }
}
