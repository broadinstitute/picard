package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SortingCollection;

import java.io.File;
import java.util.*;

/**
 * Class that takes in a set of alignment information in SAM format and merges it with the set
 * of all reads for which alignment was attempted, stored in an unmapped SAM file.  This
 * class overrides mergeAlignment in AbstractAlignmentNMerger and proceeds on the assumption that
 * the underlying alignment records are aleady in query-name sorted order (true for bwa).  If
 * they are not, the mergeAlignment method catches the IllegalStateException, forces a sort
 * of the underlying alignment records, and tries again.
 *
 * @author ktibbett@broadinstitute.org
 */
public class SamAlignmentMerger extends AbstractAlignmentMerger {

    private final Log log = Log.getInstance(SamAlignmentMerger.class);
    private final List<File> alignedSamFile;
    private final List<File> read1AlignedSamFile;
    private final List<File> read2AlignedSamFile;
    private final boolean pairedRun;
    private final int maxGaps;
    private boolean forceSort = false;

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the aligner, which will
     *                          include info on all the reads that did not map.  Required.
     * @param targetBamFile     The file to which to write the merged SAM records. Required.
     * @param referenceFasta    The reference sequence for the map files. Required.
     * @param programRecord     Program record for taget file SAMRecords created.
     * @param clipAdapters      Whether adapters marked in unmapped BAM file should be marked as
     *                          soft clipped in the merged bam. Required.
     * @param bisulfiteSequence Whether the reads are bisulfite sequence (used when calculating the
     *                          NM and UQ tags). Required.
     * @param pairedRun           Whether the run is a paired-end run. Required.
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param alignedSamFile      The SAM file(s) with alignment information.  Optional.  If this is
     *                            not provided, then read1AlignedSamFile and read2AlignedSamFile must be.
     * @param maxGaps             The maximum number of insertions or deletions permitted in an
     *                            alignment.  Alignments with more than this many gaps will be ignored.
     * @param attributesToRetain  private attributes from the alignment record that should be
     *                          included when merging.  This overrides the exclusion of
     *                          attributes whose tags start with the reserved characters
     *                          of X, Y, and Z
     * @param read1BasesTrimmed The number of bases trimmed from start of read 1 prior to alignment.  Optional.
     * @param read2BasesTrimmed The number of bases trimmed from start of read 2 prior to alignment.  Optional.
     * @param read1AlignedSamFile The alignment records for read1.  Used when the two ends of a read are
     *                            aligned separately.  This is optional, but must be specified if
     *                            alignedSamFile is not.
     * @param read2AlignedSamFile The alignment records for read1.  Used when the two ends of a read are
     *                            aligned separately.  This is optional, but must be specified if
     *                            alignedSamFile is not.
     * @param expectedOrientations A List of SamPairUtil.PairOrientations that are expected for
     *                          aligned pairs.  Used to determine the properPair flag.
     *
     */
    public SamAlignmentMerger (final File unmappedBamFile, final File targetBamFile, final File referenceFasta,
                 final SAMProgramRecord programRecord, final boolean clipAdapters, final boolean bisulfiteSequence,
                 final boolean pairedRun, final boolean alignedReadsOnly,
                 final List<File> alignedSamFile, final int maxGaps, final List<String> attributesToRetain,
                 final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                 final List<File> read1AlignedSamFile, final List<File> read2AlignedSamFile,
                 final List<SamPairUtil.PairOrientation> expectedOrientations) {

        super(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
              alignedReadsOnly, programRecord, attributesToRetain, read1BasesTrimmed,
              read2BasesTrimmed, expectedOrientations);

        if ((alignedSamFile == null || alignedSamFile.size() == 0) &&
            (read1AlignedSamFile == null || read1AlignedSamFile.size() == 0 ||
             read2AlignedSamFile == null || read2AlignedSamFile.size() == 0)) {
            throw new IllegalArgumentException("Either alignedSamFile or BOTH of read1AlignedSamFile and " +
                    "read2AlignedSamFile must be specified.");
        }

        if (alignedSamFile != null) {
            for (File f : alignedSamFile) {
                IoUtil.assertFileIsReadable(f);
            }
        }
        else {
            for (File f : read1AlignedSamFile) {
                IoUtil.assertFileIsReadable(f);
            }
            for (File f : read2AlignedSamFile) {
                IoUtil.assertFileIsReadable(f);
            }
        }

        this.alignedSamFile = alignedSamFile;
        this.read1AlignedSamFile = read1AlignedSamFile;
        this.read2AlignedSamFile = read2AlignedSamFile;
        this.pairedRun = pairedRun;
        this.maxGaps = maxGaps;
        if (programRecord == null) {
            File tmpFile = this.alignedSamFile != null && this.alignedSamFile.size() > 0
                    ? this.alignedSamFile.get(0)
                    : this.read1AlignedSamFile.get(0);
            SAMFileReader tmpReader = new SAMFileReader(tmpFile);
            tmpReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
            if (tmpReader.getFileHeader().getProgramRecords().size() == 1) {
                setProgramRecord(tmpReader.getFileHeader().getProgramRecords().get(0));
            }
            tmpReader.close();
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
        log.info("Processing SAM file(s): " + alignedSamFile != null ? alignedSamFile : read1AlignedSamFile + "," + read2AlignedSamFile);
    }

    /**
     * Merges the alignment from the map file with the non-aligned records from the source BAM file.
     * Overrides mergeAlignment in AbstractAlignmentMerger.  Tries first to proceed on the assumption
     * that the alignment records are pre-sorted.  If not, catches the exception, forces a sort, and
     * tries again.
     */
    public void mergeAlignment() {
        try {
            super.mergeAlignment();
        }
        catch(IllegalStateException e) {
            forceSort = true;
            resetRefSeqFileWalker();
            super.mergeAlignment();
        }
    }

    /**
     * Reads the aligned SAM records into a SortingCollection and returns an iterator over that collection
     */
    protected CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords() {

        CloseableIterator<SAMRecord> mergingIterator;
        SAMFileHeader header;

        // When the alignment records, including both ends of a pair, are in SAM files
        if (alignedSamFile != null && alignedSamFile.size() > 0) {
            List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(alignedSamFile.size());
            List<SAMFileReader> readers = new ArrayList<SAMFileReader>(alignedSamFile.size());
            for (File f : this.alignedSamFile) {
                SAMFileReader r = new SAMFileReader(f);
                headers.add(r.getFileHeader());
                readers.add(r);
            }

            SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);

            mergingIterator = new MergingSamRecordIterator(headerMerger, readers, true);
            header = headerMerger.getMergedHeader();

        }
        // When the ends are aligned separately and don't have firstOfPair information correctly
        // set we use this branch.
        else {
            mergingIterator = new SeparateEndAlignmentIterator(this.read1AlignedSamFile, this.read2AlignedSamFile);
            header = ((SeparateEndAlignmentIterator)mergingIterator).getHeader();
        }


        if (!forceSort) {
            return mergingIterator;
        }


        SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                    new BAMRecordCodec(header), new SAMRecordQueryNameComparator(), MAX_RECORDS_IN_RAM);

        int count = 0;
        while (mergingIterator.hasNext()) {
            alignmentSorter.add(mergingIterator.next());
            count++;
            if (count > 0 && count % 1000000 == 0) {
                log.info("Read " + count + " records from alignment SAM/BAM.");
            }
        }
        log.info("Finished reading " + count + " total records from alignment SAM/BAM.");

        mergingIterator.close();
        return alignmentSorter.iterator();
    }

    private class SeparateEndAlignmentIterator implements CloseableIterator<SAMRecord> {

        private final PeekableIterator<SAMRecord> read1Iterator;
        private final PeekableIterator<SAMRecord> read2Iterator;
        private final SAMFileHeader header;

        public SeparateEndAlignmentIterator(List<File> read1Alignments, List<File> read2Alignments) {
            final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
            final List<SAMFileReader> read1 = new ArrayList<SAMFileReader>(read1Alignments.size());
            final List<SAMFileReader> read2 = new ArrayList<SAMFileReader>(read2Alignments.size());
            for (File f : read1Alignments) {
                SAMFileReader r = new SAMFileReader(f);
                headers.add(r.getFileHeader());
                read1.add(r);
            }
            for (File f : read2Alignments) {
                SAMFileReader r = new SAMFileReader(f);
                headers.add(r.getFileHeader());
                read2.add(r);
            }

            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);
            read1Iterator = new PeekableIterator(new MergingSamRecordIterator(headerMerger, read1, true));
            read2Iterator = new PeekableIterator(new MergingSamRecordIterator(headerMerger, read2, true));

            header = headerMerger.getMergedHeader();
        }

        public void close() {
            read1Iterator.close();
            read2Iterator.close();
        }

        public boolean hasNext() {
            return read1Iterator.hasNext() || read2Iterator.hasNext();
        }

        public SAMRecord next() {
            if (read1Iterator.hasNext()) {
                if (read2Iterator.hasNext()) {
                    return (read1Iterator.peek().getReadName().compareTo(read2Iterator.peek().getReadName()) <= 0)
                        ? setPairFlags(read1Iterator.next(), true)
                        : setPairFlags(read2Iterator.next(), false);
                }
                else {
                    return setPairFlags(read1Iterator.next(), true);
                }
            }
            else {
                return setPairFlags(read2Iterator.next(), false);
            }
        }

        public void remove() {
            throw new UnsupportedOperationException("remove() not supported");
        }

        public SAMFileHeader getHeader() { return this.header; }

        private SAMRecord setPairFlags(SAMRecord sam, boolean firstOfPair) {
            sam.setReadPairedFlag(true);
            sam.setFirstOfPairFlag(firstOfPair);
            sam.setSecondOfPairFlag(!firstOfPair);
            return sam;
        }
    }

    /**
     * For now, we only ignore those alignments that have more than <code>maxGaps</code> insertions
     * or deletions.
     */
    protected boolean ignoreAlignment(SAMRecord sam) {
        int gaps = 0;
        for (CigarElement el : sam.getCigar().getCigarElements()) {
            if (el.getOperator() == CigarOperator.I || el.getOperator() == CigarOperator.D ) {
                gaps++;
            }
        }
        return gaps > maxGaps;
    }

    // Accessor for testing
    public boolean getForceSort() {return this.forceSort; }
}
