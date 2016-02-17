package picard.sam;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.OverclippedReadFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.DelegatingIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
    private final int maxGaps;
    private final int minUnclippedBases;
    private boolean forceSort = false;
    private final OverclippedReadFilter contaminationFilter;

    /**
     * Constructor
     *
     * @param unmappedBamFile                   The BAM file that was used as the input to the aligner, which will
     *                                          include info on all the reads that did not map.  Required.
     * @param targetBamFile                     The file to which to write the merged SAM records. Required.
     * @param referenceFasta                    The reference sequence for the map files. Required.
     * @param programRecord                     Program record for taget file SAMRecords created.
     * @param clipAdapters                      Whether adapters marked in unmapped BAM file should be marked as
     *                                          soft clipped in the merged bam. Required.
     * @param bisulfiteSequence                 Whether the reads are bisulfite sequence (used when calculating the
     *                                          NM and UQ tags). Required.
     * @param alignedReadsOnly                  Whether to output only those reads that have alignment data
     * @param alignedSamFile                    The SAM file(s) with alignment information.  Optional.  If this is
     *                                          not provided, then read1AlignedSamFile and read2AlignedSamFile must be.
     * @param maxGaps                           The maximum number of insertions or deletions permitted in an
     *                                          alignment.  Alignments with more than this many gaps will be ignored.
     *                                          -1 means to allow any number of gaps.
     * @param attributesToRetain                attributes from the alignment record that should be
     *                                          removed when merging.  This overrides attributesToRetain if they share
     *                                          common tags.
     * @param read1BasesTrimmed                 The number of bases trimmed from start of read 1 prior to alignment.  Optional.
     * @param read2BasesTrimmed                 The number of bases trimmed from start of read 2 prior to alignment.  Optional.
     * @param read1AlignedSamFile               The alignment records for read1.  Used when the two ends of a read are
     *                                          aligned separately.  This is optional, but must be specified if
     *                                          alignedSamFile is not.
     * @param read2AlignedSamFile               The alignment records for read1.  Used when the two ends of a read are
     *                                          aligned separately.  This is optional, but must be specified if
     *                                          alignedSamFile is not.
     * @param expectedOrientations              A List of SamPairUtil.PairOrientations that are expected for
     *                                          aligned pairs.  Used to determine the properPair flag.
     * @param sortOrder                         The order in which the merged records should be output.  If null,
     *                                          output will be coordinate-sorted
     * @param primaryAlignmentSelectionStrategy How to handle multiple alignments for a fragment or read pair,
     *                                          in which none are primary, or more than one is marked primary
     * @param addMateCigar                      True if we are to add or maintain the mate CIGAR (MC) tag, false if we are to remove or not include.
     *
     * @param unmapContaminantReads             If true, identify reads having the signature of contamination from a foreign organism (i.e. mostly clipped bases),
     *                                          and mark them as unmapped.
     * @param minUnclippedBases                 If unmapContaminantReads is set, require this many unclipped bases or else the read will be marked as contaminant.
     */
    public SamAlignmentMerger(final File unmappedBamFile, final File targetBamFile, final File referenceFasta,
                              final SAMProgramRecord programRecord, final boolean clipAdapters, final boolean bisulfiteSequence,
                              final boolean alignedReadsOnly,
                              final List<File> alignedSamFile, final int maxGaps, final List<String> attributesToRetain,
                              final List<String> attributesToRemove,
                              final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                              final List<File> read1AlignedSamFile, final List<File> read2AlignedSamFile,
                              final List<SamPairUtil.PairOrientation> expectedOrientations,
                              final SortOrder sortOrder,
                              final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy,
                              final boolean addMateCigar,
                              final boolean unmapContaminantReads,
                              final int minUnclippedBases) {

        super(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
                alignedReadsOnly, programRecord, attributesToRetain, attributesToRemove, read1BasesTrimmed,
                read2BasesTrimmed, expectedOrientations, sortOrder, primaryAlignmentSelectionStrategy, addMateCigar, unmapContaminantReads);

        if ((alignedSamFile == null || alignedSamFile.isEmpty()) &&
                (read1AlignedSamFile == null || read1AlignedSamFile.isEmpty() ||
                        read2AlignedSamFile == null || read2AlignedSamFile.isEmpty())) {
            throw new IllegalArgumentException("Either alignedSamFile or BOTH of read1AlignedSamFile and " +
                    "read2AlignedSamFile must be specified.");
        }

        if (alignedSamFile != null) {
            for (final File f : alignedSamFile) {
                IOUtil.assertFileIsReadable(f);
            }
        } else {
            for (final File f : read1AlignedSamFile) {
                IOUtil.assertFileIsReadable(f);
            }
            for (final File f : read2AlignedSamFile) {
                IOUtil.assertFileIsReadable(f);
            }
        }

        this.alignedSamFile = alignedSamFile;
        this.read1AlignedSamFile = read1AlignedSamFile;
        this.read2AlignedSamFile = read2AlignedSamFile;
        this.maxGaps = maxGaps;
        this.minUnclippedBases = minUnclippedBases;
        this.contaminationFilter = new OverclippedReadFilter(minUnclippedBases, false);

        log.info("Processing SAM file(s): " + alignedSamFile != null ? alignedSamFile : read1AlignedSamFile + "," + read2AlignedSamFile);
    }


    /**
     * Merges the alignment from the map file with the non-aligned records from the source BAM file.
     * Overrides mergeAlignment in AbstractAlignmentMerger.  Tries first to proceed on the assumption
     * that the alignment records are pre-sorted.  If not, catches the exception, forces a sort, and
     * tries again.
     */
    public void mergeAlignment(final File referenceFasta) {
        try {
            super.mergeAlignment(referenceFasta);
        } catch (final IllegalStateException ise) {
            log.warn("Exception merging bam alignment - attempting to sort aligned reads and try again: ", ise.getMessage());
            forceSort = true;
            resetRefSeqFileWalker();
            super.mergeAlignment(referenceFasta);
        }
    }

    /**
     * Reads the aligned SAM records into a SortingCollection and returns an iterator over that collection
     */
    protected CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords() {

        final CloseableIterator<SAMRecord> mergingIterator;
        final SAMFileHeader header;

        // When the alignment records, including both ends of a pair, are in SAM files
        if (alignedSamFile != null && !alignedSamFile.isEmpty()) {
            final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(alignedSamFile.size());
            final List<SamReader> readers = new ArrayList<SamReader>(alignedSamFile.size());
            for (final File f : this.alignedSamFile) {
                final SamReader r = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(f);
                headers.add(r.getFileHeader());
                readers.add(r);

                // As we're going through and opening the aligned files, if we don't have a @PG yet
                // and there is only a single one in the input file, use that!
                if (getProgramRecord() == null && r.getFileHeader().getProgramRecords().size() == 1) {
                    setProgramRecord(r.getFileHeader().getProgramRecords().iterator().next());
                }
            }

            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SortOrder.queryname, headers, false);

            mergingIterator = new MergingSamRecordIterator(headerMerger, readers, true);
            header = headerMerger.getMergedHeader();

        }
        // When the ends are aligned separately and don't have firstOfPair information correctly
        // set we use this branch.
        else {
            mergingIterator = new SeparateEndAlignmentIterator(this.read1AlignedSamFile, this.read2AlignedSamFile, referenceFasta);
            header = ((SeparateEndAlignmentIterator) mergingIterator).getHeader();

            // As we're going through and opening the aligned files, if we don't have a @PG yet
            // and there is only a single one in the input file, use that!
            if (getProgramRecord() == null && header.getProgramRecords().size() == 1) {
                setProgramRecord(header.getProgramRecords().iterator().next());
            }
        }


        if (!forceSort) {
            return mergingIterator;
        }


        final SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
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
        return new DelegatingIterator<SAMRecord>(alignmentSorter.iterator()) {
            @Override
            public void close() {
                super.close();
                alignmentSorter.cleanup();
            }
        };
    }

    private class SuffixTrimingSamRecordIterator implements CloseableIterator<SAMRecord> {
        private final CloseableIterator<SAMRecord> underlyingIterator;
        private final String suffixToTrim;

        private SuffixTrimingSamRecordIterator(final CloseableIterator<SAMRecord> underlyingIterator, final String suffixToTrim) {
            this.underlyingIterator = underlyingIterator;
            this.suffixToTrim = suffixToTrim;
        }

        @Override
        public void close() {
            underlyingIterator.close();
        }

        @Override
        public boolean hasNext() {
            return underlyingIterator.hasNext();
        }

        @Override
        public SAMRecord next() {
            final SAMRecord rec = underlyingIterator.next();
            final String readName = rec.getReadName();
            if (readName.endsWith(suffixToTrim)) {
                rec.setReadName(readName.substring(0, readName.length() - suffixToTrim.length()));
            }
            return rec;
        }

        @Override
        public void remove() {
            underlyingIterator.remove();
        }
    }

    private class SeparateEndAlignmentIterator implements CloseableIterator<SAMRecord> {

        private final PeekableIterator<SAMRecord> read1Iterator;
        private final PeekableIterator<SAMRecord> read2Iterator;
        private final SAMFileHeader header;

        public SeparateEndAlignmentIterator(final List<File> read1Alignments, final List<File> read2Alignments, File referenceFasta) {
            final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>();
            final List<SamReader> read1 = new ArrayList<SamReader>(read1Alignments.size());
            final List<SamReader> read2 = new ArrayList<SamReader>(read2Alignments.size());
            for (final File f : read1Alignments) {
                final SamReader r = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(f);
                headers.add(r.getFileHeader());
                read1.add(r);
            }
            for (final File f : read2Alignments) {
                final SamReader r = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(f);
                headers.add(r.getFileHeader());
                read2.add(r);
            }

            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate, headers, false);
            read1Iterator = new PeekableIterator<SAMRecord>(
                    new SuffixTrimingSamRecordIterator(new MergingSamRecordIterator(headerMerger, read1, true), "/1"));
            read2Iterator = new PeekableIterator<SAMRecord>(
                    new SuffixTrimingSamRecordIterator(new MergingSamRecordIterator(headerMerger, read2, true), "/2"));

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
                } else {
                    return setPairFlags(read1Iterator.next(), true);
                }
            } else {
                return setPairFlags(read2Iterator.next(), false);
            }
        }

        public void remove() {
            throw new UnsupportedOperationException("remove() not supported");
        }

        public SAMFileHeader getHeader() { return this.header; }

        private SAMRecord setPairFlags(final SAMRecord sam, final boolean firstOfPair) {
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
    protected boolean ignoreAlignment(final SAMRecord sam) {
        if (maxGaps == -1) return false;
        int gaps = 0;
        for (final CigarElement el : sam.getCigar().getCigarElements()) {
            if (el.getOperator() == CigarOperator.I || el.getOperator() == CigarOperator.D) {
                gaps++;
            }
        }
        return gaps > maxGaps;
    }

    /**
     * Criteria for contaminant reads:
     * 1. primary alignment has fewer than minUnclippedBases unclipped bases
     * 2. primary alignment has both ends clipped
     * 3. for pairs, at least one end of primary alignment meets above criteria
     */
    protected boolean isContaminant(final HitsForInsert hits) {
        boolean isContaminant = false;
        if (hits.numHits() > 0) {
            final int primaryIndex = hits.getIndexOfEarliestPrimary();
            if (primaryIndex < 0) throw new IllegalStateException("No primary alignment was found, despite having nonzero hits.");
            final SAMRecord primaryRead1 = hits.getFirstOfPair(primaryIndex);
            final SAMRecord primaryRead2 = hits.getSecondOfPair(primaryIndex);
            if (primaryRead1 != null && primaryRead2 != null) isContaminant = contaminationFilter.filterOut(primaryRead1, primaryRead2);
            else if (primaryRead1 != null) isContaminant = contaminationFilter.filterOut(primaryRead1);
            else if (primaryRead2 != null) isContaminant = contaminationFilter.filterOut(primaryRead2);
            else throw new IllegalStateException("Neither read1 or read2 exist for chosen primary alignment");
        }
        return isContaminant;
    }

    // Accessor for testing
    public boolean getForceSort() {return this.forceSort; }
}
