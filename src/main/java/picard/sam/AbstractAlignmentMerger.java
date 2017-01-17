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
package picard.sam;

import htsjdk.samtools.*;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.*;
import picard.PicardException;

import java.io.File;
import java.util.*;

/**
 * Abstract class that coordinates the general task of taking in a set of alignment information,
 * possibly in SAM format, possibly in other formats, and merging that with the set of all reads
 * for which alignment was attempted, stored in an unmapped SAM file.
 * <p/>
 * The order of processing is as follows:
 * <p/>
 * 1.  Get records from the unmapped bam and the alignment data
 * 2.  Merge the alignment information and public tags ONLY from the aligned SAMRecords
 * 3.  Do additional modifications -- handle clipping, trimming, etc.
 * 4.  Fix up mate information on paired reads
 * 5.  Do a final calculation of the NM and UQ tags (coordinate sorted only)
 * 6.  Write the records to the output file.
 * <p/>
 * Concrete subclasses which extend AbstractAlignmentMerger should implement getQueryNameSortedAlignedRecords.
 * If these records are not in queryname order, mergeAlignment will throw an IllegalStateException.
 * <p/>
 * Subclasses may optionally implement ignoreAlignment(), which can be used to skip over certain alignments.
 *
 * @author ktibbett@broadinstitute.org
 */
public abstract class AbstractAlignmentMerger {

    public static final int MAX_RECORDS_IN_RAM = 500000;

    private static final char[] RESERVED_ATTRIBUTE_STARTS = {'X', 'Y', 'Z'};
    private int crossSpeciesReads = 0;

    private final Log log = Log.getInstance(AbstractAlignmentMerger.class);
    private final ProgressLogger progress = new ProgressLogger(this.log, 1000000, "Merged", "records");

    private final File unmappedBamFile;
    private final File targetBamFile;
    private ReferenceSequenceFileWalker refSeq = null;
    private final boolean clipAdapters;
    private final boolean bisulfiteSequence;
    private SAMProgramRecord programRecord;
    private final boolean alignedReadsOnly;
    private final SAMFileHeader header;
    private final List<String> attributesToRetain = new ArrayList<>();
    private final List<String> attributesToRemove = new ArrayList<>();
    private Set<String> attributesToReverse = new TreeSet<>(SAMRecord.TAGS_TO_REVERSE);
    private Set<String> attributesToReverseComplement = new TreeSet<>(SAMRecord.TAGS_TO_REVERSE_COMPLEMENT);
    protected final File referenceFasta;
    private final Integer read1BasesTrimmed;
    private final Integer read2BasesTrimmed;
    private final List<SamPairUtil.PairOrientation> expectedOrientations;
    private final SortOrder sortOrder;
    private MultiHitAlignedReadIterator alignedIterator = null;
    private boolean clipOverlappingReads = true;
    private int maxRecordsInRam = MAX_RECORDS_IN_RAM;
    private final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy;
    private boolean keepAlignerProperPairFlags = false;
    private boolean addMateCigar = false;
    private boolean unmapContaminantReads = false;
    private UnmappingReadStrategy unmappingReadsStrategy = UnmappingReadStrategy.DO_NOT_CHANGE;


    private final SamRecordFilter alignmentFilter = new SamRecordFilter() {
        public boolean filterOut(final SAMRecord record) {
            return ignoreAlignment(record);
        }

        public boolean filterOut(final SAMRecord first, final SAMRecord second) {
            throw new UnsupportedOperationException("Paired SamRecordFilter not implemented!");
        }
    };

    private boolean includeSecondaryAlignments = true;

    /** Class that allows a Sorting Collection and a SAMFileWriter to be treated identically. */
    private static class Sink {
        private final SAMFileWriter writer;
        private final SortingCollection<SAMRecord> sorter;

        /** Constructs a sink that outputs to a SAMFileWriter. */
        public Sink(final SAMFileWriter writer) {
            this.writer = writer;
            this.sorter = null;
        }

        /** Constructs a sink that outputs to a Sorting Collection. */
        public Sink(final SortingCollection<SAMRecord> sorter) {
            this.writer = null;
            this.sorter = sorter;
        }

        /** Adds a record to the sink. */
        void add(final SAMRecord rec) {
            if (writer != null) writer.addAlignment(rec);
            if (sorter != null) sorter.add(rec);
        }

        /** Closes the underlying resource. */
        void close() {
            if (this.writer != null) this.writer.close();
            if (this.sorter != null) this.sorter.doneAdding();
        }
    }

    public enum UnmappingReadStrategy {
        // Leave on record, and copy to tag
        COPY_TO_TAG(false, true),
        // Leave on record, but do not create additional tag
        DO_NOT_CHANGE(false, false),
        // Add tag with information, and remove from standard fields in record
        MOVE_TO_TAG(true, true);

        private final boolean resetMappingInformation, populatePATag;

        UnmappingReadStrategy(final boolean resetMappingInformation, final boolean populatePATag) {
            this.resetMappingInformation = resetMappingInformation;
            this.populatePATag = populatePATag;
        }

        public boolean isResetMappingInformation() {
            return resetMappingInformation;
        }

        public boolean isPopulatePaTag() {
            return populatePATag;
        }
    }
    protected abstract SAMSequenceDictionary getDictionaryForMergedBam();

    protected abstract CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords();

    protected boolean ignoreAlignment(final SAMRecord sam) { return false; } // default implementation

    protected boolean isContaminant(final HitsForInsert hits) { return false; } // default implementation

    /** constructor with a default setting for unmappingReadsStrategy.
     *
     * see full constructor for parameters
     *
     *
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean alignedReadsOnly,
                                   final SAMProgramRecord programRecord, final List<String> attributesToRetain,
                                   final List<String> attributesToRemove,
                                   final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                                   final List<SamPairUtil.PairOrientation> expectedOrientations,
                                   final SortOrder sortOrder,
                                   final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy,
                                   final boolean addMateCigar,
                                   final boolean unmapContaminantReads) {
        this(unmappedBamFile,
             targetBamFile,
             referenceFasta,
             clipAdapters,
             bisulfiteSequence,
             alignedReadsOnly,
             programRecord,
             attributesToRetain,
             attributesToRemove,
             read1BasesTrimmed,
             read2BasesTrimmed,
             expectedOrientations,
             sortOrder,
             primaryAlignmentSelectionStrategy,
             addMateCigar,
             unmapContaminantReads,
             UnmappingReadStrategy.DO_NOT_CHANGE);
    }


    /**
     * Constructor
     *
     * @param unmappedBamFile                   The BAM file that was used as the input to the aligner, which will
     *                                          include info on all the reads that did not map.  Required.
     * @param targetBamFile                     The file to which to write the merged SAM records. Required.
     * @param referenceFasta                    The reference sequence for the map files. Required.
     * @param clipAdapters                      Whether adapters marked in unmapped BAM file should be marked as
     *                                          soft clipped in the merged bam. Required.
     * @param bisulfiteSequence                 Whether the reads are bisulfite sequence (used when calculating the
     *                                          NM and UQ tags). Required.
     * @param alignedReadsOnly                  Whether to output only those reads that have alignment data
     * @param programRecord                     Program record for target file SAMRecords created.
     * @param attributesToRetain                private attributes from the alignment record that should be
     *                                          included when merging.  This overrides the exclusion of
     *                                          attributes whose tags start with the reserved characters
     *                                          of X, Y, and Z
     * @param attributesToRemove                attributes from the alignment record that should be
     *                                          removed when merging.  This overrides attributesToRetain if they share
     *                                          common tags.
     * @param read1BasesTrimmed                 The number of bases trimmed from start of read 1 prior to alignment.  Optional.
     * @param read2BasesTrimmed                 The number of bases trimmed from start of read 2 prior to alignment.  Optional.
     * @param expectedOrientations              A List of SamPairUtil.PairOrientations that are expected for
     *                                          aligned pairs.  Used to determine the properPair flag.
     * @param sortOrder                         The order in which the merged records should be output.  If null,
     *                                          output will be coordinate-sorted
     * @param primaryAlignmentSelectionStrategy What to do when there are multiple primary alignments, or multiple
     *                                          alignments but none primary, for a read or read pair.
     * @param addMateCigar                      True if we are to add or maintain the mate CIGAR (MC) tag, false if we are to remove or not include.
     * @param unmapContaminantReads             If true, identify reads having the signature of cross-species contamination (i.e. mostly clipped bases),
     *                                          and mark them as unmapped.
     * @param unmappingReadsStrategy            An enum describing how to deal with reads whose mapping information are being removed (currently this happens due to cross-species
     *                                          contamination). Ignored unless unmapContaminantReads is true.
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean alignedReadsOnly,
                                   final SAMProgramRecord programRecord, final List<String> attributesToRetain,
                                   final List<String> attributesToRemove,
                                   final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                                   final List<SamPairUtil.PairOrientation> expectedOrientations,
                                   final SortOrder sortOrder,
                                   final PrimaryAlignmentSelectionStrategy primaryAlignmentSelectionStrategy,
                                   final boolean addMateCigar,
                                   final boolean unmapContaminantReads,
                                   final UnmappingReadStrategy unmappingReadsStrategy) {
        IOUtil.assertFileIsReadable(unmappedBamFile);
        IOUtil.assertFileIsWritable(targetBamFile);
        IOUtil.assertFileIsReadable(referenceFasta);

        this.unmappedBamFile = unmappedBamFile;
        this.targetBamFile = targetBamFile;
        this.referenceFasta = referenceFasta;

        this.refSeq = new ReferenceSequenceFileWalker(referenceFasta);

        this.clipAdapters = clipAdapters;
        this.bisulfiteSequence = bisulfiteSequence;
        this.alignedReadsOnly = alignedReadsOnly;

        this.header = new SAMFileHeader();
        this.sortOrder = sortOrder != null ? sortOrder : SortOrder.coordinate;
        header.setSortOrder(SortOrder.coordinate);
        if (programRecord != null) {
            setProgramRecord(programRecord);
        }

        if (attributesToRetain != null) {
            this.attributesToRetain.addAll(attributesToRetain);
        }
        if (attributesToRemove != null) {
            this.attributesToRemove.addAll(attributesToRemove);
            // attributesToRemove overrides attributesToRetain
            if (!this.attributesToRetain.isEmpty()) {
                this.attributesToRemove.stream()
                        .filter(this.attributesToRetain::contains)
                        .peek(a->log.info("Overriding retaining the " + a + " tag since 'remove' overrides 'retain'."))
                        .forEach(this.attributesToRetain::remove);
            }
        }
        this.read1BasesTrimmed = read1BasesTrimmed;
        this.read2BasesTrimmed = read2BasesTrimmed;
        this.expectedOrientations = expectedOrientations;

        this.primaryAlignmentSelectionStrategy = primaryAlignmentSelectionStrategy;

        this.addMateCigar = addMateCigar;
        this.unmapContaminantReads = unmapContaminantReads;
        this.unmappingReadsStrategy = unmappingReadsStrategy;
    }

    /** Gets the set of attributes to be reversed on reads marked as negative strand. */
    public Set<String> getAttributesToReverse() { return attributesToReverse; }

    /** Sets the set of attributes to be reversed on reads marked as negative strand. */
    public void setAttributesToReverse(Set<String> attributesToReverse) { this.attributesToReverse = attributesToReverse; }

    /** Gets the set of attributes to be reverse complemented on reads marked as negative strand. */
    public Set<String> getAttributesToReverseComplement() { return attributesToReverseComplement; }

    /** Sets the set of attributes to be reverse complemented on reads marked as negative strand. */
    public void setAttributesToReverseComplement(Set<String> attributesToReverseComplement) {
        this.attributesToReverseComplement = attributesToReverseComplement;
    }

    /** Allows the caller to override the maximum records in RAM. */
    public void setMaxRecordsInRam(final int maxRecordsInRam) {
        this.maxRecordsInRam = maxRecordsInRam;
    }

    /**
     * Do this unconditionally, not just for aligned records, for two reasons:
     * - An unaligned read has been processed by the aligner, so it is more truthful.
     * - When chaining additional PG records, having all the records in the output file refer to the same PG
     * record means that only one chain will need to be created, rather than a chain for the mapped reads
     * and a separate chain for the unmapped reads.
     */
    private void maybeSetPgTag(final SAMRecord rec) {
        if (this.programRecord != null) {
            rec.setAttribute(ReservedTagConstants.PROGRAM_GROUP_ID, this.programRecord.getProgramGroupId());
        }
    }

    /**
     * Merges the alignment data with the non-aligned records from the source BAM file.
     */
    public void mergeAlignment(final File referenceFasta) {
        // Open the file of unmapped records and write the read groups to the the header for the merged file
        final SamReader unmappedSam = SamReaderFactory.makeDefault().referenceSequence(referenceFasta).open(this.unmappedBamFile);

        final CloseableIterator<SAMRecord> unmappedIterator = unmappedSam.iterator();
        this.header.setSequenceDictionary(getDictionaryForMergedBam());
        this.header.setReadGroups(unmappedSam.getFileHeader().getReadGroups());

        int aligned = 0;
        int unmapped = 0;

        // Get the aligned records and set up the first one
        alignedIterator = new MultiHitAlignedReadIterator(new FilteringSamIterator(getQuerynameSortedAlignedRecords(), alignmentFilter), primaryAlignmentSelectionStrategy);
        HitsForInsert nextAligned = nextAligned();

        // Check that the program record we are going to insert is not already used in the unmapped SAM
        // Must come after calling getQuerynameSortedAlignedRecords() in case opening the aligned records
        // sets the program group
        if (getProgramRecord() != null) {
            for (final SAMProgramRecord pg : unmappedSam.getFileHeader().getProgramRecords()) {
                if (pg.getId().equals(getProgramRecord().getId())) {
                    throw new PicardException("Program Record ID already in use in unmapped BAM file.");
                }
            }
        }

        // If the output requested is coordinate order then run everything through a sorting collection
        // in order to have access to the records in coordinate order prior to outputting them. Otherwise
        // write directly to the output BAM file in queryname order.
        final Sink sink;
        if (this.sortOrder == SortOrder.coordinate) {
            final SortingCollection<SAMRecord> sorted1 = SortingCollection.newInstance(
                    SAMRecord.class, new BAMRecordCodec(header), new SAMRecordCoordinateComparator(),
                    MAX_RECORDS_IN_RAM);
            sink = new Sink(sorted1);
        }
        else { // catches queryname and unsorted
            final SAMFileHeader header = this.header.clone();
            header.setSortOrder(this.sortOrder);
            final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, this.targetBamFile);
            writer.setProgressLogger(new ProgressLogger(log, (int) 1e7, "Wrote", "records to output in queryname order"));
            sink = new Sink(writer);
        }

        while (unmappedIterator.hasNext()) {
            // Load next unaligned read or read pair.
            final SAMRecord rec = unmappedIterator.next();

            rec.setHeader(this.header);
            maybeSetPgTag(rec);

            final SAMRecord secondOfPair;
            if (rec.getReadPairedFlag()) {
                secondOfPair = unmappedIterator.next();
                secondOfPair.setHeader(this.header);
                maybeSetPgTag(secondOfPair);

                // Validate that paired reads arrive as first of pair followed by second of pair
                if (!rec.getReadName().equals(secondOfPair.getReadName()))
                    throw new PicardException("Second read from pair not found in unmapped bam: " + rec.getReadName() + ", " + secondOfPair.getReadName());

                if (!rec.getFirstOfPairFlag())
                    throw new PicardException("First record in unmapped bam is not first of pair: " + rec.getReadName());
                if (!secondOfPair.getReadPairedFlag())
                    throw new PicardException("Second record in unmapped bam is not marked as paired: " + secondOfPair.getReadName());
                if (!secondOfPair.getSecondOfPairFlag())
                    throw new PicardException("Second record in unmapped bam is not second of pair: " + secondOfPair.getReadName());
            } else {
                secondOfPair = null;
            }

            // See if there are alignments for current unaligned read or read pair.
            if (nextAligned != null && rec.getReadName().equals(nextAligned.getReadName())) {
                // If there are multiple alignments for a read (pair), then the unaligned SAMRecord must be cloned
                // before copying info from the aligned record to the unaligned.
                final boolean clone = nextAligned.numHits() > 1 || nextAligned.hasSupplementalHits();
                SAMRecord r1Primary = null, r2Primary = null;

                // by this point there should be a single chosen primary alignment, which we will use to determine whether the read is contaminant.
                // this must be done before the main iteration, since secondary / supplementary alignments will be affected by the primary.
                final boolean unmapDueToContaminant = this.unmapContaminantReads && isContaminant(nextAligned);

                if (rec.getReadPairedFlag()) {
                    for (int i = 0; i < nextAligned.numHits(); ++i) {
                        // firstAligned or secondAligned may be null, if there wasn't an alignment for the end,
                        // or if the alignment was rejected by ignoreAlignment.
                        final SAMRecord firstAligned = nextAligned.getFirstOfPair(i);
                        final SAMRecord secondAligned = nextAligned.getSecondOfPair(i);

                        final boolean isPrimaryAlignment = (firstAligned != null && !firstAligned.isSecondaryOrSupplementary()) ||
                                (secondAligned != null && !secondAligned.isSecondaryOrSupplementary());

                        final SAMRecord firstToWrite;
                        final SAMRecord secondToWrite;
                        if (clone) {
                            firstToWrite = clone(rec);
                            secondToWrite = clone(secondOfPair);
                        } else {
                            firstToWrite = rec;
                            secondToWrite = secondOfPair;
                        }

                        // If these are the primary alignments then stash them for use on any supplemental alignments
                        if (isPrimaryAlignment) {
                            r1Primary = firstToWrite;
                            r2Primary = secondToWrite;
                        }

                        transferAlignmentInfoToPairedRead(firstToWrite, secondToWrite, firstAligned, secondAligned, unmapDueToContaminant, clone);

                        // Only write unmapped read when it has the mate info from the primary alignment.
                        // this avoids the scenario of having multiple unmapped reads with the same name & pair flags
                        if (!firstToWrite.getReadUnmappedFlag() || isPrimaryAlignment) {
                            addIfNotFiltered(sink, firstToWrite);
                            if (firstToWrite.getReadUnmappedFlag()) ++unmapped;
                            else ++aligned;
                        }
                        if (!secondToWrite.getReadUnmappedFlag() || isPrimaryAlignment) {
                            addIfNotFiltered(sink, secondToWrite);
                            if (!secondToWrite.getReadUnmappedFlag()) ++aligned;
                            else ++unmapped;
                        }
                    }

                    // Take all of the supplemental reads which had been stashed and add them (as appropriate) to sorted
                    for (final boolean isRead1 : new boolean[]{true, false}) {
                        final List<SAMRecord> supplementals = isRead1 ? nextAligned.getSupplementalFirstOfPairOrFragment() : nextAligned.getSupplementalSecondOfPair();
                        final SAMRecord sourceRec = isRead1 ? rec : secondOfPair;
                        final SAMRecord matePrimary = isRead1 ? r2Primary : r1Primary;

                        for (final SAMRecord supp : supplementals) {
                            final SAMRecord out = clone(sourceRec);
                            transferAlignmentInfoToFragment(out, supp, unmapDueToContaminant, clone);
                            if (matePrimary != null) SamPairUtil.setMateInformationOnSupplementalAlignment(out, matePrimary, addMateCigar);
                            // don't write supplementary reads that were unmapped by transferAlignmentInfoToFragment
                            if (!out.getReadUnmappedFlag()) {
                                addIfNotFiltered(sink, out);
                                ++aligned;
                            } else ++unmapped;
                        }
                    }
                } else {
                    for (int i = 0; i < nextAligned.numHits(); ++i) {
                        final SAMRecord recToWrite = clone ? clone(rec) : rec;
                        final boolean isPrimary = !nextAligned.getFragment(i).isSecondaryOrSupplementary();
                        transferAlignmentInfoToFragment(recToWrite, nextAligned.getFragment(i), unmapDueToContaminant, clone);
                        // Only write unmapped read if it was originally the primary.
                        // this avoids the scenario of having multiple unmapped reads with the same name & pair flags
                        if (!recToWrite.getReadUnmappedFlag() || isPrimary) addIfNotFiltered(sink, recToWrite);
                        if (recToWrite.getReadUnmappedFlag()) ++unmapped;
                        else ++aligned;
                    }
                    // Take all of the supplemental reads which had been stashed and add them (as appropriate) to sorted
                    for (final SAMRecord supplementalRec : nextAligned.getSupplementalFirstOfPairOrFragment()) {
                        final SAMRecord recToWrite = clone(rec);
                        transferAlignmentInfoToFragment(recToWrite, supplementalRec, unmapDueToContaminant, clone);
                        // don't write supplementary reads that were unmapped by transferAlignmentInfoToFragment
                        if (!recToWrite.getReadUnmappedFlag()) {
                            addIfNotFiltered(sink, recToWrite);
                            ++aligned;
                        } else ++unmapped;
                    }
                }
                nextAligned = nextAligned();
            } else {
                // There was no alignment for this read or read pair.
                if (nextAligned != null &&
                        SAMRecordQueryNameComparator.compareReadNames(rec.getReadName(), nextAligned.getReadName()) > 0) {
                    throw new IllegalStateException("Aligned record iterator (" + nextAligned.getReadName() +
                            ") is behind the unmapped reads (" + rec.getReadName() + ")");
                }
                // No matching read from alignedIterator -- just output reads as is.
                if (!alignedReadsOnly) {
                    sink.add(rec);
                    ++unmapped;
                    if (secondOfPair != null) {
                        sink.add(secondOfPair);
                        ++unmapped;
                    }
                }
            }
        }
        unmappedIterator.close();
        if (alignedIterator.hasNext()) {
            throw new IllegalStateException("Reads remaining on alignment iterator: " + alignedIterator.next().getReadName() + "!");
        }
        alignedIterator.close();
        sink.close();

        // Write the records to the output file in specified sorted order,
        if (this.sortOrder == SortOrder.coordinate) {
            header.setSortOrder(this.sortOrder);
            final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, this.targetBamFile);
            writer.setProgressLogger(new ProgressLogger(log, (int) 1e7, "Wrote", "records from a sorting collection"));
            final ProgressLogger finalProgress = new ProgressLogger(log, 10000000, "Written in coordinate order to output", "records");

            for (final SAMRecord rec : sink.sorter) {
                if (!rec.getReadUnmappedFlag() && refSeq != null) {
                    fixNmMdAndUq(rec, refSeq, bisulfiteSequence);
                }
                writer.addAlignment(rec);
                finalProgress.record(rec);
            }
            writer.close();
            sink.sorter.cleanup();
        }

        CloserUtil.close(unmappedSam);
        log.info("Wrote " + aligned + " alignment records and " + (alignedReadsOnly ? 0 : unmapped) + " unmapped reads.");
    }

    /** Calculates and sets the NM, MD, and and UQ tags from the record and the reference
     *
     * @param record the record to be fixed
     * @param refSeqWalker a ReferenceSequenceWalker that will be used to traverse the reference
     * @param isBisulfiteSequence a flag indicating whether the sequence came from bisulfite-sequencing which would imply a different
     * calculation of the NM tag.
     *
     * No return value, modifies the provided record.
     */
    public static void fixNmMdAndUq(final SAMRecord record, final ReferenceSequenceFileWalker refSeqWalker, final boolean isBisulfiteSequence) {
        final byte[] referenceBases = refSeqWalker.get(refSeqWalker.getSequenceDictionary().getSequenceIndex(record.getReferenceName())).getBases();
        // only recalculate NM if it isn't bisulfite, since it needs to be treated specially below
        SequenceUtil.calculateMdAndNmTags(record, referenceBases, true, !isBisulfiteSequence);
        if (isBisulfiteSequence) {  // recalculate the NM tag for bisulfite data
            record.setAttribute(SAMTag.NM.name(), SequenceUtil.calculateSamNmTag(record, referenceBases, 0, isBisulfiteSequence));
        }
        if (record.getBaseQualities() != SAMRecord.NULL_QUALS) {
            record.setAttribute(SAMTag.UQ.name(), SequenceUtil.sumQualitiesOfMismatches(record, referenceBases, 0, isBisulfiteSequence));
        }
    }

    /**
     * Add record if it is primary or optionally secondary.
     */
    private void addIfNotFiltered(final Sink out, final SAMRecord rec) {
        if (includeSecondaryAlignments || !rec.getNotPrimaryAlignmentFlag()) {
            out.add(rec);
            if (this.progress.record(rec) && crossSpeciesReads > 0) {
                log.info(String.format("%d Reads have been unmapped due to being suspected of being Cross-species contamination.", crossSpeciesReads));
            }
        }
    }

    private SAMRecord clone(final SAMRecord rec) {
        try {
            return (SAMRecord) rec.clone();
        } catch (CloneNotSupportedException e) {
            throw new PicardException("Should never happen.");
        }
    }

    /**
     * @return Next read's alignment(s) from aligned input or null, if there are no more.
     * The alignments are run through ignoreAlignment() filter before being returned, which may result
     * in an entire read being skipped if all alignments for that read should be ignored.
     */
    private HitsForInsert nextAligned() {
        if (alignedIterator.hasNext()) return alignedIterator.next();
        return null;
    }

    /**
     * Copies alignment info from aligned to unaligned read, clips as appropriate, and sets PG ID.
     * May also un-map the resulting read if the alignment is bad (e.g. no unclipped bases).
     *
     * Depending on the value of unmappingReadsStrategy this function will potentially
     *
     * - reset the mapping information (MOVE_TO_TAG)
     * - copy the original mapping information to a tag "PA" (MOVE_TO_TAG, or COPY_TO_TAG)
     * a third possibility for unmappingReadsStrategy is to do nothig (DO_NOT_CHANGE)
     *
     * This is because the CRAM format will reset mapping information for reads that are flagged as unaligned.
     *
     * @param unaligned Original SAMRecord, and object into which values are copied.
     * @param aligned   Holds alignment info that will be copied into unaligned.
     * @param isContaminant Should this read be unmapped due to contamination?
     */
    private void transferAlignmentInfoToFragment(final SAMRecord unaligned, final SAMRecord aligned, final boolean isContaminant, final boolean needsSafeReverseComplement) {
        setValuesFromAlignment(unaligned, aligned, needsSafeReverseComplement);
        updateCigarForTrimmedOrClippedBases(unaligned, aligned);
        if (SAMUtils.cigarMapsNoBasesToRef(unaligned.getCigar())) {
            log.warn("Record contains no unclipped bases; making unmapped: " + aligned);
            SAMUtils.makeReadUnmapped(unaligned);
        } else if (SAMUtils.recordMapsEntirelyBeyondEndOfReference(aligned)) {
            log.warn("Record mapped off end of reference; making unmapped: " + aligned);
            SAMUtils.makeReadUnmapped(unaligned);
        } else if (isContaminant) {

            crossSpeciesReads++;

            if (unmappingReadsStrategy.isPopulatePaTag()) {
                unaligned.setAttribute("PA", encodeMappingInformation(aligned));
            }

            if (unmappingReadsStrategy.isResetMappingInformation()) {
                unaligned.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                unaligned.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                unaligned.setCigar(null);
                unaligned.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            }

            unaligned.setReadUnmappedFlag(true);
            // Unmapped read cannot have non-zero mapping quality and remain valid
            unaligned.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);

            // if there already is a comment, add second comment with a | separator:
            Optional<String> optionalComment = Optional.ofNullable(unaligned.getStringAttribute(SAMTag.CO.name()));
            unaligned.setAttribute(optionalComment.map(s -> s + " | ").orElse("") +
                    SAMTag.CO.name(), "Cross-species contamination");
        }
    }

    /**
     * Encodes mapping information from a record into a string according to the format sepcified in the
     * Sam-Spec under the SA tag. No protection against missing values (for cigar, and NM tag).
     * (Might make sense to move this to htsJDK.)
     *
     * @param rec SAMRecord whose alignment information will be encoded
     * @return String encoding rec's alignment information according to SA tag in the SAM spec
     */
    static private String encodeMappingInformation(SAMRecord rec) {
        return String.join(",",
                rec.getContig(),
                ((Integer) rec.getAlignmentStart()).toString(),
                rec.getCigarString(),
                ((Integer) rec.getMappingQuality()).toString(),
                rec.getStringAttribute(SAMTag.NM.name()))+";";
    }

    /**
     * Copies alignment info from aligned to unaligned read, if there is an alignment, and sets mate information.
     *
     * @param firstUnaligned  Original first of pair, into which alignment and pair info will be written.
     * @param secondUnaligned Original second of pair, into which alignment and pair info will be written.
     * @param firstAligned    Aligned first of pair, or null if no alignment.
     * @param secondAligned   Aligned second of pair, or null if no alignment.
     * @param isContaminant Should this pair be unmapped due to contamination?
     */
    private void transferAlignmentInfoToPairedRead(final SAMRecord firstUnaligned, final SAMRecord secondUnaligned,
                                                   final SAMRecord firstAligned, final SAMRecord secondAligned,
                                                   final boolean isContaminant, final boolean needsSafeReverseComplement) {
        if (firstAligned != null) transferAlignmentInfoToFragment(firstUnaligned, firstAligned, isContaminant, needsSafeReverseComplement);
        if (secondAligned != null) transferAlignmentInfoToFragment(secondUnaligned, secondAligned, isContaminant, needsSafeReverseComplement);
        if (isClipOverlappingReads()) clipForOverlappingReads(firstUnaligned, secondUnaligned);
        SamPairUtil.setMateInfo(secondUnaligned, firstUnaligned, addMateCigar);
        if (!keepAlignerProperPairFlags) {
            SamPairUtil.setProperPairFlags(secondUnaligned, firstUnaligned, expectedOrientations);
        }
    }

    /**
     * Checks to see whether the ends of the reads overlap and soft clips reads
     * them if necessary.
     */
    protected static void clipForOverlappingReads(final SAMRecord read1, final SAMRecord read2) {
        // If both reads are mapped, see if we need to clip the ends due to small
        // insert size
        if (!(read1.getReadUnmappedFlag() || read2.getReadUnmappedFlag())) {
            if (read1.getReadNegativeStrandFlag() != read2.getReadNegativeStrandFlag()) {
                final SAMRecord pos = (read1.getReadNegativeStrandFlag()) ? read2 : read1;
                final SAMRecord neg = (read1.getReadNegativeStrandFlag()) ? read1 : read2;

                // Innies only -- do we need to do anything else about jumping libraries?
                if (pos.getAlignmentStart() < neg.getAlignmentEnd()) {
                    final int posDiff = pos.getAlignmentEnd() - neg.getAlignmentEnd();
                    final int negDiff = pos.getAlignmentStart() - neg.getAlignmentStart();

                    if (posDiff > 0) {
                        final List<CigarElement> elems = new ArrayList<>(pos.getCigar().getCigarElements());
                        Collections.reverse(elems);
                        final int clipped = lengthOfSoftClipping(elems.iterator());
                        final int clipFrom = pos.getReadLength() - posDiff - clipped + 1;
                        CigarUtil.softClip3PrimeEndOfRead(pos, Math.min(pos.getReadLength(), clipFrom));
                        removeNmMdAndUqTags(pos); // these tags are now invalid!
                    }

                    if (negDiff > 0) {
                        final int clipped = lengthOfSoftClipping(neg.getCigar().getCigarElements().iterator());
                        final int clipFrom = neg.getReadLength() - negDiff - clipped + 1;
                        CigarUtil.softClip3PrimeEndOfRead(neg, Math.min(neg.getReadLength(), clipFrom));
                        removeNmMdAndUqTags(neg); // these tags are now invalid!
                    }
                }
            }
        }
    }

    /** Returns the number of soft-clipped bases until a non-soft-clipping element is encountered. */
    private static int lengthOfSoftClipping(Iterator<CigarElement> iterator) {
        int clipped = 0;
        while (iterator.hasNext()) {
            final CigarElement elem = iterator.next();
            if (elem.getOperator() != CigarOperator.SOFT_CLIP && elem.getOperator() != CigarOperator.HARD_CLIP) break;
            if (elem.getOperator() == CigarOperator.SOFT_CLIP) clipped = elem.getLength();
        }

        return clipped;
    }

    /**
     * Sets the values from the alignment record on the unaligned BAM record.  This
     * preserves all data from the unaligned record (ReadGroup, NoiseRead status, etc)
     * and adds all the alignment info
     *
     * @param rec       The unaligned read record
     * @param alignment The alignment record
     */
    protected void setValuesFromAlignment(final SAMRecord rec, final SAMRecord alignment, final boolean needsSafeReverseComplement) {
        for (final SAMRecord.SAMTagAndValue attr : alignment.getAttributes()) {
            // Copy over any non-reserved attributes.  attributesToRemove overrides attributesToRetain.
            if ((!isReservedTag(attr.tag) || this.attributesToRetain.contains(attr.tag)) && !this.attributesToRemove.contains(attr.tag)) {
                rec.setAttribute(attr.tag, attr.value);
            }
        }
        rec.setReadUnmappedFlag(alignment.getReadUnmappedFlag());

        // Note that it is important to get reference names rather than indices in case the sequence dictionaries
        // in the two files are in different orders.
        rec.setReferenceName(alignment.getReferenceName());

        rec.setAlignmentStart(alignment.getAlignmentStart());
        rec.setReadNegativeStrandFlag(alignment.getReadNegativeStrandFlag());
        rec.setNotPrimaryAlignmentFlag(alignment.getNotPrimaryAlignmentFlag());
        rec.setSupplementaryAlignmentFlag(alignment.getSupplementaryAlignmentFlag());
        if (!alignment.getReadUnmappedFlag()) {
            // only aligned reads should have cigar and mapping quality set
            rec.setCigar(alignment.getCigar());  // cigar may change when a
            // clipCigar called below
            rec.setMappingQuality(alignment.getMappingQuality());
        }
        if (rec.getReadPairedFlag()) {
            rec.setProperPairFlag(alignment.getProperPairFlag());
            // Mate info and alignment size will get set by the ClippedPairFixer.
        }

        // If it's on the negative strand, reverse complement the bases
        // and reverse the order of the qualities
        if (rec.getReadNegativeStrandFlag()) {
            rec.reverseComplement(attributesToReverseComplement, attributesToReverse, !needsSafeReverseComplement);
        }
    }

    private static Cigar createNewCigarIfMapsOffEndOfReference(SAMFileHeader header,
                                                               boolean isUnmapped,
                                                               int referenceIndex,
                                                               int alignmentEnd,
                                                               int readLength,
                                                               Cigar oldCigar) {
        Cigar newCigar = null;
        if (!isUnmapped) {
            final SAMSequenceRecord refseq = header.getSequence(referenceIndex);
            final int overhang = alignmentEnd - refseq.getSequenceLength();
            if (overhang > 0) {
                // 1-based index of first base in read to clip.
                int clipFrom = readLength - overhang + 1;
                // we have to check if the last element is soft-clipping, so we can subtract that from clipFrom
                final CigarElement cigarElement = oldCigar.getCigarElement(oldCigar.getCigarElements().size()-1);
                if (CigarOperator.SOFT_CLIP == cigarElement.getOperator()) clipFrom -= cigarElement.getLength();
                final List<CigarElement> newCigarElements = CigarUtil.softClipEndOfRead(clipFrom, oldCigar.getCigarElements());
                newCigar = new Cigar(newCigarElements);
            }
        }
        return newCigar;
    }

    /**
     * Soft-clip an alignment that hangs off the end of its reference sequence.  Checks both the read and its mate,
     * if available.
     *
     * @param rec
     */
    public static void createNewCigarsIfMapsOffEndOfReference(final SAMRecord rec) {
        // If the read maps off the end of the alignment, clip it
        if (!rec.getReadUnmappedFlag()) {
            final Cigar readCigar = createNewCigarIfMapsOffEndOfReference(rec.getHeader(),
                    rec.getReadUnmappedFlag(),
                    rec.getReferenceIndex(),
                    rec.getAlignmentEnd(),
                    rec.getReadLength(),
                    rec.getCigar());
            if (null != readCigar) rec.setCigar(readCigar);
        }

        // If the read's mate maps off the end of the alignment, clip it
        if (SAMUtils.hasMateCigar(rec)) {
            Cigar mateCigar = SAMUtils.getMateCigar(rec);
            mateCigar = createNewCigarIfMapsOffEndOfReference(rec.getHeader(),
                    rec.getMateUnmappedFlag(),
                    rec.getMateReferenceIndex(),
                    SAMUtils.getMateAlignmentEnd(rec), // NB: this could be computed without another call to getMateCigar
                    mateCigar.getReadLength(),
                    mateCigar);
            if (null != mateCigar) rec.setAttribute(SAMTag.MC.name(), mateCigar.toString());
        }
    }

    protected void updateCigarForTrimmedOrClippedBases(final SAMRecord rec, final SAMRecord alignment) {
        // If the read was trimmed or not all the bases were sent for alignment, clip it
        final int alignmentReadLength = alignment.getReadLength();
        final int originalReadLength = rec.getReadLength();
        final int trimmed = (!rec.getReadPairedFlag()) || rec.getFirstOfPairFlag()
                ? this.read1BasesTrimmed != null ? this.read1BasesTrimmed : 0
                : this.read2BasesTrimmed != null ? this.read2BasesTrimmed : 0;
        final int notWritten = originalReadLength - (alignmentReadLength + trimmed);

        // Update cigar if the mate maps off the reference
        createNewCigarsIfMapsOffEndOfReference(rec);

        rec.setCigar(CigarUtil.addSoftClippedBasesToEndsOfCigar(
                rec.getCigar(), rec.getReadNegativeStrandFlag(), notWritten, trimmed));

        // If the adapter sequence is marked and clipAdapter is true, clip it
        if (this.clipAdapters && rec.getAttribute(ReservedTagConstants.XT) != null) {
            CigarUtil.softClip3PrimeEndOfRead(rec, rec.getIntegerAttribute(ReservedTagConstants.XT));
            removeNmMdAndUqTags(rec); // these tags are now invalid!
        }
    }

    protected SAMProgramRecord getProgramRecord() { return this.programRecord; }

    protected void setProgramRecord(final SAMProgramRecord pg) {
        if (this.programRecord != null) {
            throw new IllegalStateException("Cannot set program record more than once on alignment merger.");
        }
        this.programRecord = pg;
        this.header.addProgramRecord(pg);
        SAMUtils.chainSAMProgramRecord(header, pg);
    }

    protected boolean isReservedTag(final String tag) {
        final char firstCharOfTag = tag.charAt(0);

        // All tags that start with a lower-case letter are user defined and should not be overridden by aligner
        // unless explicitly specified in attributesToRetain.
        if (Character.isLowerCase(firstCharOfTag)) return true;

        for (final char c : RESERVED_ATTRIBUTE_STARTS) {
            if (firstCharOfTag == c) return true;
        }
        return false;
    }

    protected SAMFileHeader getHeader() { return this.header; }

    protected void resetRefSeqFileWalker() {
        this.refSeq = new ReferenceSequenceFileWalker(referenceFasta);
    }

    public boolean isClipOverlappingReads() {
        return clipOverlappingReads;
    }

    public void setClipOverlappingReads(final boolean clipOverlappingReads) {
        this.clipOverlappingReads = clipOverlappingReads;
    }

    public boolean isKeepAlignerProperPairFlags() {
        return keepAlignerProperPairFlags;
    }

    /**
     * If true, keep the aligner's idea of proper pairs rather than letting alignment merger decide.
     */
    public void setKeepAlignerProperPairFlags(final boolean keepAlignerProperPairFlags) {
        this.keepAlignerProperPairFlags = keepAlignerProperPairFlags;
    }

    public void setIncludeSecondaryAlignments(final boolean includeSecondaryAlignments) {
        this.includeSecondaryAlignments = includeSecondaryAlignments;
    }

    public void close() {
        CloserUtil.close(this.refSeq);
    }


    /** Removes the NM, MD, and UQ tags.  This is useful if we modify the read and are not able to recompute these tags,
     * for example when no reference is available.
     * @param rec the record to modify.
     */
    private static void removeNmMdAndUqTags(final SAMRecord rec) {
        rec.setAttribute(SAMTag.NM.name(), null);
        rec.setAttribute(SAMTag.MD.name(), null);
        rec.setAttribute(SAMTag.UQ.name(), null);
    }
}
