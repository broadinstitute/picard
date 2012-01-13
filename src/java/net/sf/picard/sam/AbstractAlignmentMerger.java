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

import net.sf.picard.PicardException;
import net.sf.picard.filter.SamRecordFilter;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.CigarUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.SortingCollection;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Abstract class that coordinates the general task of taking in a set of alignment information,
 * possibly in SAM format, possibly in other formats, and merging that with the set of all reads
 * for which alignment was attempted, stored in an unmapped SAM file.
 *
 * The order of processing is as follows:
 *
 *   1.  Get records from the unmapped bam and the alignment data
 *   2.  Merge the alignment information and public tags ONLY from the aligned SAMRecords
 *   3.  Do additional modifications -- handle clipping, trimming, etc.
 *   4.  Fix up mate information on paired reads
 *   5.  Do a final calculation of the NM and UQ tags.
 *   6.  Write the records to the output file.
 *
 * Concrete subclasses which extend AbstractAlignmentMerger should implement getQueryNameSortedAlignedRecords.
 * If these records are not in queryname order, mergeAlignment will throw an IllegalStateException.
 *
 * Subclasses may optionally implement ignoreAlignment(), which can be used to skip over certain alignments.
 *
 *
 * @author ktibbett@broadinstitute.org
 */
public abstract class AbstractAlignmentMerger {

    public static final int MAX_RECORDS_IN_RAM = 500000;

    private static final char[] RESERVED_ATTRIBUTE_STARTS = {'X','Y', 'Z'};
    private final NumberFormat FMT = new DecimalFormat("#,###");

    private final Log log = Log.getInstance(AbstractAlignmentMerger.class);
    private final File unmappedBamFile;
    private final File targetBamFile;
    private SAMSequenceDictionary sequenceDictionary = null;
    private ReferenceSequenceFileWalker refSeq = null;
    private final boolean clipAdapters;
    private final boolean bisulfiteSequence;
    private SAMProgramRecord programRecord;
    private final boolean alignedReadsOnly;
    private final SAMFileHeader header;
    private final List<String> attributesToRetain = new ArrayList<String>();
    private final File referenceFasta;
    private final Integer read1BasesTrimmed;
    private final Integer read2BasesTrimmed;
    private final List<SamPairUtil.PairOrientation> expectedOrientations;
    private final SortOrder sortOrder;
    private MultiHitAlignedReadIterator alignedIterator = null;
    private boolean clipOverlappingReads = true;
    private int maxRecordsInRam = MAX_RECORDS_IN_RAM;

    private SamRecordFilter alignmentFilter = new SamRecordFilter() {
        public boolean filterOut(SAMRecord record) {
            return ignoreAlignment(record);
        }
        public boolean filterOut(final SAMRecord first, final SAMRecord second) {
            throw new UnsupportedOperationException("Paired SamRecordFilter not implemented!");
        }
    };

    protected abstract CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords();
    
    protected boolean ignoreAlignment(SAMRecord sam) { return false; } // default implementation

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the aligner, which will
     *                          include info on all the reads that did not map.  Required.
     * @param targetBamFile     The file to which to write the merged SAM records. Required.
     * @param referenceFasta    The reference sequence for the map files. Required.
     * @param clipAdapters      Whether adapters marked in unmapped BAM file should be marked as
     *                          soft clipped in the merged bam. Required.
     * @param bisulfiteSequence Whether the reads are bisulfite sequence (used when calculating the
     *                          NM and UQ tags). Required.
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param programRecord     Program record for taget file SAMRecords created.
     * @param attributesToRetain  private attributes from the alignment record that should be
     *                          included when merging.  This overrides the exclusion of
     *                          attributes whose tags start with the reserved characters
     *                          of X, Y, and Z
     * @param read1BasesTrimmed The number of bases trimmed from start of read 1 prior to alignment.  Optional.
     * @param read2BasesTrimmed The number of bases trimmed from start of read 2 prior to alignment.  Optional.
     * @param expectedOrientations A List of SamPairUtil.PairOrientations that are expected for
     *                          aligned pairs.  Used to determine the properPair flag.
     * @param sortOrder           The order in which the merged records should be output.  If null,
     *                            output will be coordinate-sorted
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean alignedReadsOnly,
                                   final SAMProgramRecord programRecord, final List<String> attributesToRetain,
                                   final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                                   final List<SamPairUtil.PairOrientation> expectedOrientations,
                                   final SAMFileHeader.SortOrder sortOrder) {
        IoUtil.assertFileIsReadable(unmappedBamFile);
        IoUtil.assertFileIsWritable(targetBamFile);
        IoUtil.assertFileIsReadable(referenceFasta);

        this.unmappedBamFile = unmappedBamFile;
        this.targetBamFile = targetBamFile;
        this.referenceFasta = referenceFasta;

        this.refSeq = new ReferenceSequenceFileWalker(referenceFasta);
        this.sequenceDictionary = refSeq.getSequenceDictionary();

        this.clipAdapters = clipAdapters;
        this.bisulfiteSequence = bisulfiteSequence;
        this.alignedReadsOnly = alignedReadsOnly;

        this.header = new SAMFileHeader();
        this.sortOrder = sortOrder != null ? sortOrder : SortOrder.coordinate;
        header.setSortOrder(SortOrder.coordinate);
        if (programRecord != null) {
            setProgramRecord(programRecord);
        }
        header.setSequenceDictionary(this.sequenceDictionary);
        if (attributesToRetain != null) {
            this.attributesToRetain.addAll(attributesToRetain);
        }
        this.read1BasesTrimmed = read1BasesTrimmed;
        this.read2BasesTrimmed = read2BasesTrimmed;
        this.expectedOrientations = expectedOrientations;
    }

    /**
     * Constructor retained for backwards compatibility
     *
     * @deprecated  Use construct that specifies sortOrder
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean alignedReadsOnly,
                                   final SAMProgramRecord programRecord, final List<String> attributesToRetain,
                                   final Integer read1BasesTrimmed, final Integer read2BasesTrimmed,
                                   final List<SamPairUtil.PairOrientation> expectedOrientations) {

        this(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
             alignedReadsOnly, programRecord, attributesToRetain, read1BasesTrimmed, read2BasesTrimmed,
             expectedOrientations, SortOrder.coordinate);
    }

    /** Allows the caller to override the maximum records in RAM. */
    public void setMaxRecordsInRam(final int maxRecordsInRam) {
        this.maxRecordsInRam = maxRecordsInRam;
    }

    /**
     * Merges the alignment data with the non-aligned records from the source BAM file.
     */
    public void mergeAlignment() {
        // Open the file of unmapped records and write the read groups to the the header for the merged file
        final SAMFileReader unmappedSam = new SAMFileReader(this.unmappedBamFile);
        final CloseableIterator<SAMRecord> unmappedIterator = unmappedSam.iterator();
        this.header.setReadGroups(unmappedSam.getFileHeader().getReadGroups());

        int aligned = 0;
        int unmapped = 0;

        // Get the aligned records and set up the first one
        alignedIterator = new MultiHitAlignedReadIterator(getQuerynameSortedAlignedRecords());
        MultiHitAlignedReadIterator.HitsForInsert nextAligned = nextAligned();

        // Create the sorting collection that will write the records in the coordinate order
        // to the final bam file
        final SortingCollection<SAMRecord> sorted = SortingCollection.newInstance(
            SAMRecord.class, new BAMRecordCodec(header), new SAMRecordCoordinateComparator(),
            MAX_RECORDS_IN_RAM);

        while (unmappedIterator.hasNext()) {
            // Load next unaligned read or read pair.
            final SAMRecord rec = unmappedIterator.next();
            rec.setHeader(this.header);

            final SAMRecord secondOfPair;
            if (rec.getReadPairedFlag()) {
                secondOfPair = unmappedIterator.next();
                secondOfPair.setHeader(this.header);

                // Validate that paired reads arrive as first of pair followed by second of pair
                if (!rec.getReadName().equals(secondOfPair.getReadName()))
                    throw new PicardException("Second read from pair not found in unmapped bam: " + rec.getReadName() + ", " + secondOfPair.getReadName());

                if (!rec.getFirstOfPairFlag()) throw new PicardException("First record in unmapped bam is not first of pair: " + rec.getReadName());
                if (!secondOfPair.getReadPairedFlag())  throw new PicardException("Second record in unmapped bam is not marked as paired: " + secondOfPair.getReadName());
                if (!secondOfPair.getSecondOfPairFlag())  throw new PicardException("Second record in unmapped bam is not second of pair: " + secondOfPair.getReadName());
            }
            else {
                secondOfPair = null;
            }

            // See if there are alignments for current unaligned read or read pair.
            if (nextAligned != null && rec.getReadName().equals(nextAligned.getReadName())) {
                // If there are multiple alignments for a read (pair), then the unaligned SAMRecord must be cloned
                // before copying info from the aligned record to the unaligned.
                final boolean clone = nextAligned.numHits() > 1;

                if (rec.getReadPairedFlag()) {
                    for (int i = 0; i < nextAligned.numHits(); ++i) {
                        final SAMRecord firstToWrite;
                        final SAMRecord secondToWrite;
                        if (clone) {
                            firstToWrite = clone(rec);
                            secondToWrite = clone(secondOfPair);
                        } else {
                            firstToWrite = rec;
                            secondToWrite = secondOfPair;
                        }
                        // firstAligned or secondAligned may be null, if there wasn't an alignment for the end,
                        // or if the alignment was rejected by ignoreAlignment.
                        final SAMRecord firstAligned = nextAligned.getFirstOfPair(i);
                        final SAMRecord secondAligned = nextAligned.getSecondOfPair(i);

                        final boolean isPrimaryAlignment = (firstAligned != null && !firstAligned.getNotPrimaryAlignmentFlag()) ||
                                (secondAligned != null && !secondAligned.getNotPrimaryAlignmentFlag());

                        transferAlignmentInfoToPairedRead(firstToWrite, secondToWrite, firstAligned, secondAligned);

                        // Only write unmapped read when it has the mate info from the primary alignment.
                        if (!firstToWrite.getReadUnmappedFlag() || isPrimaryAlignment) {
                            sorted.add(firstToWrite);
                            if (firstToWrite.getReadUnmappedFlag()) ++unmapped;
                            else ++aligned;
                        }
                        if (!secondToWrite.getReadUnmappedFlag() || isPrimaryAlignment) {
                            sorted.add(secondToWrite);
                            if (!secondToWrite.getReadUnmappedFlag()) ++aligned;
                            else ++unmapped;
                        }
                    }
                } else {
                    for (int i = 0; i < nextAligned.numHits(); ++i) {
                        final SAMRecord recToWrite = clone ? clone(rec) : rec;
                        transferAlignmentInfoToFragment(recToWrite, nextAligned.getFragment(i));
                        sorted.add(recToWrite);
                        if (recToWrite.getReadUnmappedFlag()) ++unmapped;
                        else ++aligned;
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
                    sorted.add(rec);
                    ++unmapped;
                    if (secondOfPair != null) {
                        sorted.add(secondOfPair);
                        ++unmapped;
                    }
                }
            }

            if ((aligned + unmapped) % 1000000 == 0) {
                log.info("Processed " + FMT.format(aligned + unmapped) + " records in query name order.");
            }
        }
        unmappedIterator.close();
        if (alignedIterator.hasNext()) {
            throw new IllegalStateException("Reads remaining on alignment iterator: " +
                    alignedIterator.next().getReadName() + "!");
        }
        alignedIterator.close();

        // Write the records to the output file in specified sorted order,
        header.setSortOrder(this.sortOrder);
        boolean presorted = this.sortOrder == SortOrder.coordinate;
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, presorted, this.targetBamFile);
        int count = 0;
        for (final SAMRecord rec : sorted) {
            if (!rec.getReadUnmappedFlag()) {
                if (refSeq != null) {
                    byte referenceBases[] = refSeq.get(sequenceDictionary.getSequenceIndex(rec.getReferenceName())).getBases();
                    rec.setAttribute(SAMTag.NM.name(),
                        SequenceUtil.calculateSamNmTag(rec, referenceBases, 0, bisulfiteSequence));
                    if (rec.getBaseQualities() != SAMRecord.NULL_QUALS) {
                        rec.setAttribute(SAMTag.UQ.name(),
                            SequenceUtil.sumQualitiesOfMismatches(rec, referenceBases, 0, bisulfiteSequence));
                    }
                }
            }
            writer.addAlignment(rec);
            if (++count % 1000000 == 0) {
                log.info(FMT.format(count) + " SAMRecords written to " + targetBamFile.getName());
            }
        }
        writer.close();
        sorted.cleanup();

        log.info("Wrote " + aligned + " alignment records and " + (alignedReadsOnly ? 0 : unmapped) + " unmapped reads.");
    }

    private SAMRecord clone(final SAMRecord rec) {
        try {
            return (SAMRecord)rec.clone();
        } catch (CloneNotSupportedException e) {
            throw new PicardException("Should never happen.");
        }
    }
    /**
     * @return Next read's alignment(s) from aligned input or null, if there are no more.
     * The alignments are run through ignoreAlignment() filter before being returned, which may result
     * in an entire read being skipped if all alignments for that read should be ignored.
     */
    private MultiHitAlignedReadIterator.HitsForInsert nextAligned() {
        while (alignedIterator.hasNext()) {
            MultiHitAlignedReadIterator.HitsForInsert hits = alignedIterator.next();
            hits.filterReads(alignmentFilter);
            if (hits.numHits() > 0) {
                return hits;
            }
        }
        return null;
    }

    /**
     * Copies alignment info from aligned to unaligned read, clips as appropriate, and sets PG ID.
     * @param unaligned Original SAMRecord, and object into which values are copied.
     * @param aligned Holds alignment info that will be copied into unaligned.
     */
    private void transferAlignmentInfoToFragment(SAMRecord unaligned, SAMRecord aligned) {
        setValuesFromAlignment(unaligned, aligned);
        updateCigarForTrimmedOrClippedBases(unaligned, aligned);
        if (this.programRecord != null) {
            unaligned.setAttribute(ReservedTagConstants.PROGRAM_GROUP_ID,
                this.programRecord.getProgramGroupId());
        }
    }

    /**
     * Copies alignment info from aligned to unaligned read, if there is an alignment, and sets mate information.
     * @param firstUnaligned Original first of pair, into which alignment and pair info will be written.
     * @param secondUnaligned Original second of pair, into which alignment and pair info will be written.
     * @param firstAligned Aligned first of pair, or null if no alignment.
     * @param secondAligned Aligned second of pair, or null if no alignment.
     */
    private void transferAlignmentInfoToPairedRead(SAMRecord firstUnaligned, SAMRecord secondUnaligned, SAMRecord firstAligned, SAMRecord secondAligned) {
        if (firstAligned != null) transferAlignmentInfoToFragment(firstUnaligned, firstAligned);
        if (secondAligned != null) transferAlignmentInfoToFragment(secondUnaligned, secondAligned);
        if (isClipOverlappingReads()) clipForOverlappingReads(firstUnaligned, secondUnaligned);
        SamPairUtil.setProperPairAndMateInfo(secondUnaligned, firstUnaligned, header, expectedOrientations);
    }



    /**
     * Checks to see whether the ends of the reads overlap and soft clips reads
     * them if necessary.
     */
    protected void clipForOverlappingReads(final SAMRecord read1, final SAMRecord read2) {
        // If both reads are mapped, see if we need to clip the ends due to small
        // insert size
        if (!(read1.getReadUnmappedFlag() || read2.getReadUnmappedFlag())) {

            if (read1.getReadNegativeStrandFlag() != read2.getReadNegativeStrandFlag())
            {
                final SAMRecord pos = (read1.getReadNegativeStrandFlag()) ? read2 : read1;
                final SAMRecord neg = (read1.getReadNegativeStrandFlag()) ? read1 : read2;

                // Innies only -- do we need to do anything else about jumping libraries?
                if (pos.getAlignmentStart() < neg.getAlignmentEnd()) {
                    final int posDiff = pos.getAlignmentEnd() - neg.getAlignmentEnd();
                    final int negDiff = pos.getAlignmentStart() - neg.getAlignmentStart();

                    if (posDiff > 0) {
                        CigarUtil.softClip3PrimeEndOfRead(pos, Math.min(pos.getReadLength(),
                                pos.getReadLength() - posDiff + 1));
                    }

                    if (negDiff > 0) {
                        CigarUtil.softClip3PrimeEndOfRead(neg, Math.min(neg.getReadLength(),
                                neg.getReadLength() - negDiff + 1));
                    }

                }
            }
            else {
                // TODO: What about RR/FF pairs?
            }
         }

    }

    /**
     * Sets the values from the alignment record on the unaligned BAM record.  This
     * preserves all data from the unaligned record (ReadGroup, NoiseRead status, etc)
     * and adds all the alignment info 
     *
     * @param rec           The unaligned read record
     * @param alignment     The alignment record
     */
    protected void setValuesFromAlignment(final SAMRecord rec, final SAMRecord alignment) {
        for (final SAMRecord.SAMTagAndValue attr : alignment.getAttributes()) {
            // Copy over any non-reserved attributes.
            if (!isReservedTag(attr.tag) || this.attributesToRetain.contains(attr.tag)) {
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
            SAMRecordUtil.reverseComplement(rec);
        }

    }

    protected void updateCigarForTrimmedOrClippedBases(final SAMRecord rec, final SAMRecord alignment) {

        // If the read maps off the end of the alignment, clip it
        SAMSequenceRecord refseq = rec.getHeader().getSequence(rec.getReferenceIndex());
        if (rec.getAlignmentEnd() > refseq.getSequenceLength()) {
            // 1-based index of first base in read to clip.
            int clipFrom = refseq.getSequenceLength() - rec.getAlignmentStart() + 1;
            List<CigarElement> newCigarElements  = CigarUtil.softClipEndOfRead(clipFrom, rec.getCigar().getCigarElements());
            rec.setCigar(new Cigar(newCigarElements));
        }

        // If the read was trimmed or not all the bases were sent for alignment, clip it
        int alignmentReadLength = alignment.getReadLength();
        int originalReadLength = rec.getReadLength();
        int trimmed = (!rec.getReadPairedFlag()) || rec.getFirstOfPairFlag()
                ? this.read1BasesTrimmed != null ? this.read1BasesTrimmed : 0
                : this.read2BasesTrimmed != null ? this.read2BasesTrimmed : 0;
        int notWritten = originalReadLength - (alignmentReadLength + trimmed);

        rec.setCigar(CigarUtil.addSoftClippedBasesToEndsOfCigar(
            rec.getCigar(), rec.getReadNegativeStrandFlag(), notWritten, trimmed));

        // If the adapter sequence is marked and clipAdapter is true, clip it
        if (this.clipAdapters && rec.getAttribute(ReservedTagConstants.XT) != null){
            CigarUtil.softClip3PrimeEndOfRead(rec, rec.getIntegerAttribute(ReservedTagConstants.XT));
        }
    }


    protected SAMSequenceDictionary getSequenceDictionary() { return this.sequenceDictionary; }

    protected SAMProgramRecord getProgramRecord() { return this.programRecord; }

    protected void setProgramRecord(SAMProgramRecord pg ) {
        if (this.programRecord != null) {
            throw new IllegalStateException("Cannot set program record more than once on alignment merger.");
        }
        this.programRecord = pg;
        this.header.addProgramRecord(pg);
        SAMUtils.chainSAMProgramRecord(header, pg);
    }

    protected boolean isReservedTag(String tag) {
        final char firstCharOfTag = tag.charAt(0);

        // All tags that start with a lower-case letter are user defined and should not be overridden by aligner
        // unless explicitly specified in attributesToRetain.
        if (Character.isLowerCase(firstCharOfTag)) return true;

        for (char c : RESERVED_ATTRIBUTE_STARTS) {
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

    public void setClipOverlappingReads(boolean clipOverlappingReads) {
        this.clipOverlappingReads = clipOverlappingReads;
    }
}
