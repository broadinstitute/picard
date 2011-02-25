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
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.CigarUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.SortingCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

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

    private static final String RESERVED_ATTRIBUTE_STARTS = "XYZ";

    private final Log log = Log.getInstance(AbstractAlignmentMerger.class);
    private final File unmappedBamFile;
    private final File targetBamFile;
    private SAMSequenceDictionary sequenceDictionary = null;
    private ReferenceSequenceFileWalker refSeq = null;
    private final boolean clipAdapters;
    private final boolean bisulfiteSequence;
    private SAMProgramRecord programRecord;
    private final boolean jumpingLibrary;
    private final boolean alignedReadsOnly;
    private final SAMFileHeader header;
    private final List<String> attributesToRetain = new ArrayList<String>();
    private final File referenceFasta;
    private final Integer read1BasesTrimmed;
    private final Integer read2BasesTrimmed;


    protected abstract CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords();
    
    protected boolean ignoreAlignment(SAMRecord sam) { return false; } // default implementation

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map.  Required.
     * @param targetBamFile     The file to which to write the merged SAM records. Required.
     * @param referenceFasta    The reference sequence for the map files. Optional.  If the referenceFasta
     *                          is present, then it is required that a sequence dictionary with the same
     *                          name but a ".dict" extension be present in the same directory as the
     *                          fasta file.
     * @param clipAdapters      Whether adapters marked in unmapped BAM file are clipped from the read.
     *                          Required.
     * @param bisulfiteSequence Whether the reads are bisulfite sequence (used when calculating the
     *                          NM and UQ tags). Required.
     * @param jumpingLibrary    Whether this is a jumping library.  Required.
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param programRecord     Program record for taget file SAMRecords created. Required.
     * @param attributesToRetain  private attributes from the alignment record that should be
     *                          included when merging.  This overrides the exclusion of
     *                          attributes whose tags start with the reserved characters
     *                          of X, Y, and Z
     * @param read1BasesTrimmed The number of bases trimmed from read 1 prior to alignment.  Optional.
     * @param read2BasesTrimmed The number of bases trimmed from read 2 prior to alignment.  Optional.
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean jumpingLibrary,
                                   final boolean alignedReadsOnly, final SAMProgramRecord programRecord,
                                   final List<String> attributesToRetain, final Integer read1BasesTrimmed,
                                   final Integer read2BasesTrimmed) {
        this.unmappedBamFile = unmappedBamFile;
        this.targetBamFile = targetBamFile;
        this.referenceFasta = referenceFasta;

        IoUtil.assertFileIsReadable(unmappedBamFile);
        IoUtil.assertFileIsWritable(targetBamFile);
        if (referenceFasta != null) {
            setRefSeqFileWalker();
            String fastaPath = referenceFasta.getAbsolutePath();
            File sd = new File(fastaPath.substring(0, fastaPath.lastIndexOf(".")) + ".dict");
            if (!sd.exists()) {
                throw new PicardException("No sequence dictionary was found for the reference " +
                    fastaPath + ".  A sequence dictionary is required for alignment merging.");
            }
            sequenceDictionary =
                    new SAMFileReader(IoUtil.openFileForReading(sd)).getFileHeader().getSequenceDictionary();
        }
        this.clipAdapters = clipAdapters;
        this.bisulfiteSequence = bisulfiteSequence;
        this.jumpingLibrary = jumpingLibrary;
        this.alignedReadsOnly = alignedReadsOnly;
        this.programRecord = programRecord;

        this.header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        if (programRecord != null) {
            header.addProgramRecord(programRecord);
        }
        header.setSequenceDictionary(this.sequenceDictionary);
        if (attributesToRetain != null) {
            this.attributesToRetain.addAll(attributesToRetain);
        }
        this.read1BasesTrimmed = read1BasesTrimmed;
        this.read2BasesTrimmed = read2BasesTrimmed;

    }

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map.  Required.
     * @param targetBamFile     The file to which to write the merged SAM records. Required.
     * @param referenceFasta    The reference sequence for the map files. Optional.  If the referenceFasta
     *                          is present, then it is required that a sequence dictionary with the same
     *                          name but a ".dict" extension be present in the same directory as the
     *                          fasta file.
     * @param clipAdapters      Whether adapters marked in unmapped BAM file are clipped from the read.
     *                          Required.
     * @param bisulfiteSequence Whether the reads are bisulfite sequence (used when calculating the
     *                          NM and UQ tags). Required.
     * @param jumpingLibrary    Whether this is a jumping library.  Required.
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param programRecord     Program record for taget file SAMRecords created. Required.
     * @param attributesToRetain  private attributes from the alignment record that should be
     *                          included when merging.  This overrides the exclusion of
     *                          attributes whose tags start with the reserved characters
     *                          of X, Y, and Z
     * @deprecated              Retained for backwards-compatibility only.
     */
    @Deprecated
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean jumpingLibrary,
                                   final boolean alignedReadsOnly, final SAMProgramRecord programRecord,
                                   final List<String> attributesToRetain) {
        this(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
             jumpingLibrary, alignedReadsOnly, programRecord, attributesToRetain, null, null);

    }

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map.  Required.
     * @param targetBamFile     The file to which to write the merged SAM records. Required.
     * @param referenceFasta    The reference sequence for the map files. Optional.  If the referenceFasta
     *                          is present, then it is required that a sequence dictionary with the same
     *                          name but a ".dict" extension be present in the same directory as the
     *                          fasta file.
     * @param clipAdapters      Whether adapters marked in unmapped BAM file are clipped from the read.
     *                          Required.
     * @param bisulfiteSequence Whether the reads are bisulfite sequence (used when calculating the
     *                          NM and UQ tags). Required.
     * @param jumpingLibrary    Whether this is a jumping library.  Required.
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param programRecord     Program record for taget file SAMRecords created. Required.
     * @deprecated              Retained for backwards-compatibility only.
     */
    @Deprecated
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean jumpingLibrary,
                                   final boolean alignedReadsOnly, final SAMProgramRecord programRecord) {
        this(unmappedBamFile, targetBamFile, referenceFasta, clipAdapters, bisulfiteSequence,
             jumpingLibrary, alignedReadsOnly, programRecord, null, null, null);
    }


    /**
     * Merges the alignment data with the non-aligned records from the source BAM file.
     */
    public void mergeAlignment() {

        final SAMRecordQueryNameComparator comparator = new SAMRecordQueryNameComparator();

        // Open the file of unmapped records and write the read groups to the the header for the merged file
        final SAMFileReader unmappedSam = new SAMFileReader(IoUtil.openFileForReading(this.unmappedBamFile));
        final CloseableIterator<SAMRecord> unmappedIterator = unmappedSam.iterator();
        header.setReadGroups(unmappedSam.getFileHeader().getReadGroups());

        int aligned = 0;
        int unmapped = 0;

        // Get the aligned records and set up the first one
        final CloseableIterator<SAMRecord> alignedIterator = getQuerynameSortedAlignedRecords();
        SAMRecord nextAligned = alignedIterator.hasNext() ? alignedIterator.next() : null;

        // Create the sorting collection that will write the records in coordinate order
        // to the final bam file
        final SortingCollection<SAMRecord> coordinateSorted = SortingCollection.newInstance(
            SAMRecord.class, new BAMRecordCodec(header), new SAMRecordCoordinateComparator(),
            MAX_RECORDS_IN_RAM);

        // A wrapper around the sorting collection that makes sure that reads get their mate
        // information properly set
        final ClippedPairFixer pairFixer = new ClippedPairFixer(coordinateSorted, header, alignedReadsOnly);


        while (unmappedIterator.hasNext()) {
            final SAMRecord rec = unmappedIterator.next();

            if (nextAligned != null && comparator.compare(rec, nextAligned) > 0) {
                throw new IllegalStateException("Aligned record iterator (" + nextAligned.getReadName() +
                        ") is behind the unmapped reads (" + rec.getReadName() + ")");
            }
            rec.setReadName(cleanReadName(rec.getReadName()));
            rec.setHeader(header);

            // If the next record is a match and is an acceptable alignment, pull the info over to the unmapped record
            if (isMatch(rec, nextAligned)) {
                if (!(nextAligned.getReadUnmappedFlag() || ignoreAlignment(nextAligned))) {
                    setValuesFromAlignment(rec, nextAligned);
                    if (programRecord != null) {
                        rec.setAttribute(ReservedTagConstants.PROGRAM_GROUP_ID,
                            programRecord.getProgramGroupId());
                    }
                    aligned++;
                }
                else {
                    unmapped++;
                }
                nextAligned = alignedIterator.hasNext() ? alignedIterator.next() : null;
            }
            else {
                unmapped++;
            }

            // Add it if either the read or its mate are mapped, unless we are adding aligned reads only
            final boolean eitherReadMapped = !rec.getReadUnmappedFlag() || (rec.getReadPairedFlag() && !rec.getMateUnmappedFlag());

            if (eitherReadMapped || !alignedReadsOnly) {
                pairFixer.add(rec);
            }
        }
        unmappedIterator.close();
        if (alignedIterator.hasNext()) {
            throw new IllegalStateException("Reads remaining on alignment iterator: " + alignedIterator.next().getReadName() + "!");
        }
        alignedIterator.close();

        // Write the records to the output file in coordinate sorted order,
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, this.targetBamFile);
        int count = 0;
        CloseableIterator<SAMRecord> it = coordinateSorted.iterator();
        while (it.hasNext()) {
            SAMRecord rec = it.next();
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
                log.info(count + " SAMRecords written to " + targetBamFile.getName());
            }
        }
        writer.close();


        log.info("Wrote " + aligned + " alignment records and " + (alignedReadsOnly ? 0 : unmapped) + " unmapped reads.");
    }


    /**
     * Determines whether two SAMRecords represent the same read
     */
    protected boolean isMatch(final SAMRecord unaligned, final SAMRecord aligned) {
        return (aligned != null &&
                aligned.getReadName().equals(unaligned.getReadName()) &&
                (unaligned.getReadPairedFlag() == false ||
                aligned.getFirstOfPairFlag() == unaligned.getFirstOfPairFlag()));
    }

    /**
     * Sets the values from the alignment record on the unaligned BAM record.  This
     * preserves all data from the unaligned record (ReadGroup, NoiseRead status, etc)
     * and adds all the alignment info from Maq
     *
     * @param rec           The unaligned read record
     * @param alignment     The alignment record
     */
    protected void setValuesFromAlignment(final SAMRecord rec, final SAMRecord alignment) {
        for (final SAMRecord.SAMTagAndValue attr : alignment.getAttributes()) {
            // Copy over any non-reserved attributes.
            if (RESERVED_ATTRIBUTE_STARTS.indexOf(attr.tag.charAt(0)) == -1
                    || attributesToRetain.contains(attr.tag)) {
                rec.setAttribute(attr.tag, attr.value);
            }
        }
        rec.setReadUnmappedFlag(alignment.getReadUnmappedFlag());
        rec.setReferenceName(alignment.getReferenceName());
        rec.setAlignmentStart(alignment.getAlignmentStart());
        rec.setReadNegativeStrandFlag(alignment.getReadNegativeStrandFlag());
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

        // Factor the trimmed bases into the cigar
        int basesAligned = alignment.getReadLength();
        int readLength = rec.getReadLength();
        int trimmed = (!rec.getReadPairedFlag()) || rec.getFirstOfPairFlag()
                ? read1BasesTrimmed != null ? read1BasesTrimmed : 0
                : read2BasesTrimmed != null ? read2BasesTrimmed : 0;
        int notWritten = readLength - (basesAligned + trimmed);

        rec.setCigar(CigarUtil.addSoftClippedBasesToEndsOfCigar(
            rec.getCigar(), rec.getReadNegativeStrandFlag(), notWritten, trimmed));

        if (clipAdapters && rec.getAttribute(ReservedTagConstants.XT) != null){
            CigarUtil.softClip3PrimeEndOfRead(rec, rec.getIntegerAttribute(ReservedTagConstants.XT));
        }
    }


    /**
     * Wrapper around the sorting collection to make sure reads get their mate information
     * set properly after clipping.
     */
    private static class ClippedPairFixer {

        private final SortingCollection<SAMRecord> collection;
        private final SAMFileHeader header;
        private final boolean alignedOnly;
        private SAMRecord pending = null;

        public ClippedPairFixer(final SortingCollection<SAMRecord> collection, final SAMFileHeader header,
                                boolean alignedReadsOnly) {
            this.collection = collection;
            this.header = header;
            this.alignedOnly = alignedReadsOnly;
        }

        public void add(final SAMRecord record) {
            if (!record.getReadPairedFlag()) {
                collection.add(record);
            }
            else if (pending == null) {
                pending = record;
            }
            else if (alignedOnly && !pending.getReadName().equals(record.getReadName())) {
                collection.add(pending);
                pending = record;
            }
            else {
                if (!record.getReadName().equals(pending.getReadName())) {
                    throw new PicardException("Non-paired reads: " + record.getReadName() +
                            ", " + pending.getReadName());
                }

                // If both reads are mapped, see if we need to clip the ends due to small
                // insert size
                if (!(record.getReadUnmappedFlag() || pending.getReadUnmappedFlag())) {

                    if (record.getReadNegativeStrandFlag() != pending.getReadNegativeStrandFlag())
                    {
                        final SAMRecord pos = (record.getReadNegativeStrandFlag()) ? pending : record;
                        final SAMRecord neg = (record.getReadNegativeStrandFlag()) ? record : pending;

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


                SamPairUtil.setMateInfo(record, pending, header);
                collection.add(pending);
                collection.add(record);
                pending = null;
            }
        }


    }

    /**
     * Strips read name of extraneous read number extensions
     */
    protected String cleanReadName(String readName) {
        if (readName.endsWith("/1") || readName.endsWith("/2")) {
            readName = readName.substring(0, readName.length()-2);
        }
        return readName;
    }

    protected SAMSequenceDictionary getSequenceDictionary() { return this.sequenceDictionary; }
    protected SAMProgramRecord getProgramRecord() { return this.programRecord; }
    protected void setProgramRecord(SAMProgramRecord pg ) {
        this.programRecord = pg;
        header.addProgramRecord(pg);
    }
    protected boolean isJumpingLibrary() { return this.jumpingLibrary; }
    protected SAMFileHeader getHeader() { return this.header; }

    protected void setRefSeqFileWalker() {
        if (this.referenceFasta != null && this.referenceFasta.exists()) {
            refSeq = new ReferenceSequenceFileWalker(referenceFasta);
        }
    }
}
