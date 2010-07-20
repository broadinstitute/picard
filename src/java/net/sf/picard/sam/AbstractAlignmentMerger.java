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
import java.util.TreeSet;

/**
 * Abstract class that coordinates the general task of taking in a set of alignment information,
 * possibly in SAM format, possibly in other formats, and merging that with the set of all reads
 * for which alignment was attempted, stored in an unmapped SAM file.
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


    protected abstract CloseableIterator<SAMRecord> getQuerynameSortedAlignedRecords();
    

    /**
     * Constructor
     *
     * @param unmappedBamFile   The BAM file that was used as the input to the Maq aligner, which will
     *                          include info on all the reads that did not map
     * @param targetBamFile     The file to which to write the merged SAM records
     * @param referenceFasta    The reference sequence for the map files
     * @param clipAdapters      Whether adapters marked in unmapped BAM file are clipped from the read
     * @param bisulfiteSequence Whether the reads are bisulfite sequence
     * @param jumpingLibrary    Whether this is a jumping library
     * @param alignedReadsOnly  Whether to output only those reads that have alignment data
     * @param programRecord     Program record for taget file SAMRecords created.
     */
    public AbstractAlignmentMerger(final File unmappedBamFile, final File targetBamFile,
                                   final File referenceFasta, final boolean clipAdapters,
                                   final boolean bisulfiteSequence, final boolean jumpingLibrary,
                                   final boolean alignedReadsOnly, final SAMProgramRecord programRecord) {
        this.unmappedBamFile = unmappedBamFile;
        this.targetBamFile = targetBamFile;

        IoUtil.assertFileIsReadable(unmappedBamFile);
        IoUtil.assertFileIsWritable(targetBamFile);
        if (referenceFasta != null) {
            if (referenceFasta.exists()) {
                refSeq = new ReferenceSequenceFileWalker(referenceFasta);
            }
            String fastaPath = referenceFasta.getAbsolutePath();
            File sd = new File(fastaPath.substring(0, fastaPath.lastIndexOf(".")) + ".dict");
            sequenceDictionary = sd.exists()
                ? new SAMFileReader(IoUtil.openFileForReading(sd)).getFileHeader().getSequenceDictionary()
                : null;
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

    }

    /**
     * Merges the alignment from the map file with the non-aligned records from the source BAM file.
     */
    public void mergeAlignment() {

        SAMFileReader unmappedSam = null;
        if (this.unmappedBamFile != null) {
            unmappedSam = new SAMFileReader(IoUtil.openFileForReading(this.unmappedBamFile));
        }

        // Write the read groups to the header
        if (unmappedSam != null) header.setReadGroups(unmappedSam.getFileHeader().getReadGroups());

        int aligned = 0;
        int unmapped = 0;

        final PeekableIterator<SAMRecord> alignedIterator = 
                new PeekableIterator(getQuerynameSortedAlignedRecords());
        final SortingCollection<SAMRecord> alignmentSorted = SortingCollection.newInstance(
            SAMRecord.class, new BAMRecordCodec(header), new SAMRecordCoordinateComparator(),
            MAX_RECORDS_IN_RAM);
        final ClippedPairFixer pairFixer = new ClippedPairFixer(alignmentSorted, header, alignedReadsOnly);

        final CloseableIterator<SAMRecord> unmappedIterator = unmappedSam.iterator();
        SAMRecord nextAligned = alignedIterator.hasNext() ? alignedIterator.next() : null;
        SAMRecord lastAligned = null;

        final UnmappedReadSorter unmappedSorter = new UnmappedReadSorter(unmappedIterator);
        while (unmappedSorter.hasNext()) {
            final SAMRecord rec = unmappedSorter.next();
            rec.setReadName(cleanReadName(rec.getReadName()));
            if (nextAligned != null && rec.getReadName().compareTo(nextAligned.getReadName()) > 0) {
                throw new PicardException("Aligned Record iterator (" + nextAligned.getReadName() +
                        ") is behind the umapped reads (" + rec.getReadName() + ")");
            }
            rec.setHeader(header);

            if (isMatch(rec, nextAligned)) {
                if (!ignoreAlignment(nextAligned)) {
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
            throw new PicardException("Reads remaining on alignment iterator!");
        }
        alignedIterator.close();

        final SAMFileWriter writer =
                new SAMFileWriterFactory().makeBAMWriter(header, true, this.targetBamFile);
        int count = 0;
        CloseableIterator<SAMRecord> it = alignmentSorted.iterator();
        while (it.hasNext()) {
            SAMRecord rec = it.next();
            if (!rec.getReadUnmappedFlag()) {
                if (refSeq != null) {
                    byte referenceBases[] = refSeq.get(rec.getReferenceIndex()).getBases();
                    rec.setAttribute(SAMTag.NM.name(),
                        SequenceUtil.calculateSamNmTag(rec, referenceBases, 0, bisulfiteSequence));
                    rec.setAttribute(SAMTag.UQ.name(),
                        SequenceUtil.sumQualitiesOfMismatches(rec, referenceBases, 0, bisulfiteSequence));
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
     * @param alignment     The Maq alignment record
     */
    protected void setValuesFromAlignment(final SAMRecord rec, final SAMRecord alignment) {
        for (final SAMRecord.SAMTagAndValue attr : alignment.getAttributes()) {
            // Copy over any non-reserved attributes.
            if (RESERVED_ATTRIBUTE_STARTS.indexOf(attr.tag.charAt(0)) == -1) {
                rec.setAttribute(attr.tag, attr.value);
            }
        }
        rec.setReadUnmappedFlag(alignment.getReadUnmappedFlag());
        rec.setReferenceIndex(alignment.getReferenceIndex());
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
            rec.setInferredInsertSize(alignment.getInferredInsertSize());
            rec.setMateUnmappedFlag(alignment.getMateUnmappedFlag());
            rec.setMateReferenceIndex(alignment.getMateReferenceIndex());
            rec.setMateAlignmentStart(alignment.getMateAlignmentStart());
            rec.setMateNegativeStrandFlag(alignment.getMateNegativeStrandFlag());
        }

        // If it's on the negative strand, reverse complement the bases
        // and reverse the order of the qualities
        if (rec.getReadNegativeStrandFlag()) {
            SAMRecordUtil.reverseComplement(rec);
        }

        if (clipAdapters && rec.getAttribute(ReservedTagConstants.XT) != null){
            CigarUtil.softClip3PrimeEndOfRead(rec, rec.getIntegerAttribute(ReservedTagConstants.XT));
        }
    }

    /**
     * Temporary class to wrap around unaligned BAM files which may be sorted using the *old*
     * SAMRecordQueryNameComparator (which did not use read ordering, and so first and second read
     * are not necessarily in order).  Keeps the next three records in a sorted set and returns the
     * next one from it rather than the order un the unmapped reads iterator.
     */
    protected static class UnmappedReadSorter  {

        private final TreeSet<SAMRecord> reallySorted = new TreeSet<SAMRecord>(new SAMRecordQueryNameComparator());
        private final CloseableIterator<SAMRecord> unmappedReads;

        public UnmappedReadSorter(final CloseableIterator<SAMRecord> unmappedReads) {
            this.unmappedReads = unmappedReads;
            // Add the first three records to our sorter
            for (int i = 0; i < 3; i++) {
                if (unmappedReads.hasNext()) {
                    reallySorted.add(unmappedReads.next());
                }
            }
        }

        public boolean hasNext() {
            return reallySorted.size() > 0;
        }

        public SAMRecord next() {
            if (reallySorted.size() == 0) {
                throw new ArrayIndexOutOfBoundsException("No such element in UnmappedReadSorter!");
            }
            final SAMRecord next = reallySorted.first();
            reallySorted.remove(next);
            if (unmappedReads.hasNext()) {
                reallySorted.add(unmappedReads.next());
            }
            return next;
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
    protected boolean ignoreAlignment(SAMRecord sam) { return false; } // default implementation

}
