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

import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CoordMath;
import net.sf.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
/**
 * Factory class for creating SAMRecords for testing purposes. Various methods can be called
 * to add new SAM records (or pairs of records) to a list which can then be returned at
 * any point. The records must reference human chromosomes (excluding randoms etc.).
 *
 * Although this is a class for testing, it is in the src tree because it is included in the sam jarfile.
 *
 * @author Tim Fennell
 */
public class SAMRecordSetBuilder implements Iterable<SAMRecord> {
    private static final String[] chroms = {
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
            "chr21", "chr22", "chrX", "chrY", "chrM"
    };
    private static final byte[] BASES = {'A','C','G','T'};
    private static final String READ_GROUP_ID = "1";
    private static final String SAMPLE = "FREE_SAMPLE";
    private final Random random = new Random();

    private final SAMFileHeader header;
    private final Collection<SAMRecord> records;

    private int readLength = 36 ;

    private SAMProgramRecord programRecord = null;
    private SAMReadGroupRecord readGroup = null;

    private static final int DEFAULT_CHROMOSOME_LENGTH = 100000000;

    /**
     * Constructs a new SAMRecordSetBuilder with all the data needed to keep the records
     * sorted in coordinate order.
     */
    public SAMRecordSetBuilder() {
        this(true, SAMFileHeader.SortOrder.coordinate);
    }

    /**
     * Construct a new SAMRecordSetBuilder.
     * @para
     * m sortForMe If true, keep the records created in sorted order.
     * @param sortOrder If sortForMe, defines the sort order.
     */
    public SAMRecordSetBuilder(final boolean sortForMe, final SAMFileHeader.SortOrder sortOrder) {
        this(sortForMe, sortOrder, true) ;
    }

    public SAMRecordSetBuilder(final boolean sortForMe, final SAMFileHeader.SortOrder sortOrder, final boolean addReadGroup) {
        this(sortForMe, sortOrder, addReadGroup, DEFAULT_CHROMOSOME_LENGTH);
    }

    public SAMRecordSetBuilder(final boolean sortForMe, final SAMFileHeader.SortOrder sortOrder, final boolean addReadGroup, final int defaultChromosomeLength) {
        final List<SAMSequenceRecord> sequences = new ArrayList<SAMSequenceRecord>();
        for (final String chrom : chroms) {
            final SAMSequenceRecord sequenceRecord = new SAMSequenceRecord(chrom, defaultChromosomeLength);
            sequences.add(sequenceRecord);
        }

        this.header = new SAMFileHeader();
        this.header.setSequenceDictionary(new SAMSequenceDictionary(sequences));
        this.header.setSortOrder(sortOrder);
        if (sortForMe) {
            final SAMRecordComparator comparator;
            if (sortOrder == SAMFileHeader.SortOrder.queryname) {
                comparator = new SAMRecordQueryNameComparator();
            } else {
                comparator = new SAMRecordCoordinateComparator();
            }
            this.records = new TreeSet<SAMRecord>(comparator);
        } else {
            this.records = new ArrayList<SAMRecord>();
        }

        if (addReadGroup) {
            final SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord(READ_GROUP_ID);
            readGroupRecord.setSample(SAMPLE);
            readGroupRecord.setPlatform("ILLUMINA");
            final List<SAMReadGroupRecord> readGroups = new ArrayList<SAMReadGroupRecord>();
            readGroups.add(readGroupRecord);
            this.header.setReadGroups(readGroups);
        }
    }

    /**
     * Set the seed of the random number generator for cases in which repeatable result is desired.
     * @param seed
     */
    public void setRandomSeed(final long seed) {
        random.setSeed(seed);
    }

    /**
     * Adds the given program record to the header, and assigns the PG tag to any SAMRecords
     * created after it has been added. May be called multiple times in order to assign different
     * PG IDs to different SAMRecords.  programRecord may be null to stop assignment of PG tag.
     * It is up to the caller to ensure that program record IDs do not collide.
     */
    public void setProgramRecord(final SAMProgramRecord programRecord) {
        this.programRecord = programRecord;
        if (programRecord != null) {
            this.header.addProgramRecord(programRecord);
        }
    }

    public void setReadGroup(final SAMReadGroupRecord readGroup) {
        this.readGroup = readGroup;
        if (readGroup != null) {
            this.header.addReadGroup(readGroup);
        }
    }

    /** Returns the accumulated list of sam records. */
    public Collection<SAMRecord> getRecords() { return this.records; }

    /** Returns a CloseableIterator over the collection of SAMRecords. */
    public CloseableIterator<SAMRecord> iterator() {
        return new CloseableIterator<SAMRecord>() {
            private final Iterator<SAMRecord> iterator = records.iterator();
            public void close() { /** Do nothing. */  }
            public boolean hasNext() { return this.iterator.hasNext(); }
            public SAMRecord next() { return this.iterator.next(); }
            public void remove() { this.iterator.remove(); }
        };
    }

    /**
     * Adds a fragment record (mapped or unmapped) to the set using the provided contig start and optionally the strand,
     * cigar string, quality string or default quality score.  This does not modify the flag field, which should be updated
     * if desired before adding the return to the list of records.
     */
    private SAMRecord createReadNoFlag(final String name, final int contig, final int start, final boolean negativeStrand,
                                       final boolean recordUnmapped, final String cigar, final String qualityString,
                                       final int defaultQuality) throws SAMException {
        final SAMRecord rec = new SAMRecord(this.header);
        rec.setReadName(name);
        if (chroms.length <= contig) {
            throw new SAMException("Contig too big [" + chroms.length + " < " + contig);
        }
        if (0 <= contig) {
            rec.setReferenceIndex(contig);
            rec.setReferenceName(chroms[contig]);
            rec.setAlignmentStart(start);
        }
        if (!recordUnmapped) {
            rec.setReadNegativeStrandFlag(negativeStrand);
            if (null != cigar) {
                rec.setCigarString(cigar);
                readLength = rec.getCigar().getReadLength();
            } else if (!rec.getReadUnmappedFlag()) {
                rec.setCigarString(readLength + "M");
            }
            rec.setMappingQuality(255);
        } else {
            rec.setReadUnmappedFlag(true);
        }
        rec.setAttribute(SAMTag.RG.name(), READ_GROUP_ID);

        if (programRecord != null) {
            rec.setAttribute(SAMTag.PG.name(), programRecord.getProgramGroupId());
        }

        if (readGroup != null) {
            rec.setAttribute(SAMTag.RG.name(), readGroup.getReadGroupId());
        }

        fillInBasesAndQualities(rec, qualityString, defaultQuality);

        return rec;
    }

    /**
     * Adds a skeletal fragment (non-PE) record to the set using the provided
     * contig start and strand information.
     */
    public void addFrag(final String name, final int contig, final int start, final boolean negativeStrand) {
        addFrag(name, contig, start, negativeStrand, false, null, null, -1);
    }

    /**
     * Adds a fragment record (mapped or unmapped) to the set using the provided contig start and optionally the strand,
     * cigar string, quality string or default quality score.
     */
    public SAMRecord addFrag(final String name, final int contig, final int start, final boolean negativeStrand,
                             final boolean recordUnmapped, final String cigar, final String qualityString,
                             final int defaultQuality) throws SAMException {
        final SAMRecord rec = createReadNoFlag(name, contig, start, negativeStrand, recordUnmapped, cigar, qualityString, defaultQuality);
        this.records.add(rec);
        return rec;
    }

    /**
     * Fills in the bases and qualities for the given record. Quality data is randomly generated if the defaultQuality
     * is set to -1. Otherwise all qualities will be set to defaultQuality. If a quality string is provided that string
     * will be used instead of the defaultQuality.
     */
    private void fillInBasesAndQualities(final SAMRecord rec, final String qualityString, final int defaultQuality) {

        if (null == qualityString) {
            fillInBasesAndQualities(rec, defaultQuality);
        } else {
            fillInBases(rec);
            rec.setBaseQualityString(qualityString);
        }
    }

    /**
     * Randomly fills in the bases for the given record.
     */
    private void fillInBases(final SAMRecord rec){
        final int length = this.readLength;
        final byte[] bases = new byte[length];

        for (int i = 0; i < length; ++i) {
            bases[i] = BASES[this.random.nextInt(BASES.length)];
        }

        rec.setReadBases(bases);
    }

    /**
     * Adds an unmapped fragment read to the builder.
     */
    public void addUnmappedFragment(final String name) {
        addFrag(name, -1, -1, false, true, null, null, -1);
    }


    /**
     * Adds a skeletal pair of records to the set using the provided
     * contig starts.  The pair is assumed to be a well
     * formed pair sitting on a single contig.
     */
    public void addPair(final String name, final int contig, final int start1, final int start2) {
        final SAMRecord end1 = new SAMRecord(this.header);
        final SAMRecord end2 = new SAMRecord(this.header);
        final boolean end1IsFirstOfPair = this.random.nextBoolean();

        end1.setReadName(name);
        end1.setReferenceIndex(contig);
        end1.setAlignmentStart(start1);
        end1.setReadNegativeStrandFlag(false);
        end1.setCigarString(readLength + "M");
        end1.setMappingQuality(255);
        end1.setReadPairedFlag(true);
        end1.setProperPairFlag(true);
        end1.setMateReferenceIndex(contig);
        end1.setAttribute(SAMTag.MC.name(), readLength + "M");
        end1.setMateAlignmentStart(start2);
        end1.setMateNegativeStrandFlag(true);
        end1.setFirstOfPairFlag(end1IsFirstOfPair);
        end1.setSecondOfPairFlag(!end1IsFirstOfPair);
        end1.setInferredInsertSize((int) CoordMath.getLength(start1, CoordMath.getEnd(start2, this.readLength)));
        end1.setAttribute(SAMTag.RG.name(), READ_GROUP_ID);
        if (programRecord != null) {
            end1.setAttribute(SAMTag.PG.name(), programRecord.getProgramGroupId());
        }
        if (readGroup != null) {
            end1.setAttribute(SAMTag.RG.name(), readGroup.getReadGroupId());
        }
        fillInBasesAndQualities(end1);

        end2.setReadName(name);
        end2.setReferenceIndex(contig);
        end2.setAlignmentStart(start2);
        end2.setReadNegativeStrandFlag(true);
        end2.setCigarString(readLength + "M");
        end2.setMappingQuality(255);
        end2.setReadPairedFlag(true);
        end2.setProperPairFlag(true);
        end2.setMateReferenceIndex(contig);
        end2.setAttribute(SAMTag.MC.name(), readLength + "M");
        end2.setMateAlignmentStart(start1);
        end2.setMateNegativeStrandFlag(false);
        end2.setFirstOfPairFlag(!end1IsFirstOfPair);
        end2.setSecondOfPairFlag(end1IsFirstOfPair);
        end2.setInferredInsertSize(end1.getInferredInsertSize());
        end2.setAttribute(SAMTag.RG.name(), READ_GROUP_ID);
        if (programRecord != null) {
            end2.setAttribute(SAMTag.PG.name(), programRecord.getProgramGroupId());
        }
        if (readGroup != null) {
            end2.setAttribute(SAMTag.RG.name(), readGroup.getReadGroupId());
        }
        fillInBasesAndQualities(end2);

        this.records.add(end1);
        this.records.add(end2);
    }

    /**
     * Adds a pair of records (mapped or unmmapped) to the set using the provided contig starts.
     * The pair is assumed to be a well formed pair sitting on a single contig.
     */
    public List<SAMRecord> addPair(final String name, final int contig, final int start1, final int start2,
                                   final boolean record1Unmapped, final boolean record2Unmapped, final String cigar1,
                                   final String cigar2, final boolean strand1, final boolean strand2, final int defaultQuality) {
        final List<SAMRecord> recordsList = new LinkedList<SAMRecord>();

        final SAMRecord end1 = createReadNoFlag(name, contig, start1, strand1, record1Unmapped, cigar1, null, defaultQuality);
        final SAMRecord end2 = createReadNoFlag(name, contig, start2, strand2, record2Unmapped, cigar2, null, defaultQuality);

        end1.setReadPairedFlag(true);
        end1.setFirstOfPairFlag(true);

        if (!record1Unmapped && !record2Unmapped) {
            end1.setProperPairFlag(true);
            end2.setProperPairFlag(true);
        }
        end2.setReadPairedFlag(true);
        end2.setSecondOfPairFlag(true);

        // set mate info
        SamPairUtil.setMateInfo(end1, end2, header, true);

        recordsList.add(end1);
        recordsList.add(end2);

        records.add(end1);
        records.add(end2);

        return recordsList;
    }

    /**
     * Adds a pair with both ends unmapped to the builder.
     */
    public void addUnmappedPair(final String name) {
        final SAMRecord end1 = new SAMRecord(this.header);
        final SAMRecord end2 = new SAMRecord(this.header);
        final boolean end1IsFirstOfPair = this.random.nextBoolean();

        end1.setReadName(name);
        end1.setReadPairedFlag(false);
        end1.setReadUnmappedFlag(true);
        end1.setAttribute(SAMTag.MC.name(), null);
        end1.setProperPairFlag(false);
        end1.setFirstOfPairFlag(end1IsFirstOfPair);
        end1.setSecondOfPairFlag(!end1IsFirstOfPair);
        end1.setAttribute(SAMTag.RG.name(), READ_GROUP_ID);
        if (programRecord != null) {
            end1.setAttribute(SAMTag.PG.name(), programRecord.getProgramGroupId());
        }
        fillInBasesAndQualities(end1);

        end2.setReadName(name);
        end2.setReadPairedFlag(false);
        end2.setReadUnmappedFlag(true);
        end2.setAttribute(SAMTag.MC.name(), null);
        end2.setProperPairFlag(false);
        end2.setFirstOfPairFlag(!end1IsFirstOfPair);
        end2.setSecondOfPairFlag(end1IsFirstOfPair);
        end2.setAttribute(SAMTag.RG.name(), READ_GROUP_ID);
        if (programRecord != null) {
            end2.setAttribute(SAMTag.PG.name(), programRecord.getProgramGroupId());
        }
        fillInBasesAndQualities(end2);

        this.records.add(end1);
        this.records.add(end2);
    }

    /**
     * Fills in bases and qualities with randomly generated data.
     * Relies on the alignment start and end having been set to get read length.
     */
    private void fillInBasesAndQualities(final SAMRecord rec) {
        fillInBasesAndQualities(rec, -1);
    }

    /**
     * Fills in bases and qualities with a set default quality. If the defaultQuality is set to -1 quality scores will
     * be randomly generated.
     * Relies on the alignment start and end having been set to get read length.
     */
    private void fillInBasesAndQualities(final SAMRecord rec, final int defaultQuality) {
        final int length = this.readLength;
        final byte[] quals = new byte[length];

        if (-1 != defaultQuality) {
            Arrays.fill(quals, (byte) defaultQuality);
        } else {
            for (int i = 0; i < length; ++i) {
                quals[i] = (byte) this.random.nextInt(50);
            }
        }
        rec.setBaseQualities(quals);
        fillInBases(rec);
    }

    /**
     * Creates samFileReader from the data in instance of this class
     * @return SAMFileReader
     */
    public SAMFileReader getSamReader() {

        final File tempFile;

        try {
            tempFile = File.createTempFile("temp", ".sam");
        } catch (final IOException e) {
            throw new RuntimeIOException("problems creating tempfile", e);
        }



        this.header.setAttribute("VN", "1.0");
        final SAMFileWriter w = new SAMFileWriterFactory().makeBAMWriter(this.header, true, tempFile);
        for (final SAMRecord r:this.getRecords()){
            w.addAlignment(r);
        }


        w.close();

        final SAMFileReader reader = new SAMFileReader(tempFile);
        tempFile.deleteOnExit();

        return reader;
    }

    public SAMFileHeader getHeader() {
        return header;
    }
    public void setReadLength(final int readLength) { this.readLength = readLength; }

}
