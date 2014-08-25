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

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMRecordUtil;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.filter.FilteringIterator;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.FastqQualityFormat;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.QualityEncodingDetector;
import htsjdk.samtools.util.SolexaQualityConverter;
import htsjdk.samtools.util.SortingCollection;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Reverts a SAM file by optionally restoring original quality scores and by removing
 * all alignment information.
 */
public class RevertSam extends CommandLineProgram {
    @Usage public String USAGE = getStandardUsagePreamble() +
            "Reverts SAM or BAM files to a previous state by removing certain types of information and/or " +
            "substituting in the original quality scores when available.";

    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM/BAM file to revert the state of.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output SAM/BAM file to create.")
    public File OUTPUT;

    @Option(shortName="SO", doc="The sort order to create the reverted output file with.")
    public SortOrder SORT_ORDER = SortOrder.queryname;

    @Option(shortName=StandardOptionDefinitions.USE_ORIGINAL_QUALITIES_SHORT_NAME, doc="True to restore original qualities from the OQ field to the QUAL field if available.")
    public boolean RESTORE_ORIGINAL_QUALITIES = true;

    @Option(doc="Remove duplicate read flags from all reads.  Note that if this is true and REMOVE_ALIGNMENT_INFORMATION==false, " +
            " the output may have the unusual but sometimes desirable trait of having unmapped reads that are marked as duplicates.")
    public boolean REMOVE_DUPLICATE_INFORMATION = true;

    @Option(doc="Remove all alignment information from the file.")
    public boolean REMOVE_ALIGNMENT_INFORMATION = true;

    @Option(doc="When removing alignment information, the set of optional tags to remove.")
    public List<String> ATTRIBUTE_TO_CLEAR = new ArrayList<String>() {{
        add(SAMTag.NM.name());
        add(SAMTag.UQ.name());
        add(SAMTag.PG.name());
        add(SAMTag.MD.name());
        add(SAMTag.MQ.name());
        add(SAMTag.SA.name()); // Supplementary alignment metadata
        add(SAMTag.MC.name());      // Mate Cigar
    }};

    @Option(doc="WARNING: This option is potentially destructive. If enabled will discard reads in order to produce " +
            "a consistent output BAM. Reads discarded include (but are not limited to) paired reads with missing " +
            "mates, duplicated records, records with mismatches in length of bases and qualities. This option can " +
            "only be enabled if the output sort order is queryname and will always cause sorting to occur.")
    public boolean SANITIZE = false;

    @Option(doc="If SANITIZE=true and higher than MAX_DISCARD_FRACTION reads are discarded due to sanitization then" +
            "the program will exit with an Exception instead of exiting cleanly. Output BAM will still be valid.")
    public double MAX_DISCARD_FRACTION = 0.01;

    @Option(doc="The sample alias to use in the reverted output file.  This will override the existing " +
            "sample alias in the file and is used only if all the read groups in the input file have the " +
            "same sample alias ", shortName=StandardOptionDefinitions.SAMPLE_ALIAS_SHORT_NAME, optional=true)
    public String SAMPLE_ALIAS;

    @Option(doc="The library name to use in the reverted output file.  This will override the existing " +
            "sample alias in the file and is used only if all the read groups in the input file have the " +
            "same sample alias ", shortName=StandardOptionDefinitions.LIBRARY_NAME_SHORT_NAME, optional=true)
    public String LIBRARY_NAME;

    private final static Log log = Log.getInstance(RevertSam.class);

    /** Default main method impl. */
    public static void main(final String[] args) {
        new RevertSam().instanceMainWithExit(args);
    }

    /**
     * Enforce that output ordering is queryname when sanitization is turned on since it requires a queryname sort.
     */
    @Override protected String[] customCommandLineValidation() {
        if (SANITIZE && SORT_ORDER != SortOrder.queryname) {
            return new String[] {"SORT_ORDER must be queryname when sanitization is enabled with SANITIZE=true."};
        }

        return null;
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final boolean sanitizing = SANITIZE;
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileHeader inHeader = in.getFileHeader();

        // If we are going to override SAMPLE_ALIAS or LIBRARY_NAME, make sure all the read
        // groups have the same values.
        final List<SAMReadGroupRecord> rgs = inHeader.getReadGroups();
        if (SAMPLE_ALIAS != null || LIBRARY_NAME != null) {
            boolean allSampleAliasesIdentical = true;
            boolean allLibraryNamesIdentical = true;
            for (int i = 1; i < rgs.size(); i++) {
                if (!rgs.get(0).getSample().equals(rgs.get(i).getSample())) {
                    allSampleAliasesIdentical = false;
                }
                if (!rgs.get(0).getLibrary().equals(rgs.get(i).getLibrary())) {
                    allLibraryNamesIdentical = false;
                }
            }
            if (SAMPLE_ALIAS != null && !allSampleAliasesIdentical) {
                throw new PicardException("Read groups have multiple values for sample.  " +
                        "A value for SAMPLE_ALIAS cannot be supplied." );
            }
            if (LIBRARY_NAME != null && !allLibraryNamesIdentical) {
                throw new PicardException("Read groups have multiple values for library name.  " +
                        "A value for library name cannot be supplied." );
            }
        }

        ////////////////////////////////////////////////////////////////////////////
        // Build the output writer with an appropriate header based on the options
        ////////////////////////////////////////////////////////////////////////////
        final boolean presorted = (inHeader.getSortOrder() == SORT_ORDER) || (SORT_ORDER == SortOrder.queryname && SANITIZE);
        final SAMFileHeader outHeader = new SAMFileHeader();
        for (final SAMReadGroupRecord rg : inHeader.getReadGroups()) {
            if (SAMPLE_ALIAS != null) {
                rg.setSample(SAMPLE_ALIAS);
            }
            if (LIBRARY_NAME != null) {
                rg.setLibrary(LIBRARY_NAME);
            }
            outHeader.addReadGroup(rg);
        }
        outHeader.setSortOrder(SORT_ORDER);
        if (!REMOVE_ALIGNMENT_INFORMATION) {
            outHeader.setSequenceDictionary(inHeader.getSequenceDictionary());
            outHeader.setProgramRecords(inHeader.getProgramRecords());
        }

        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, presorted, OUTPUT);


        ////////////////////////////////////////////////////////////////////////////
        // Build a sorting collection to use if we are sanitizing
        ////////////////////////////////////////////////////////////////////////////
        final SortingCollection<SAMRecord> sorter;
        if (sanitizing) {
            sorter = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(outHeader), new SAMRecordQueryNameComparator(), MAX_RECORDS_IN_RAM);
        }
        else {
            sorter = null;
        }

        final ProgressLogger progress = new ProgressLogger(log, 1000000, "Reverted");
        for (final SAMRecord rec : in) {
            // Weed out non-primary and supplemental read as we don't want duplicates in the reverted file!
            if (rec.isSecondaryOrSupplementary()) continue;

            // Actually to the reverting of the remaining records
            revertSamRecord(rec);

            if (sanitizing) sorter.add(rec);
            else out.addAlignment(rec);
            progress.record(rec);
        }

        ////////////////////////////////////////////////////////////////////////////
        // Now if we're sanitizing, clean up the records and write them to the output
        ////////////////////////////////////////////////////////////////////////////
        if (!sanitizing) {
            out.close();
        }
        else {

            long total = 0, discarded = 0;
            final PeekableIterator<SAMRecord> iterator;
            FastqQualityFormat targetFormat = null;
            final int readGroupCount = inHeader.getReadGroups().size();
            final Map<SAMReadGroupRecord, FastqQualityFormat> readGroupToFormat = new HashMap<SAMReadGroupRecord, FastqQualityFormat>();

            if (readGroupCount > 1) {
                final BufferingIterator buff = new BufferingIterator(sorter.iterator(), ((int)QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE) * readGroupCount * 20);

                for (final SAMReadGroupRecord rg : inHeader.getReadGroups()) {
                    final SamRecordFilter filter = new SamRecordFilter() {
                        public boolean filterOut(final SAMRecord rec) {
                            return !rec.getReadGroup().getId().equals(rg.getId());
                        }
                        public boolean filterOut(final SAMRecord first, final SAMRecord second) {
                            throw new UnsupportedOperationException();
                        }
                    };
                    readGroupToFormat.put(rg, QualityEncodingDetector.detect(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, new FilteringIterator(buff.iterator(), filter)));
                }
                final Set<FastqQualityFormat> formats = new HashSet<FastqQualityFormat>(readGroupToFormat.values());
                for(final SAMReadGroupRecord r : readGroupToFormat.keySet()) {
                    log.info("Detected quality format for " + r.getReadGroupId() + ": " + readGroupToFormat.get(r));
                }
                if (formats.contains(FastqQualityFormat.Solexa)) {
                    log.error("No quality score encoding harmonization implemented for " + FastqQualityFormat.Solexa);
                    return -1;
                }
                //  Standardize to Phred scores if there's more than 1
                targetFormat = (formats.size() == 1) ? formats.iterator().next() : FastqQualityFormat.Standard;
                iterator = buff.iterator();

            }
            else {
                iterator = new PeekableIterator<SAMRecord>(sorter.iterator());
            }



            final ProgressLogger sanitizerProgress = new ProgressLogger(log, 1000000, "Sanitized");

            readNameLoop: while (iterator.hasNext()) {
                final List<SAMRecord> recs = fetchByReadName(iterator);
                total += recs.size();

                // Check that all the reads have bases and qualities of the same length
                for (final SAMRecord rec : recs) {
                    if (rec.getReadBases().length != rec.getBaseQualities().length) {
                        log.debug("Discarding " + recs.size() + " reads with name " + rec.getReadName() + " for mismatching bases and quals length.");
                        discarded += recs.size();
                        continue readNameLoop;
                    }
                }

                // Check that if the first read is marked as unpaired that there is in fact only one read
                if (!recs.get(0).getReadPairedFlag() && recs.size() > 1) {
                    log.debug("Discarding " + recs.size() + " reads with name " + recs.get(0).getReadName() + " because they claim to be unpaired.");
                    discarded += recs.size();
                    continue readNameLoop;
                }

                // Check that if we have paired reads there is exactly one first of pair and one second of pair
                if (recs.get(0).getReadPairedFlag()) {
                    int firsts=0, seconds=0, unpaired=0;
                    for (final SAMRecord rec : recs) {
                        if (!rec.getReadPairedFlag())  ++unpaired;
                        if (rec.getFirstOfPairFlag())  ++firsts;
                        if (rec.getSecondOfPairFlag()) ++seconds;
                    }

                    if (unpaired > 0 || firsts != 1 || seconds != 1) {
                        log.debug("Discarding " + recs.size() + " reads with name " + recs.get(0).getReadName() + " because pairing information in corrupt.");
                        discarded += recs.size();
                        continue readNameLoop;
                    }
                }

                // If we've made it this far spit the records into the output!
                for (final SAMRecord rec : recs) {
                    // Harmonize the quality scores  encoding scores -- either it's all one quality or we're moving them all to Phred
                    final FastqQualityFormat recordFormat = readGroupToFormat.get(rec.getReadGroup());
                    if (!recordFormat.equals(targetFormat)) {
                        final byte quals[] = rec.getBaseQualities();
                        for (int i = 0; i < quals.length; i++) {
                            quals[i] -= SolexaQualityConverter.ILLUMINA_TO_PHRED_SUBTRAHEND;
                        }
                    }
                    out.addAlignment(rec);
                    sanitizerProgress.record(rec);
                }
            }

            out.close();

            final double discardRate = discarded / (double) total;
            final NumberFormat fmt = new DecimalFormat("0.000%");
            log.info("Discarded " + discarded + " out of " + total + " (" + fmt.format(discardRate) + ") reads in order to sanitize output.");

            if (discarded / (double) total > MAX_DISCARD_FRACTION) {
                throw new PicardException("Discarded " + fmt.format(discardRate) + " which is above MAX_DISCARD_FRACTION of " + fmt.format(MAX_DISCARD_FRACTION));
            }
        }

        return 0;
    }

    /**
     * Generates a list by consuming from the iterator in order starting with the first available
     * read and continuing while subsequent reads share the same read name. If there are no reads
     * remaining returns an empty list.
     */
    private List<SAMRecord> fetchByReadName(final PeekableIterator<SAMRecord> iterator) {
        final List<SAMRecord> out = new LinkedList<SAMRecord>();

        if (iterator.hasNext()) {
            final SAMRecord first = iterator.next();
            out.add(first);

            while (iterator.hasNext() && iterator.peek().getReadName().equals(first.getReadName())) {
                out.add(iterator.next());
            }
        }

        return out;
    }

    /**
     * Takes an individual SAMRecord and applies the set of changes/reversions to it that
     * have been requested by program level options.
     */
    public void revertSamRecord(final SAMRecord rec) {
        if (RESTORE_ORIGINAL_QUALITIES) {
            final byte[] oq = rec.getOriginalBaseQualities();
            if (oq != null) {
                rec.setBaseQualities(oq);
                rec.setOriginalBaseQualities(null);
            }
        }

        if (REMOVE_DUPLICATE_INFORMATION) {
            rec.setDuplicateReadFlag(false);
        }

        if (REMOVE_ALIGNMENT_INFORMATION) {
            if (rec.getReadNegativeStrandFlag()) {
                SAMRecordUtil.reverseComplement(rec);
                rec.setReadNegativeStrandFlag(false);
            }

            // Remove all alignment based information about the read itself
            rec.setReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
            rec.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
            rec.setMappingQuality(SAMRecord.NO_MAPPING_QUALITY);

            if (!rec.getReadUnmappedFlag()) {
                rec.setInferredInsertSize(0);
                rec.setNotPrimaryAlignmentFlag(false);
                rec.setProperPairFlag(false);
                rec.setReadUnmappedFlag(true);

            }

            // Then remove any mate flags and info related to alignment
            if (rec.getReadPairedFlag()) {
                rec.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
                rec.setMateNegativeStrandFlag(false);
                rec.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
                rec.setMateUnmappedFlag(true);
            }

            // And then remove any tags that are calculated from the alignment
            for (final String tag : ATTRIBUTE_TO_CLEAR) {
                rec.setAttribute(tag, null);
            }
        }
    }

    /**
     * In order to do quality score encoding detection when sanitizing, we need the ability to read
     * ahead thousands of records into the sorting collection before we start processing it.  This
     * class buffers the first <code>bufferSize</code> records from the underlying iterator, permitting them
     * to be looped through several times, and allows one final iteration that goes past the end of the
     * buffer.  If iterator() is called after passing the end of the buffer, an exception is thrown.
     */
    private static class BufferingIterator implements Iterator<SAMRecord>
    {
        private final ArrayList<SAMRecord> buffer;
        private int bufferIndex = 0;
        private final Iterator<SAMRecord> iterator;

        public BufferingIterator(final Iterator<SAMRecord> it, final int bufferSize) {
            this.iterator = it;
            buffer = new ArrayList<SAMRecord>();
            for (int i=0; i< bufferSize && this.iterator.hasNext(); i++) {
                buffer.add(this.iterator.next());
            }
        }

        /** True if there are more items */
        public boolean hasNext() {
            return bufferIndex < buffer.size() || this.iterator.hasNext();
        }

        /** Returns the next object and advances the iterator.*/
        public SAMRecord next() {
            if (bufferIndex < buffer.size()) {
                return buffer.get(bufferIndex++);
            }
            else {
                buffer.clear();
                return this.iterator.next();
            }
        }

        public PeekableIterator<SAMRecord> iterator() {
            if (bufferIndex >= buffer.size()) {
                throw new IllegalStateException("iterator() called after scrolling past end of buffer: " + buffer.size());
            }
            bufferIndex = 0;
            return new PeekableIterator<SAMRecord>(this);
        }

        /** Unsupported Operation. */
        public void remove() {
            throw new UnsupportedOperationException("Not supported: remove");
        }


    }
}
