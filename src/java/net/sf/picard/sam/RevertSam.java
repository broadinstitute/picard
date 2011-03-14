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
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

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
        add("NM");
        add("UQ");
        add("PG");
        add("MD");
        add("MQ");
    }};

    @Option(doc="The sample alias to use in the reverted output file.  This will override the existing " +
            "sample alias in the file and is used only if all the read groups in the input file have the " +
            "same sample alias ", shortName=StandardOptionDefinitions.SAMPLE_ALIAS_SHORT_NAME, optional=true)
    public String SAMPLE_ALIAS;

    @Option(doc="The library name to use in the reverted output file.  This will override the existing " +
            "sample alias in the file and is used only if all the read groups in the input file have the " +
            "same sample alias ", shortName=StandardOptionDefinitions.LIBRARY_NAME_SHORT_NAME, optional=true)
    public String LIBRARY_NAME;


    /** Default main method impl. */
    public static void main(final String[] args) {
        System.exit(new RevertSam().instanceMain(args));
    }

    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        final SAMFileReader in = new SAMFileReader(INPUT, true);
        final SAMFileHeader inHeader = in.getFileHeader();

        // If we are going to override SAMPLE_ALIAS or LIBRARY_NAME, make sure all the read
        // groups have the same values.
        List<SAMReadGroupRecord> rgs = inHeader.getReadGroups();
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


        // Build the output writer with an appropriate header based on the options
        final boolean presorted = inHeader.getSortOrder() == SORT_ORDER;
        final SAMFileHeader outHeader = new SAMFileHeader();
        for (SAMReadGroupRecord rg : inHeader.getReadGroups()) {
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

        for (final SAMRecord rec : in) {
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

            out.addAlignment(rec);
        }

        out.close();

        return 0;
    }
}
