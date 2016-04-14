/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;

/**
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = FixNmAndUqTags.USAGE_SUMMARY + FixNmAndUqTags.USAGE_DETAILS,
        usageShort = FixNmAndUqTags.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class FixNmAndUqTags extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Fixes the UQ and NM tags in a SAM file.  ";
    static final String USAGE_DETAILS = "This tool takes in a SAM or BAM file (sorted by coordinate) and calculates the NM and UQ tags by comparing with the reference."+
            "<br />" +
                    "This may be needed when MergeBamAlignment was run with SORT_ORDER different from 'coordinate' and thus could not fix\n"+
                    "these tags then.<br />"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar FixNmAndUqTags \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=sorted.bam \\<br />"+
            "</pre>" +
            "<hr />";
    @Option(doc = "The BAM or SAM file to fix.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "The fixed BAM or SAM output file. ", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Option(doc = "If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = false;

    @Option(doc = "Whether the file contains bisulfite sequence (used when calculating the NM tag).")
    public boolean IS_BISULFITE_SEQUENCE = false;

    @Override
    protected String[] customCommandLineValidation() {
        if (REFERENCE_SEQUENCE == null)
            return new String[]
                    {"Must have a non-null REFERENCE_SEQUENCE"};
        return super.customCommandLineValidation();
    }

    private final Log log = Log.getInstance(FixNmAndUqTags.class);

    public static void main(final String[] argv) {
        new FixNmAndUqTags().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        if( !ASSUME_SORTED && reader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new SAMException("Input must be coordinate-sorted, or ASSUME_SORTED must be true.");
        }
        if( ASSUME_SORTED && reader.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            log.warn("File indicates sort-order " + reader.getFileHeader().getSortOrder() + " but ASSUME_SORTED is true, so continuing.");
        }
        reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
        writer.setProgressLogger(
                new ProgressLogger(log, (int) 1e7, "Wrote", "records"));

        final ReferenceSequenceFileWalker refSeq = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);

        for (final SAMRecord rec : reader) {
            if (!rec.getReadUnmappedFlag()) {
                AbstractAlignmentMerger.fixNMandUQ(rec, refSeq, IS_BISULFITE_SEQUENCE);
            }
            writer.addAlignment(rec);
        }

        CloserUtil.close(reader);
        writer.close();
        return 0;
    }
}
