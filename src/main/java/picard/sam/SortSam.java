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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = SortSam.USAGE_SUMMARY + SortSam.USAGE_DETAILS,
        oneLineSummary = SortSam.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class SortSam extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Sorts a SAM or BAM file.  ";
    static final String USAGE_DETAILS = "This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property " +
            "of the SAM record. The SortOrder of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled SO.  " +
            "" +
            "<p>For a coordinate sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME) field using the " +
            "reference sequence dictionary (@SQ tag).  Alignments within these subgroups are secondarily sorted using the left-most mapping " +
            "position of the read (POS).  Subsequent to this sorting scheme, alignments are listed arbitrarily.</p>" +
            "" +
            "For queryname-sorted alignments, all alignments are grouped using the queryname field but the alignments are not necessarily sorted within these groups.  " +
            "Reads having the same queryname are derived from the same template.  " +

            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SortSam \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=sorted.bam \\<br />" +
            "      SORT_ORDER=coordinate" +
            "</pre>" +

            "<hr />";
    @Argument(doc = "The BAM or SAM file to sort.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The sorted BAM or SAM output file. ", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc = "Sort order of output file")
    public SAMFileHeader.SortOrder SORT_ORDER;

    private final Log log = Log.getInstance(SortSam.class);

    public static void main(final String[] argv) {
        new SortSam().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        ;
        reader.getFileHeader().setSortOrder(SORT_ORDER);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, OUTPUT);
        writer.setProgressLogger(
                new ProgressLogger(log, (int) 1e7, "Wrote", "records from a sorting collection"));

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");
        for (final SAMRecord rec : reader) {
            writer.addAlignment(rec);
            progress.record(rec);
        }

        log.info("Finished reading inputs, merging and writing to output now.");

        CloserUtil.close(reader);
        writer.close();
        return 0;
    }
}
