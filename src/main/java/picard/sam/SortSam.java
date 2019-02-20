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
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;


/**
 * Sorts a SAM or BAM file.
 *
 * <h3> Summary </h3>
 * This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property of the SAM
 * record. The SortOrder of a SAM/BAM file is found in the SAM file header tag labeled SO.
 * <p> For a coordinate sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME)
 * field using the reference sequence dictionary tag labeled SQ.
 * Alignments within these subgroups are secondarily sorted using the left-most mapping position of the read (POS).
 * Subsequent to this sorting scheme, alignments are listed arbitrarily.</p> For queryname-sorted alignments, the tool
 * orders records deterministically by queryname field followed by record strand orientation flag, primary record flag, and secondary
 * alignment flag. (See {@link htsjdk.samtools.SAMRecordQueryNameComparator#compare(SAMRecord, SAMRecord)}} for details)
 *
 * <h3> Inputs</h3>
 * <ul>
 *     <li> Input BAM or SAM file to sort </li>
 *     <li> Sorted BAM or SAM output file </li>
 *     <li> Sort order of output file </li>
 * </ul>
 *
 * <h3>Usage example:</h3>
 * <pre>
 *     java -jar picard.jar SortSam \
 *     INPUT=input.bam \
 *     OUTPUT=sorted.bam \
 *     SORT_ORDER=coordinate
 * </pre>
 *
 *
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = SortSam.USAGE_DETAILS,
        oneLineSummary = SortSam.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class SortSam extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Sorts a SAM or BAM file";
    static final String USAGE_DETAILS = "This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or some other property " +
            "of the SAM record. The SortOrder of a SAM/BAM file is found in the SAM file header tag @HD in the field labeled SO.  " +
            "<p> For a coordinate sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME) field using the " +
            "reference sequence dictionary (@SQ tag).  Alignments within these subgroups are secondarily sorted using the left-most mapping " +
            "position of the read (POS).  Subsequent to this sorting scheme, alignments are listed arbitrarily.</p>" +
            "<p> For queryname-sorted alignments, the tool orders records deterministically by queryname field followed by " +
            "record strand orientation flag, primary record flag, and secondary alignment flag.</p>" +
            "<hr /><hr />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SortSam \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=sorted.bam \\<br />" +
            "      SORT_ORDER=coordinate" +
            "</pre>";
    @Argument(doc = "Input BAM or SAM file to sort.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "Sorted BAM or SAM output file.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    // note that SortOrder here is a local enum, not the SamFileHeader version.
    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc = "Sort order of output file. ")
    public SortOrder SORT_ORDER;


    private final Log log = Log.getInstance(SortSam.class);

    /** a SortOrder class intended to expose the various options available as inputs to SortSam
     *
     * In particular this enables to add a description and also to not expose "unsorted" and "unknown"
     * as they are not appropriate values to sort a file into.
     *
     */
    private enum SortOrder implements CommandLineParser.ClpEnum {
        queryname("Sorts according to the readname. This will place read-pairs and other derived reads (secondary and " +
                "supplementary) adjacent to each other. Note that the readnames are compared lexicographically, even though " +
                "they may include numbers. In paired reads, Read1 sorts before Read2."),
        coordinate("Sorts primarily according to the SEQ and POS fields of the record. The sequence will sorted according to " +
                "the order in the sequence dictionary, taken from from the header of the file. Within each reference sequence, the " +
                "reads are sorted by the position. Unmapped reads whose mates are mapped will be placed near their mates. " +
                "Unmapped read-pairs are placed after all the mapped reads and their mates."),
        duplicate("Sorts the reads so that duplicates reads are adjacent. Required that the mate-cigar (MC) tag is present. " +
                "The resulting will be sorted by library, unclipped 5-prime position, orientation, and mate's unclipped " +
                "5-prime position.")
        ;

        private String description;
        SortOrder(String description) {
            this.description=description;
        }

        public SAMFileHeader.SortOrder getSortOrder(){
            return SAMFileHeader.SortOrder.valueOf(this.name());
        }

        @Override
        public String getHelpDoc() {
            return description;
        }
    }


    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        ;
        reader.getFileHeader().setSortOrder(SORT_ORDER.getSortOrder());
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
