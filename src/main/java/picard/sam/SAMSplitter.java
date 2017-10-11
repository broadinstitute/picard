/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * <p/>
 * Splits the input queryname sorted or querygrouped SAM/BAM file and writes it into
 * multiple BAM files of approximately equal size. This will retain the sort order
 * within each output BAM and if the BAMs are concatenated in order (output files are named
 * numerically) the order of the reads will match the original BAM. As this tool is
 * intended for unmapped input files, all secondary and supplementary reads will be removed.
 */
@CommandLineProgramProperties(
        summary = SAMSplitter.USAGE_SUMMARY + SAMSplitter.USAGE_DETAILS,
        oneLineSummary = SAMSplitter.USAGE_SUMMARY,
        programGroup = SamOrBam.class)
@DocumentedFeature
public class SAMSplitter extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Splits a SAM or BAM file to multiple BAMs.";
    static final String USAGE_DETAILS = "This tool splits the input SAM/BAM file into multiple BAM files " +
            "This can be used to split a large unmapped BAM in order to parallelize alignment."+
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SAMSplitter \\<br />" +
            "     I=paired_unmapped_input.bam \\<br />" +
            "     OUTPUT_DIR=out_dir \\ <br />" +
            "     TOTAL_READS_IN_INPUT=800000000" +
            "     SPLIT_TO_N_TEMPLATES=24000000" +
            "</pre>" +
            "<hr />";
    @Argument(doc = "Input SAM/BAM file to split", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = "N", doc = "Split to have approximately N templates per output file.")
    public int SPLIT_TO_N_TEMPLATES;

    @Argument(shortName = "TOTAL_READS", doc = "Total number of reads in the input file.")
    public int TOTAL_READS_IN_INPUT;

    @Argument(shortName = "ODIR", doc = "Directory in which to output the split BAM files.")
    public File OUTPUT_DIR;

    @Argument(shortName = "PAIR", doc = "If true the input bam contains paired reads and " +
            "each template is assumed to contain two reads (rather than one for unpaired inputs).")
    public boolean PAIRED_INPUT = true;

    private final Log log = Log.getInstance(SAMSplitter.class);

    public static void main(final String[] argv) {
        System.exit(new SAMSplitter().instanceMain(argv));
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileHeader header = reader.getFileHeader();
        if (header.getSortOrder() == SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("Cannot split coordinate sorted bam. Only queryname sorted (or querygrouped) " +
                    "input will ensure that templates are not split into multiple files.");
        }

        final Map<String, SAMRecord> firstSeenMates = new HashMap<>();
        final SAMFileWriterFactory factory = new SAMFileWriterFactory();
        double unitConverter = PAIRED_INPUT ? 2.0 : 1.0;

        int splitToNFiles = (int) Math.round((TOTAL_READS_IN_INPUT / unitConverter) / SPLIT_TO_N_TEMPLATES);
        if (splitToNFiles < 1) {
            throw new PicardException("Splitting " + TOTAL_READS_IN_INPUT + " reads into " + SPLIT_TO_N_TEMPLATES +
                    " templates per file would result in fewer than 1 output file.");
        }

        int templatesPerFile = (int) Math.ceil((TOTAL_READS_IN_INPUT / (double) splitToNFiles) / unitConverter);
        SAMFileWriter[] writers = generateWriters(factory, header, splitToNFiles);

        int templateIndex = 0;
        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord currentRecord : reader) {
            if (currentRecord.isSecondaryOrSupplementary())
                continue;

            int writerIndex = templateIndex / templatesPerFile;
            if (PAIRED_INPUT) {
                final String currentReadName = currentRecord.getReadName();
                final SAMRecord firstRecord = firstSeenMates.remove(currentReadName);
                if (firstRecord == null) {
                    firstSeenMates.put(currentReadName, currentRecord);
                } else {
                    final SAMRecord read1 =
                            currentRecord.getFirstOfPairFlag() ? currentRecord : firstRecord;
                    final SAMRecord read2 =
                            currentRecord.getFirstOfPairFlag() ? firstRecord : currentRecord;
                    writers[writerIndex].addAlignment(read1);
                    writers[writerIndex].addAlignment(read2);
                    templateIndex++;
                }
            } else {
                writers[writerIndex].addAlignment(currentRecord);
                templateIndex++;
            }

            progress.record(currentRecord);
        }

        for (final SAMFileWriter w : writers) {
            w.close();
        }

        if (!firstSeenMates.isEmpty()) {
            SAMUtils.processValidationError(new SAMValidationError(SAMValidationError.Type.MATE_NOT_FOUND,
                    "Found " + firstSeenMates.size() + " unpaired mates", null), VALIDATION_STRINGENCY);
        }
        return 0;
    }

    private SAMFileWriter[] generateWriters(final SAMFileWriterFactory factory, final SAMFileHeader header, int splitToNFiles) {
        int digits = String.valueOf(splitToNFiles).length();
        final DecimalFormat fileNameFormatter = new DecimalFormat(String.format("%0" + digits + "d", 0));
        int fileIndex = 1;
        SAMFileWriter[] writers = new SAMFileWriter[splitToNFiles];
        for (int i = 0; i < splitToNFiles; i++) {
            writers[i] = factory.makeSAMOrBAMWriter(header, true, new File(OUTPUT_DIR.getAbsolutePath() + "/" +
                    fileNameFormatter.format(fileIndex++) + BamFileIoUtils.BAM_FILE_EXTENSION));
        }
        return writers;
    }

    protected String[] customCommandLineValidation() {
        if (SPLIT_TO_N_TEMPLATES < 1) {
            return new String[]{
                    "Cannot set SPLIT_TO_N_TEMPLATES to a number less than 1."
            };
        }

        if (TOTAL_READS_IN_INPUT < 1) {
            return new String[]{
                    "Cannot set TOTAL_READS_IN_INPUT to a number less than 1."
            };
        }

        return null;
    }
}