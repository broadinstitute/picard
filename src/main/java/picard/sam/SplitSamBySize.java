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
import java.util.HashMap;
import java.util.Map;

/**
 * <p/>
 * Splits the input queryname sorted or query-grouped SAM/BAM file and writes it into
 * multiple BAM files, each with an approximately equal number of reads. This will retain the sort order
 * within each output BAM and if the BAMs are concatenated in order (output files are named
 * numerically) the order of the reads will match the original BAM.
 */
@CommandLineProgramProperties(
        summary = SplitSamBySize.USAGE_SUMMARY + SplitSamBySize.USAGE_DETAILS,
        oneLineSummary = SplitSamBySize.USAGE_SUMMARY,
        programGroup = SamOrBam.class)
@DocumentedFeature
public class SplitSamBySize extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Splits a SAM or BAM file to multiple BAMs.";
    static final String USAGE_DETAILS = "This tool splits the input query-grouped SAM/BAM file into multiple BAM files " +
            "while maintaining the sort order. This can be used to split a large unmapped BAM in order to parallelize alignment."+
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SplitSamBySize \\<br />" +
            "     I=paired_unmapped_input.bam \\<br />" +
            "     OUTPUT_DIR=out_dir \\ <br />" +
            "     TOTAL_READS_IN_INPUT=800000000 \\ <br />" +
            "     SPLIT_TO_N_READS=48000000" +
            "</pre>" +
            "<hr />";
    @Argument(doc = "Input SAM/BAM file to split", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = "N_READS", doc = "Split to have approximately N reads per output file.", optional = true,
            mutex = {"SPLIT_TO_N_FILES"})
    public int SPLIT_TO_N_READS;

    @Argument(shortName = "N_FILES", doc = "Split to N files.", optional = true, mutex = {"SPLIT_TO_N_READS"})
    public int SPLIT_TO_N_FILES;

    @Argument(shortName = "TOTAL_READS", doc = "Total number of reads in the input file.")
    public int TOTAL_READS_IN_INPUT;

    @Argument(shortName = "O", doc = "Directory in which to output the split BAM files.")
    public File OUTPUT_DIR;

    @Argument(shortName = "OUT_PREFIX", doc = "Prefix to name output files.")
    public String OUT_PREFIX = "shard";

    private final Log log = Log.getInstance(SplitSamBySize.class);

    public static void main(final String[] argv) {
        System.exit(new SplitSamBySize().instanceMain(argv));
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileHeader header = reader.getFileHeader();
        if (header.getSortOrder() == SAMFileHeader.SortOrder.coordinate) {
            log.warn("Splitting a coordinate sorted bam may result in invalid bams " +
                    "that do not always contain each read's mate in the same bam.");
        }

        final SAMFileWriterFactory factory = new SAMFileWriterFactory();

        int splitToNFiles = SPLIT_TO_N_FILES != 0 ? SPLIT_TO_N_FILES : (int) Math.ceil(TOTAL_READS_IN_INPUT / SPLIT_TO_N_READS);
        int readsPerFile = (int) Math.ceil(TOTAL_READS_IN_INPUT / (double) splitToNFiles);
        SAMFileWriter[] writers = generateWriters(factory, header, splitToNFiles);

        int recordsWritten = 0;
        int writerIndex = 0;
        String lastReadName = "";
        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord currentRecord : reader) {
            if (recordsWritten >= readsPerFile && !lastReadName.equals(currentRecord.getReadName())) {
                writerIndex++;
                recordsWritten = 0;
            }
            writers[writerIndex].addAlignment(currentRecord);
            lastReadName = currentRecord.getReadName();
            recordsWritten++;
            progress.record(currentRecord);
        }

        for (final SAMFileWriter w : writers) {
            w.close();
        }

        return 0;
    }

    private SAMFileWriter[] generateWriters(final SAMFileWriterFactory factory, final SAMFileHeader header, final int splitToNFiles) {
        final int digits = String.valueOf(splitToNFiles).length();
        final DecimalFormat fileNameFormatter = new DecimalFormat(OUT_PREFIX + "_" + String.format("%0" + digits + "d", 0) + BamFileIoUtils.BAM_FILE_EXTENSION);
        int fileIndex = 1;
        SAMFileWriter[] writers = new SAMFileWriter[splitToNFiles];
        for (int i = 0; i < splitToNFiles; i++) {
            writers[i] = factory.makeSAMOrBAMWriter(header, true, new File(OUTPUT_DIR, fileNameFormatter.format(fileIndex++)));
        }
        return writers;
    }

    protected String[] customCommandLineValidation() {
        if (TOTAL_READS_IN_INPUT < 1) {
            return new String[]{
                    "Cannot set TOTAL_READS_IN_INPUT to a number less than 1."
            };
        }

        return null;
    }
}
