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
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;
import java.nio.file.Files;
import java.util.function.Function;

/**
 * <p/>
 * Splits the input queryname sorted or query-grouped SAM/BAM file and writes it into
 * multiple BAM files, each with an approximately equal number of reads. This will retain the sort order
 * within each output BAM and if the BAMs are concatenated in order (output files are named
 * numerically) the order of the reads will match the original BAM. It will traverse the bam twice unless
 * TOTAL_READS_IN_INPUT is provided.
 */
@CommandLineProgramProperties(
        summary = SplitSamByNumberOfReads.USAGE_SUMMARY + SplitSamByNumberOfReads.USAGE_DETAILS,
        oneLineSummary = SplitSamByNumberOfReads.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class SplitSamByNumberOfReads extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Splits a SAM/BAM/CRAM file to multiple files.";
    static final String USAGE_DETAILS = "This tool splits the input query-grouped SAM/BAM/CRAM file into multiple files " +
            "while maintaining the sort order. This can be used to split a large unmapped input in order to parallelize alignment. " +
            "It will traverse the input twice unless TOTAL_READS_IN_INPUT is provided." +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SplitSamByNumberOfReads \\<br />" +
            "     I=paired_unmapped_input.bam \\<br />" +
            "     OUTPUT=out_dir \\ <br />" +
            "     TOTAL_READS_IN_INPUT=800000000 \\ <br />" +
            "     SPLIT_TO_N_READS=48000000" +
            "</pre>" +
            "<hr />";
    @Argument(doc = "Input SAM/BAM/CRAM file to split", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = "N_READS", doc = "Split to have approximately N reads per output file. The actual number of " +
            "reads per output file will vary by no more than the number of output files * (the maximum " +
            "number of reads with the same queryname - 1).", mutex = {"SPLIT_TO_N_FILES"})
    public int SPLIT_TO_N_READS;

    @Argument(shortName = "N_FILES", doc = "Split to N files.", mutex = {"SPLIT_TO_N_READS"})
    public int SPLIT_TO_N_FILES;

    @Argument(shortName = "TOTAL_READS", doc = "Total number of reads in the input file. If this is not provided, the " +
            "input will be read twice, the first time to get a count of the total reads.", optional = true)
    public long TOTAL_READS_IN_INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Directory in which to output the split BAM files.")
    public File OUTPUT;

    @Argument(shortName = "OUT_PREFIX", doc = "Output files will be named <OUT_PREFIX>_N.EXT, where N enumerates the output file and EXT is the " +
            "same as that of the input.")
    public String OUT_PREFIX = "shard";

    private final Log log = Log.getInstance(SplitSamByNumberOfReads.class);

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        if (TOTAL_READS_IN_INPUT == 0 && !Files.isRegularFile(INPUT.toPath())) {
            log.error(String.format("INPUT is not a regular file: %s. " +
                    "If TOTAL_READS_IN_INPUT is not supplied, INPUT cannot be a stream.", INPUT));
            return 1;
        }
        IOUtil.assertDirectoryIsWritable(OUTPUT);
        final SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
        final SamReader reader = readerFactory.referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final String extension = reader.type().fileExtension();

        final SAMFileHeader header = reader.getFileHeader();

        if (header.getSortOrder() == SAMFileHeader.SortOrder.coordinate) {
            log.warn("Splitting a coordinate sorted bam may result in invalid bams " +
                    "that do not always contain each read's mate in the same bam.");
        }
        if (!header.getVersion().equals(SAMFileHeader.CURRENT_VERSION)) {
            log.warn(String.format("Input file's version is %s, but the current SAM format version is %s. Outputs will be written " +
                    "with current version.", header.getVersion(), SAMFileHeader.CURRENT_VERSION));
        }

        final ProgressLogger firstPassProgress = new ProgressLogger(log, 1000000, "Counted");
        if (TOTAL_READS_IN_INPUT == 0) {
            final SamReader firstPassReader = readerFactory.referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
            log.info("First pass traversal to count number of reads is beginning. If number of reads " +
                    "is known, use TOTAL_READS_IN_INPUT to skip first traversal.");
            for (SAMRecord rec : firstPassReader) {
                firstPassProgress.record(rec);
            }
            CloserUtil.close(firstPassReader);
            log.info(String.format("First pass traversal to count number of reads ended, found %d total reads.", firstPassProgress.getCount()));
        }
        final long totalReads = TOTAL_READS_IN_INPUT == 0 ? firstPassProgress.getCount() : TOTAL_READS_IN_INPUT;

        final SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();

        final int splitToNFiles = SPLIT_TO_N_FILES != 0 ? SPLIT_TO_N_FILES : (int) Math.ceil(totalReads / (double) SPLIT_TO_N_READS);
        final int readsPerFile = (int) Math.ceil(totalReads / (double) splitToNFiles);

        int readsWritten = 0;
        int fileIndex = 1;

        Function<Integer, SAMFileWriter> createWriter = (index) ->
                writerFactory.makeSAMOrBAMWriter(header, true,
                        new File(OUTPUT, String.format("%s_%04d.%s", OUT_PREFIX,index,extension)));

        SAMFileWriter currentWriter = createWriter.apply(fileIndex++);
        String lastReadName = "";
        final ProgressLogger progress = new ProgressLogger(log);
        for (SAMRecord currentRecord : reader) {
            if (readsWritten >= readsPerFile && !lastReadName.equals(currentRecord.getReadName())) {
                currentWriter.close();

                currentWriter = createWriter.apply(fileIndex++);
                readsWritten = 0;
            }
            currentWriter.addAlignment(currentRecord);
            lastReadName = currentRecord.getReadName();
            readsWritten++;
            progress.record(currentRecord);
        }
        currentWriter.close();
        CloserUtil.close(reader);

        if (progress.getCount() != totalReads) {
            log.warn(String.format("The totalReads (%d) provided does not match the reads found in the " +
                            "input file (%d). Files may not be split evenly or number of files may not " +
                            "match what was requested. There were %d files generated each with around %d " +
                            "reads except the last file which contained %d reads.",
                    totalReads, progress.getCount(), fileIndex - 1, readsPerFile, readsWritten)
            );
        }

        return 0;
    }

    protected String[] customCommandLineValidation() {
        if (TOTAL_READS_IN_INPUT < 0) {
            return new String[]{
                    String.format("Cannot set TOTAL_READS_IN_INPUT to a number less than 1, found %d.", TOTAL_READS_IN_INPUT)
            };
        }

        if (SPLIT_TO_N_FILES <= 1 && SPLIT_TO_N_READS <= 1) {
            return new String[]{
                    String.format("One of SPLIT_TO_N_FILES or SPLIT_TO_N_READS must be greater than 0. " +
                            "Found SPLIT_TO_N_FILES is %d and SPLIT_TO_N_READS is %d.", SPLIT_TO_N_FILES, SPLIT_TO_N_READS)
            };
        }

        return null;
    }
}
