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
import htsjdk.samtools.seekablestream.SeekablePathStream;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

/**
 * <p/>
 * Splits the input queryname sorted or query-grouped SAM/BAM file and writes it into
 * multiple BAM files, each approximately the same size in bytes. This will retain the sort order
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
            "while maintaining the sort order. This can be used to split a large unmapped BAM in order to parallelize alignment." +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SplitSamBySize \\<br />" +
            "     I=paired_unmapped_input.bam \\<br />" +
            "     OUTPUT=out_dir \\ <br />" +
            "     TOTAL_READS_IN_INPUT=800000000 \\ <br />" +
            "     SPLIT_TO_N_READS=48000000" +
            "</pre>" +
            "<hr />";
    @Argument(doc = "Input SAM/BAM file to split", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = "N_FILES", doc = "Split to N files.")
    public int SPLIT_TO_N_FILES;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Directory in which to output the split BAM files.")
    public File OUTPUT;

    @Argument(shortName = "OUT_PREFIX", doc = "Output files will be named <OUT_PREFIX>_N.bam, where N enumerates the output file.")
    public String OUT_PREFIX = "shard";

    private final Log log = Log.getInstance(SplitSamBySize.class);

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertDirectoryIsWritable(OUTPUT);
        SeekablePathStream inStream;
        try {
            inStream = new SeekablePathStream(INPUT.toPath());
        } catch (IOException e) {
            e.printStackTrace();
            return 1;
        }
        final SamInputResource inResource = SamInputResource.of(inStream);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(inResource);
        final SAMFileHeader header = reader.getFileHeader();

        if (header.getSortOrder() == SAMFileHeader.SortOrder.coordinate) {
            log.warn("Splitting a coordinate sorted bam may result in invalid bams " +
                    "that do not always contain each read's mate in the same bam.");
        }
        if (!header.getVersion().equals(SAMFileHeader.CURRENT_VERSION)) {
            log.warn(String.format("Input file's version is %s, but the current SAM format version is %s. Outputs will be written " +
                    "with current version.", header.getVersion(), SAMFileHeader.CURRENT_VERSION));
        }

        final SAMFileWriterFactory factory = new SAMFileWriterFactory();

        final long totalCompressedBytes = INPUT.length();
        final int inputBytesPerFile = (int) Math.ceil(totalCompressedBytes / (double) SPLIT_TO_N_FILES);
        final SAMFileWriter[] writers = generateWriters(factory, header, SPLIT_TO_N_FILES);

        int writerIndex = 0;
        String lastReadName = "";
        ArrayList<SAMRecord> blockRecords = new ArrayList();
        long bufferedBlockEnd = 0;
        long bufferedBlockStart = 0;
        final SAMRecordIterator iterator = reader.iterator();
        final ProgressLogger progress = new ProgressLogger(log);
        // Iterate through the input file, buffer each block's worth of reads so we can split within a block
        while (iterator.hasNext()) {
            final SAMRecord currentRecord = iterator.next();
            final long currentBlockEnd;
            try {
                currentBlockEnd = inStream.position();
            } catch (IOException e) {
                e.printStackTrace();
                return 1;
            }
            if (!iterator.hasNext()) {
                // This is the last record of the bam, add it to the previous block so it gets written out
                blockRecords.add(currentRecord);
            }

            // We finished the block, write it out
            if (currentBlockEnd != bufferedBlockEnd || !iterator.hasNext()) {
                int idxInBlock = 0;
                final long bufferedBlockSize = bufferedBlockEnd - bufferedBlockStart;
                for (final SAMRecord record : blockRecords) {
                    // Estimate the current position in bytes by subdividing the buffered block by the number of records in it
                    final long currentPosition = bufferedBlockStart + ((idxInBlock * bufferedBlockSize) / blockRecords.size());
                    if (currentPosition >= inputBytesPerFile * (writerIndex + 1) && !lastReadName.equals(record.getReadName())) {
                        writerIndex++;
                    }
                    writers[writerIndex].addAlignment(record);
                    progress.record(record);
                    lastReadName = record.getReadName();
                    idxInBlock++;
                }
                blockRecords.clear();
                bufferedBlockStart = bufferedBlockEnd;
                bufferedBlockEnd = currentBlockEnd;
            }
            blockRecords.add(currentRecord);
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
        final SAMFileWriter[] writers = new SAMFileWriter[splitToNFiles];
        for (int i = 0; i < splitToNFiles; i++) {
            writers[i] = factory.makeSAMOrBAMWriter(header, true, new File(OUTPUT, fileNameFormatter.format(fileIndex++)));
        }
        return writers;
    }

    protected String[] customCommandLineValidation() {
        if (SPLIT_TO_N_FILES < 1) {
            return new String[]{
                    String.format("Cannot set SPLIT_TO_N_FILES to a number less than 1, found %s.", SPLIT_TO_N_FILES)
            };
        }

        return null;
    }
}
