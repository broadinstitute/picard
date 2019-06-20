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
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;

import java.io.File;

/**
 * Converts a BAM file to human-readable SAM output or vice versa
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = "Convert a BAM file to a SAM file, or SAM to BAM.\n" +
                "Input and output formats are determined by file extension.",
        oneLineSummary = "Convert a BAM file to a SAM file, or a SAM to a BAM",
        programGroup = ReadDataManipulationProgramGroup.class)
@DocumentedFeature
public class SamFormatConverter extends CommandLineProgram {

    // The following attributes define the command-line arguments
    @Argument(doc = "The BAM or SAM file to parse.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(doc = "The BAM or SAM output file. ", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    protected int doWork() {
        convert(INPUT, OUTPUT, REFERENCE_SEQUENCE, CREATE_INDEX);
        return 0;
    }

    /**
     * Convert a file from one of sam/bam/cram format to another based on the extension of output.
     *
     * @param input             input file in one of sam/bam/cram format
     * @param output            output to write converted file to, the conversion is based on the extension of this filename
     * @param referenceSequence the reference sequence to use, necessary when reading/writing cram
     * @param createIndex       whether or not an index should be written alongside the output file
     */
    public static void convert(final File input, final File output, final File referenceSequence, final Boolean createIndex) {
        IOUtil.assertFileIsReadable(input);
        IOUtil.assertFileIsWritable(output);
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(referenceSequence).open(input);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(reader.getFileHeader(), true, output, referenceSequence);
        if (createIndex && writer.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate) {
            throw new PicardException("Can't CREATE_INDEX unless sort order is coordinate");
        }

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(SamFormatConverter.class));
        for (final SAMRecord rec : reader) {
            writer.addAlignment(rec);
            progress.record(rec);
        }
        CloserUtil.close(reader);
        writer.close();
    }
}
