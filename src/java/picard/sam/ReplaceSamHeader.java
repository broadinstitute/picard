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

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
@CommandLineProgramProperties(
        usage = "Replace the SAMFileHeader in a SAM file with the given header. " +
                "Validation is minimal.  It is up to the user to ensure that all the elements referred to in the SAMRecords " +
                "are present in the new header.  Sort order of the two input files must be the same.",
        usageShort = "Replace the SAMFileHeader in a SAM file with the given header",
        programGroup = SamOrBam.class
)
public class ReplaceSamHeader extends CommandLineProgram {

    @Option(doc = "SAM file from which SAMRecords will be read.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc = "SAM file from which SAMFileHeader will be read.")
    public File HEADER;

    @Option(doc = "SAMFileHeader from HEADER file will be written to this file, followed by SAMRecords from INPUT file",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    public static void main(final String[] argv) {
        new ReplaceSamHeader().instanceMainWithExit(argv);
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(HEADER);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader headerReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(HEADER);
        final SAMFileHeader replacementHeader = headerReader.getFileHeader();
        CloserUtil.close(headerReader);


        if (BamFileIoUtils.isBamFile(INPUT)) {
            blockCopyReheader(replacementHeader);
        } else {
            standardReheader(replacementHeader);
        }

        return 0;
    }

    private void standardReheader(final SAMFileHeader replacementHeader) {
        final SamReader recordReader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).validationStringency(ValidationStringency.SILENT).open(INPUT);
        if (replacementHeader.getSortOrder() != recordReader.getFileHeader().getSortOrder()) {
            throw new PicardException("Sort orders of INPUT (" + recordReader.getFileHeader().getSortOrder().name() +
                    ") and HEADER (" + replacementHeader.getSortOrder().name() + ") do not agree.");
        }
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(replacementHeader, true, OUTPUT);

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(ReplaceSamHeader.class));
        for (final SAMRecord rec : recordReader) {
            rec.setHeader(replacementHeader);
            writer.addAlignment(rec);
            progress.record(rec);
        }
        writer.close();
        CloserUtil.close(recordReader);
    }

    private void blockCopyReheader(final SAMFileHeader replacementHeader) {
        BamFileIoUtils.reheaderBamFile(replacementHeader, INPUT, OUTPUT, CREATE_MD5_FILE, CREATE_INDEX);
    }
}
