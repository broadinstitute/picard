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
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.*;

import java.io.File;
import java.util.Iterator;

/**
 * Converts a BAM file to human-readable SAM output or vice versa
 *
 * @author ktibbett@broadinstitute.org
 */
public class SamFormatConverter extends CommandLineProgram {

    private static final String PROGRAM_VERSION = "1.0";

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Convert a BAM file to a SAM file, or BAM to SAM.\n" + "" +
            "Input and output formats are determined by file extension.";

    @Option(doc="The BAM or SAM file to parse.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME) public File INPUT;
    @Option(doc="The BAM or SAM output file. ", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME) public File OUTPUT;

    public static void main(final String[] argv) {
        new SamFormatConverter().instanceMainWithExit(argv);
    }


    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(INPUT));
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);

        if  (CREATE_INDEX && writer.getFileHeader().getSortOrder() != SAMFileHeader.SortOrder.coordinate){
            throw new PicardException("Can't CREATE_INDEX unless sort order is coordinate");
        }

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(SamFormatConverter.class));
        for (final SAMRecord rec : reader) {
            writer.addAlignment(rec);
            progress.record(rec);
        }
        reader.close();
        writer.close();
        return 0;
    }
}
