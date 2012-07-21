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

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.*;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class SortSam extends CommandLineProgram {
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Sorts the input SAM or BAM.\n" + "" +
            "Input and output formats are determined by file extension.";

    @Option(doc="The BAM or SAM file to sort.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(doc="The sorted BAM or SAM output file. ", shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc="Sort order of output file")
    public SAMFileHeader.SortOrder SORT_ORDER;

    private final Log log = Log.getInstance(SortSam.class);

    public static void main(final String[] argv) {
        new SortSam().instanceMainWithExit(argv);
    }

    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(INPUT));
        reader.getFileHeader().setSortOrder(SORT_ORDER);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), false, OUTPUT);

        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Read");
        for (final SAMRecord rec: reader) {
            writer.addAlignment(rec);
            progress.record(rec);
        }

        log.info("Finished reading inputs, merging and writing to output now.");

        reader.close();
        writer.close();
        return 0;
    }
}
