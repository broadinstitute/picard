/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;

/**
 * @author alecw@broadinstitute.org
 */
public class CleanSam extends CommandLineProgram {
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Read SAM and perform various fix-ups.  " +
            "Currently, the only fix-ups are 1: to soft-clip an alignment that hangs off the end of its reference sequence; " +
            "and 2: to set MAPQ to 0 if a read is unmapped.";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM to be cleaned.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write cleaned SAM.")
    public File OUTPUT;

    public static void main(final String[] argv) {
        new CleanSam().instanceMainWithExit(argv);
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
        IOUtil.assertFileIsWritable(OUTPUT);
        final ValidationStringency originalStringency = SAMFileReader.getDefaultValidationStringency();
        if (VALIDATION_STRINGENCY == ValidationStringency.STRICT) {
            SAMFileReader.setDefaultValidationStringency(ValidationStringency.LENIENT);
        }
        try {
            final SAMFileReader reader = new SAMFileReader(INPUT);
            final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
            final CloseableIterator<SAMRecord> it = reader.iterator();
            final ProgressLogger progress = new ProgressLogger(Log.getInstance(CleanSam.class));

            // If the read (or its mate) maps off the end of the alignment, clip it
            while(it.hasNext()) {
                final SAMRecord rec = it.next();

                // If the read (or its mate) maps off the end of the alignment, clip it
                AbstractAlignmentMerger.createNewCigarsIfMapsOffEndOfReference(rec);

                // check the read's mapping quality
                if (rec.getReadUnmappedFlag() && 0 != rec.getMappingQuality()) {
                    rec.setMappingQuality(0);
                }

                writer.addAlignment(rec);
                progress.record(rec);
            }

            writer.close();
            it.close();
        }
        finally {
            SAMFileReader.setDefaultValidationStringency(originalStringency);
        }
        return 0;
    }
}
