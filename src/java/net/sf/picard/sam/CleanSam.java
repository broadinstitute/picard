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
package net.sf.picard.sam;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.CigarUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;

import java.io.File;
import java.util.List;

/**
 * @author alecw@broadinstitute.org
 */
public class CleanSam extends CommandLineProgram {
    @Usage
    public String USAGE = getStandardUsagePreamble() + "Read SAM and perform various fix-ups.  " +
            "Currently, the only fix-up it to soft-clip an alignment that hangs off the end of its reference sequence.";

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
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileReader.ValidationStringency originalStringency = SAMFileReader.getDefaultValidationStringency();
        if (VALIDATION_STRINGENCY == SAMFileReader.ValidationStringency.STRICT) {
            SAMFileReader.setDefaultValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
        }
        try {
            final SAMFileReader reader = new SAMFileReader(INPUT);
            final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);
            final CloseableIterator<SAMRecord> it = reader.iterator();
            final ProgressLogger progress = new ProgressLogger(Log.getInstance(CleanSam.class));

            // If the read maps off the end of the alignment, clip it
            while(it.hasNext()) {
                final SAMRecord rec = it.next();
                if (!rec.getReadUnmappedFlag()) {
                    final SAMSequenceRecord refseq = rec.getHeader().getSequence(rec.getReferenceIndex());
                    if (rec.getAlignmentEnd() > refseq.getSequenceLength()) {
                        // 1-based index of first base in read to clip.
                        final int clipFrom = refseq.getSequenceLength() - rec.getAlignmentStart() + 1;
                        final List<CigarElement> newCigarElements  = CigarUtil.softClipEndOfRead(clipFrom, rec.getCigar().getCigarElements());
                        rec.setCigar(new Cigar(newCigarElements));
                    }
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
