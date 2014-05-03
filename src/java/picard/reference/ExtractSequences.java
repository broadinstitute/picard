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

package net.sf.picard.reference;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.util.SequenceUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 * Simple command line program that allows sub-sequences represented by an interval
 * list to be extracted from a reference sequence file.
 *
 * @author Tim Fennell
 */
public class ExtractSequences extends CommandLineProgram {
    @Usage public final String USAGE = "Extracts one or more intervals described in an interval_list file " +
            "from a given reference sequence and writes them out in FASTA format. Requires a fasta index " +
            "file to be present.";

    @Option(doc="Interval list describing intervals to be extracted from the reference sequence.")
    public File INTERVAL_LIST;

    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference sequence file.")
    public File REFERENCE_SEQUENCE;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output fasta file.")
    public File OUTPUT;

    @Option(doc="Maximum line length for sequence data.")
    public int LINE_LENGTH = 80;

    public static void main(final String[] args) {
        new ExtractSequences().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INTERVAL_LIST);
        IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IoUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervals = IntervalList.fromFile(INTERVAL_LIST);
        final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        SequenceUtil.assertSequenceDictionariesEqual(intervals.getHeader().getSequenceDictionary(), ref.getSequenceDictionary());

        final BufferedWriter out = IoUtil.openFileForBufferedWriting(OUTPUT);

        for (final Interval interval : intervals) {
            final ReferenceSequence seq = ref.getSubsequenceAt(interval.getSequence(), interval.getStart(), interval.getEnd());
            final byte[] bases = seq.getBases();
            if (interval.isNegativeStrand()) SequenceUtil.reverseComplement(bases);

            try {
                out.write(">");
                out.write(interval.getName());
                out.write("\n");

                for (int i=0; i<bases.length; ++i) {
                    if (i > 0 && i % LINE_LENGTH == 0) out.write("\n");
                    out.write(bases[i]);
                }

                out.write("\n");
            }
            catch (IOException ioe) {
                throw new PicardException("Error writing to file " + OUTPUT.getAbsolutePath(), ioe);

            }
        }

        CloserUtil.close(out);

        return 0;
    }
}
