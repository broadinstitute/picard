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

package picard.reference;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 * Simple command line program that allows sub-sequences represented by an interval
 * list to be extracted from a reference sequence file.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = ExtractSequences.USAGE_SUMMARY + ExtractSequences.USAGE_DETAILS,
        oneLineSummary = ExtractSequences.USAGE_SUMMARY,
        programGroup = ReferenceProgramGroup.class
)
@DocumentedFeature
public class ExtractSequences extends CommandLineProgram {
    static final String USAGE_SUMMARY ="Subsets intervals from a reference sequence to a new FASTA file.";
    static final String USAGE_DETAILS ="This tool takes a list of intervals, reads the corresponding subsquences from a reference " +
            "FASTA file and writes them to a new FASTA file as separate records. Note that the reference FASTA file must be " +
            "accompanied by an index file and the interval list must be provided in Picard list format. The names provided for the " +
            "intervals will be used to name the corresponding records in the output file." +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar ExtractSequences \\<br />" +
            "      INTERVAL_LIST=regions_of_interest.interval_list \\<br />" +
            "      R=reference.fasta \\<br />" +
            "      O=extracted_IL_sequences.fasta" +
            "</pre>" +
            "<hr />";
    @Argument(doc="Interval list describing intervals to be extracted from the reference sequence.")
    public File INTERVAL_LIST;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output FASTA file.")
    public File OUTPUT;

    @Argument(doc="Maximum line length for sequence data.")
    public int LINE_LENGTH = 80;

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INTERVAL_LIST);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervals = IntervalList.fromFile(INTERVAL_LIST);
        final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
        SequenceUtil.assertSequenceDictionariesEqual(intervals.getHeader().getSequenceDictionary(), ref.getSequenceDictionary());

        final BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);

        for (final Interval interval : intervals) {
            final ReferenceSequence seq = ref.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
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
