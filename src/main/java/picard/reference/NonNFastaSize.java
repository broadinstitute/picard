/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalListReferenceSequenceMask;
import htsjdk.samtools.util.ReferenceSequenceMask;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.samtools.util.WholeGenomeReferenceSequenceMask;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fasta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

/**
 * A tool to count the number of non-N bases in a fasta file
 */
@CommandLineProgramProperties(
        summary = NonNFastaSize.USAGE_SUMMARY + NonNFastaSize.USAGE_DETAILS,
        oneLineSummary = NonNFastaSize.USAGE_SUMMARY,
        programGroup = Fasta.class
)
public class NonNFastaSize extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Counts the number of non-N bases in a fasta file.";

    static final String USAGE_DETAILS = "This tool takes any FASTA-formatted file and counts the number of non-N bases in it." +
            "Note that it requires that the fasta file have associated index (.fai) and dictionary (.dict) files.<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar NonNFastaSize \\<br />" +
            "      I=input_sequence.fasta \\<br />" +
            "      O=count.txt" +
            "</pre>" +
            "<hr />"
            ;
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input FASTA file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output file in which to record the count.")
    public File OUTPUT;

    @Argument(shortName = "INTERVALS", doc = "An interval list file that contains the locations of the positions to assess.  If not provided, the entire reference will be used", optional = true)
    public File INTERVALS = null;

    public static void main(final String[] args) {
        new NonNFastaSize().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        // set up the reference and a mask so that we only count the positions requested by the user
        final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(INPUT);
        final ReferenceSequenceMask referenceSequenceMask;
        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            final IntervalList intervalList = IntervalList.fromFile(INTERVALS);
            referenceSequenceMask = new IntervalListReferenceSequenceMask(intervalList);
        } else {
            final SAMFileHeader header = new SAMFileHeader();
            header.setSequenceDictionary(ref.getSequenceDictionary());
            referenceSequenceMask = new WholeGenomeReferenceSequenceMask(header);
        }

        long nonNbases = 0L;

        for (final SAMSequenceRecord rec : ref.getSequenceDictionary().getSequences()) {
            // pull out the contig and set up the bases
            final ReferenceSequence sequence = ref.getSequence(rec.getSequenceName());
            final byte[] bases = sequence.getBases();
            StringUtil.toUpperCase(bases);

            for (int i = 0; i < bases.length; i++) {
                // only investigate this position if it's within our mask
                if (referenceSequenceMask.get(sequence.getContigIndex(), i+1)) {
                    nonNbases += bases[i] == SequenceUtil.N ? 0 : 1;
                }
            }
        }

        try {
            final BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
            out.write(nonNbases + "\n");
            out.close();
        }
        catch (IOException ioe) {
            throw new PicardException("Error writing to file " + OUTPUT.getAbsolutePath(), ioe);
        }

        return 0;
    }
}
