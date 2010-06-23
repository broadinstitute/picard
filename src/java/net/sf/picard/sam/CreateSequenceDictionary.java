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
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.*;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Create a SAM/BAM file from a fasta containing reference sequence.  The output SAM file contains a header but no
 * SAMRecords, and the header contains only sequence records.
 */
public class CreateSequenceDictionary extends CommandLineProgram {

    // The following attributes define the command-line arguments
    @Usage
    public String USAGE =
            "Usage: " + getClass().getName() + " [options]\n\n" +
                    "Read fasta or fasta.gz containing reference sequences, and write as a SAM or BAM file with only sequence dictionary.\n";

    @Option(doc = "Input reference fasta or fasta.gz", shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;

    @Option(doc = "Output SAM or BAM file containing only the sequence dictionary",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Option(doc = "Put into AS field of sequence dictionary entry if supplied", optional = true)
    public String GENOME_ASSEMBLY;

    @Option(doc = "Put into UR field of sequence dictionary entry.  If not supplied, input reference file is used",
            optional = true)
    public String URI;

    @Option(doc = "Put into SP field of sequence dictionary entry", optional = true)
    public String SPECIES;

    @Option(doc = "Make sequence name the first word from the > line in the fasta file.  " +
            "By default the entire contents of the > line is used, excluding leading and trailing whitespace.")
    public boolean TRUNCATE_NAMES_AT_WHITESPACE = true;

    @Option(doc = "Stop after writing this many sequences.  For testing.")
    public int NUM_SEQUENCES = Integer.MAX_VALUE;

    private final MessageDigest md5;

    public CreateSequenceDictionary() {
        try {
            md5 = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException e) {
            throw new PicardException("MD5 algorithm not found", e);
        }
    }

    public static void main(final String[] argv) {
        System.exit(new CreateSequenceDictionary().instanceMain(argv));
    }

    /**
     * Use reference filename to create URI to go into header if URI was not passed on cmd line.
     */
    protected String[] customCommandLineValidation() {
        if (URI == null) {
            URI = "file:" + REFERENCE.getAbsolutePath();
        }
        return null;
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    protected int doWork() {
        final SAMSequenceDictionary sequences = makeSequenceDictionary(REFERENCE);
        final SAMFileHeader samHeader = new SAMFileHeader();
        samHeader.setSequenceDictionary(sequences);
        final SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(samHeader, false, OUTPUT);
        samWriter.close();
        return 0;
    }


    /**
     * Read all the sequences from the given reference file, and convert into SAMSequenceRecords
     * @param referenceFile fasta or fasta.gz
     * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
     */
    SAMSequenceDictionary makeSequenceDictionary(final File referenceFile) {
        final ReferenceSequenceFile refSeqFile =
                ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, TRUNCATE_NAMES_AT_WHITESPACE);
        ReferenceSequence refSeq;
        final List<SAMSequenceRecord> ret = new ArrayList<SAMSequenceRecord>();
        final Set<String> sequenceNames = new HashSet<String>();
        for (int numSequences = 0; numSequences < NUM_SEQUENCES && (refSeq = refSeqFile.nextSequence()) != null; ++numSequences) {
            if (sequenceNames.contains(refSeq.getName())) {
                throw new PicardException("Sequence name appears more than once in reference: " + refSeq.getName());
            }
            sequenceNames.add(refSeq.getName());
            ret.add(makeSequenceRecord(refSeq));
        }
        return new SAMSequenceDictionary(ret);
    }

    /**
     * Create one SAMSequenceRecord from a single fasta sequence
     */
    private SAMSequenceRecord makeSequenceRecord(final ReferenceSequence refSeq) {
        final SAMSequenceRecord ret = new SAMSequenceRecord(refSeq.getName(), refSeq.length());

        // Compute MD5 of upcased bases
        final byte[] bases = refSeq.getBases();
        for (int i = 0; i < bases.length; ++i) {
                bases[i] = StringUtil.toUpperCase(bases[i]);
            }

        ret.setAttribute(SAMSequenceRecord.MD5_TAG, md5Hash(bases));
        if (GENOME_ASSEMBLY != null) {
            ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, GENOME_ASSEMBLY);
        }
        ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
        if (SPECIES != null) {
                ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, SPECIES);
            }
        return ret;
    }

    private String md5Hash(final byte[] bytes) {
        md5.reset();
        md5.update(bytes);
        String s = new BigInteger(1, md5.digest()).toString(16);
        if (s.length() != 32) {
            final String zeros = "00000000000000000000000000000000";
            s = zeros.substring(0, 32 - s.length()) + s;
        }
        return s;
    }
}
