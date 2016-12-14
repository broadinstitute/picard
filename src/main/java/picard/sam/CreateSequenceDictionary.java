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

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Md5CalculatingOutputStream;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Fasta;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

/**
 * Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no
 * SAMRecords, and the header contains only sequence records.
 */
@CommandLineProgramProperties(
        usage = CreateSequenceDictionary.USAGE_SUMMARY + CreateSequenceDictionary.USAGE_DETAILS,
        usageShort = CreateSequenceDictionary.USAGE_SUMMARY,
        programGroup = Fasta.class
)
public class CreateSequenceDictionary extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Creates a sequence dictionary for a reference sequence.  ";
    static final String USAGE_DETAILS = "This tool creates a sequence dictionary file (with \".dict\" extension) from a reference " +
            "sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a " +
            "header but no SAMRecords, and the header contains only sequence records." +
            "<br /><br />" +
            "The reference sequence can be gzipped (both .fasta and .fasta.gz are supported)."  +
            "" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CreateSequenceDictionary \\ <br />" +
            "      R=reference.fasta \\ <br />" +
            "      O=reference.dict" +
            "" +
            "</pre>" +
            "<hr />";
    // The following attributes define the command-line arguments

    private static final Log logger = Log.getInstance(CreateSequenceDictionary.class);

    @Option(doc = "Input reference fasta or fasta.gz", shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME)
    public File REFERENCE;

    @Option(doc = "Output SAM file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
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
     * Read all the sequences from the given reference file, and convert into SAMSequenceRecords
     * @param referenceFile fasta or fasta.gz
     * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
     */
    @Deprecated
    public SAMSequenceDictionary makeSequenceDictionary(final File referenceFile) {
        final ReferenceSequenceFile refSeqFile =
                ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, TRUNCATE_NAMES_AT_WHITESPACE);
        ReferenceSequence refSeq;
        final List<SAMSequenceRecord> ret = new ArrayList<>();
        final Set<String> sequenceNames = new HashSet<>();
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
     * Use reference filename to create URI to go into header if URI was not passed on cmd line.
     */
    protected String[] customCommandLineValidation() {
        if (URI == null) {
            URI = "file:" + REFERENCE.getAbsolutePath();
        }
        if (OUTPUT == null) {
            // TODO: use the htsjdk method implemented in https://github.com/samtools/htsjdk/pull/774
            OUTPUT = getDefaultDictionaryForReferenceSequence(REFERENCE);
            logger.info("Output dictionary will be written in ", OUTPUT);
        }
        return null;
    }

    // TODO: this method will be in htsjdk (https://github.com/samtools/htsjdk/pull/774)
    @VisibleForTesting
    static File getDefaultDictionaryForReferenceSequence(final File fastaFile) {
        final String name = fastaFile.getName();
        final String extension = ReferenceSequenceFileFactory.FASTA_EXTENSIONS.stream().filter(name::endsWith).findFirst()
                .orElseGet(() -> {throw new IllegalArgumentException("File is not a supported reference file type: " + fastaFile.getAbsolutePath());});
        final int extensionIndex = name.length() - extension.length();
        return new File(fastaFile.getParentFile(), name.substring(0, extensionIndex) + IOUtil.DICT_FILE_EXTENSION);
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    protected int doWork() {
        if (OUTPUT.exists()) {
            throw new PicardException(OUTPUT.getAbsolutePath() +
                    " already exists.  Delete this file and try again, or specify a different output file.");
        }

        // SortingCollection is used to check uniqueness of sequence names
        final SortingCollection<String> sequenceNames = makeSortingCollection();
        try (BufferedWriter writer = makeWriter()) {
            final ReferenceSequenceFile refSeqFile = ReferenceSequenceFileFactory.
                    getReferenceSequenceFile(REFERENCE, TRUNCATE_NAMES_AT_WHITESPACE);
            SAMSequenceDictionaryCodec samDictCodec = new SAMSequenceDictionaryCodec(writer);

            samDictCodec.encodeHeaderLine(false);
            // read reference sequence one by one and write its metadata
            for (ReferenceSequence refSeq = refSeqFile.nextSequence(); refSeq != null; refSeq = refSeqFile.nextSequence()) {
                final SAMSequenceRecord samSequenceRecord = makeSequenceRecord(refSeq);
                samDictCodec.encodeSequenceRecord(samSequenceRecord);
                sequenceNames.add(refSeq.getName());
            }
        } catch (FileNotFoundException e) {
            throw new PicardException("File " + OUTPUT.getAbsolutePath() + " not found");
        } catch (IOException e) {
            throw new PicardException("Can't write to or close output file " + OUTPUT.getAbsolutePath());
        }

        // check uniqueness of sequences names
        final CloseableIterator<String> iterator = sequenceNames.iterator();

        if(!iterator.hasNext()) return 0;

        String current = iterator.next();
        while (iterator.hasNext()) {
            final String next = iterator.next();
            if (current.equals(next)) {
                OUTPUT.delete();
                throw new PicardException("Sequence name " + current +
                        " appears more than once in reference file");
            }
            current = next;
        }
        return 0;
    }

    private BufferedWriter makeWriter() throws FileNotFoundException {
        return new BufferedWriter(
                new AsciiWriter(this.CREATE_MD5_FILE ?
                        new Md5CalculatingOutputStream(
                                new FileOutputStream(OUTPUT, false),
                                new File(OUTPUT.getAbsolutePath() + ".md5")
                        )
                        : new FileOutputStream(OUTPUT)
                )
        );
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

    private SortingCollection<String> makeSortingCollection() {
        final String name = getClass().getSimpleName();
        final File tmpDir = IOUtil.createTempDir(name, null);
        tmpDir.deleteOnExit();
        // 256 byte for one name, and 1/10 part of all memory for this, rough estimate
        long maxNamesInRam = Runtime.getRuntime().maxMemory() / 256 / 10;
        return SortingCollection.newInstance(
                String.class,
                new StringCodec(),
                String::compareTo,
                (int) Math.min(maxNamesInRam, Integer.MAX_VALUE),
                tmpDir
        );
    }

    private static class StringCodec implements SortingCollection.Codec<String> {
        private DataInputStream dis;
        private DataOutputStream dos;

        public StringCodec clone() {
            return new StringCodec();
        }

        public void setOutputStream(final OutputStream os) {
            dos = new DataOutputStream(os);
        }

        public void setInputStream(final InputStream is) {
            dis = new DataInputStream(is);
        }

        public void encode(final String str) {
            try {
                dos.writeUTF(str);
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

        public String decode() {
            try {
                return dis.readUTF();
            } catch (EOFException e) {
                return null;
            } catch (IOException e) {
                throw new PicardException("Exception reading sequence name from temporary file.", e);
            }
        }
    }
}
