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

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceDictionaryCodec;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import picard.cmdline.StandardOptionDefinitions;

import java.io.*;
import java.math.BigInteger;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Create a SAM/BAM file from a fasta containing reference sequence. The output SAM file contains a header but no
 * SAMRecords, and the header contains only sequence records.
 */
@DocumentedFeature
@CommandLineProgramProperties(
        summary = CreateSequenceDictionary.USAGE_SUMMARY + CreateSequenceDictionary.USAGE_DETAILS,
        oneLineSummary = CreateSequenceDictionary.USAGE_SUMMARY,
        programGroup = ReferenceProgramGroup.class
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

    @Argument(doc = "Output SAM file containing only the sequence dictionary. By default it will use the base name of the input reference with the .dict extension",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true)
    public File OUTPUT;

    @Argument(shortName = "AS", doc = "Put into AS field of sequence dictionary entry if supplied", optional = true)
    public String GENOME_ASSEMBLY;

    @Argument(shortName = "UR", doc = "Put into UR field of sequence dictionary entry.  If not supplied, input reference file is used",
            optional = true)
    public String URI;

    @Argument(shortName = "SP", doc = "Put into SP field of sequence dictionary entry", optional = true)
    public String SPECIES;

    @Argument(doc = "Make sequence name the first word from the > line in the fasta file.  " +
            "By default the entire contents of the > line is used, excluding leading and trailing whitespace.")
    public boolean TRUNCATE_NAMES_AT_WHITESPACE = true;

    @Argument(doc = "Stop after writing this many sequences.  For testing.")
    public int NUM_SEQUENCES = Integer.MAX_VALUE;

    @Argument(shortName = "AN", doc = "Optional file containing the alternative names for the contigs. "
    		+ "Tools may use this information to consider different contig notations as identical (e.g: 'chr1' and '1'). "
    		+ "The alternative names will be put into the appropriate @AN annotation for each contig. "
    		+ "No header. "
    		+ "First column is the original name, the second column is an alternative name. "
            + "One contig may have more than one alternative name." ,
            optional=true)
    public File ALT_NAMES = null;

    private final MessageDigest md5;

    /**
     * 'AN' attribute in the dictionary
     * TODO: replace "AN" with a constant : see https://github.com/samtools/htsjdk/pull/956/files
     */
    private static final String AN_ATTRIBUTE = "AN";

    public CreateSequenceDictionary() {
        try {
            md5 = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException e) {
            throw new PicardException("MD5 algorithm not found", e);
        }
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
            URI = "file:" + referenceSequence.getReferenceFile().getAbsolutePath();
        }
        if (OUTPUT == null) {
            OUTPUT = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(referenceSequence.getReferenceFile());
            logger.info("Output dictionary will be written in ", OUTPUT);
        }
        return super.customCommandLineValidation();
    }

    // return a custom argument collection because this tool uses the argument name
    // "REFERENCE" instead of "REFERENCE_SEQUENCE"
    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new CreateSeqDictReferenceArgumentCollection();
    }

    public static class CreateSeqDictReferenceArgumentCollection implements ReferenceArgumentCollection {
        @Argument(doc = "Input reference fasta or fasta.gz", shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME)
        public File REFERENCE;

        @Override
        public File getReferenceFile() {
            return REFERENCE;
        };
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

        // map for aliases mapping a contig to its aliases
        final Map<String, Set<String>> aliasesByContig = loadContigAliasesMap();

        // SortingCollection is used to check uniqueness of sequence names
        final SortingCollection<String> sequenceNames = makeSortingCollection();
        try (BufferedWriter writer = makeWriter()) {
            final ReferenceSequenceFile refSeqFile = ReferenceSequenceFileFactory.
                    getReferenceSequenceFile(REFERENCE_SEQUENCE, TRUNCATE_NAMES_AT_WHITESPACE);
            SAMSequenceDictionaryCodec samDictCodec = new SAMSequenceDictionaryCodec(writer);

            samDictCodec.encodeHeaderLine(false);
            // read reference sequence one by one and write its metadata
            for (ReferenceSequence refSeq = refSeqFile.nextSequence(); refSeq != null; refSeq = refSeqFile.nextSequence()) {
                final SAMSequenceRecord samSequenceRecord = makeSequenceRecord(refSeq);
                // retrieve aliases, if any
                final Set<String> aliases = aliasesByContig.get(samSequenceRecord.getSequenceName());
                if (aliases != null) {
                    // "Alternative names is a comma separated list of alternative names"
                    samSequenceRecord.setAttribute(AN_ATTRIBUTE, String.join(",", aliases));
                }
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

        if (!iterator.hasNext()) return 0;

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

        ret.setAttribute(SAMSequenceRecord.MD5_TAG, SequenceUtil.calculateMD5String(bases));

        if (GENOME_ASSEMBLY != null) {
            ret.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, GENOME_ASSEMBLY);
        }
        ret.setAttribute(SAMSequenceRecord.URI_TAG, URI);
        if (SPECIES != null) {
                ret.setAttribute(SAMSequenceRecord.SPECIES_TAG, SPECIES);
            }
        return ret;
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

    /**
     * Load the file ALT_NAMES containing the alternative contig names
     * @author Pierre Lindenbaum
     * @return a <code>Map&lt;src_contig,Set&lt;new_names&gt;&gt;</code>. Never null. May be empty if ALT_NAMES is null.
     * @throws PicardException if there is any error in the file ALT_NAMES
     */
    private Map<String, Set<String>> loadContigAliasesMap() throws PicardException {
        // return an empty map if no mapping file was provided
        if (this.ALT_NAMES == null) {
            return Collections.emptyMap();
        }
        // the map returned by the function
        final Map<String, Set<String>> aliasesByContig = new HashMap<>();
        try {
            for (final String line : IOUtil.slurpLines(this.ALT_NAMES)) {
                if (StringUtil.isBlank(line)) {
                        continue;
                    }
                final int tab = line.indexOf('\t');
                if (tab == -1) {
                    throw new IOException("tabulation missing in " + line);
                }
                final String contigName = line.substring(0, tab);
                final String altName = line.substring(tab + 1);
                // check for empty values
                if (StringUtil.isBlank(contigName)) {
                    throw new IOException("empty contig in " + line);
                }
                if (StringUtil.isBlank(altName)) {
                    throw new IOException("empty alternative name in " + line);
                }
                if (altName.equals(contigName)) {
                    continue;
                }
                try {
                    SAMSequenceRecord.validateSequenceName(altName);
                } catch (final SAMException exception) {
                    throw new IOException("Illegal alternative reference sequence name in " + line, exception);
                }
                // check alias not previously defined as contig
                if (aliasesByContig.containsKey(altName)) {
                    throw new IOException("alternate name " + altName +
                            " previously defined as a contig in " + line);
                }
                // check contig not previously defined as alias
                if (aliasesByContig.keySet().stream().
                        // not an error if defined twice for same contig
                        filter(K -> !K.equals(contigName)). 
                        anyMatch(K -> aliasesByContig.get(K).contains(contigName))) {
                            throw new IOException("contig " + contigName +
                                " previously defined as an alternate name in " + line);
                }
                // add alias
                if (!aliasesByContig.containsKey(contigName)) {
                    aliasesByContig.put(contigName, new HashSet<>());
                }
                aliasesByContig.get(contigName).add(altName);
            }
            return aliasesByContig;  
        } catch (final IOException e) {
            throw new PicardException("Can't read alias file " + ALT_NAMES, e);
        }
   
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
