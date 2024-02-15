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
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Md5CalculatingOutputStream;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.ReferenceProgramGroup;
import picard.nio.PicardBucketUtils;
import picard.nio.PicardHtsPath;
import picard.util.SequenceDictionaryUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
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
            "The reference sequence can be gzipped (both .fasta and .fasta.gz are supported)." +
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
    public PicardHtsPath OUTPUT;

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
            + "One contig may have more than one alternative name.",
            optional = true)
    public File ALT_NAMES = null;

    private final MessageDigest md5;
    /**
     * Regular expression defined in the sam spec. Any alternative contig should match this regular expression
     * TODO: replace the pattern with a constant : see https://github.com/samtools/htsjdk/pull/956/files
     */
    private static final Pattern ALTERNATIVE_CONTIG_NAME_PATTERN = Pattern.compile("[0-9A-Za-z][0-9A-Za-z\\*\\+@\\|\\-]*");

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
     *
     * @param referenceFile fasta or fasta.gz
     * @return SAMSequenceRecords containing info from the fasta, plus from cmd-line arguments.
     */
    public SAMSequenceDictionary makeSequenceDictionary(final File referenceFile) {
        final Iterable<SAMSequenceRecord> samSequenceRecordsIterable = getSamSequenceRecordsIterable();

        final List<SAMSequenceRecord> ret = new ArrayList<>();
        final Set<String> sequenceNames = new HashSet<>();
        for (SAMSequenceRecord rec : samSequenceRecordsIterable) {

            if (sequenceNames.contains(rec.getSequenceName())) {
                throw new PicardException("Sequence name appears more than once in reference: " + rec.getSequenceName());
            }
            sequenceNames.add(rec.getSequenceName());
            ret.add(rec);
        }
        return new SAMSequenceDictionary(ret);
    }

    private Iterable<SAMSequenceRecord> getSamSequenceRecordsIterable() {
        return () -> {
            final SequenceDictionaryUtils.SamSequenceRecordsIterator iterator =
                    new SequenceDictionaryUtils.SamSequenceRecordsIterator(referenceSequence.getReferencePath(),
                            TRUNCATE_NAMES_AT_WHITESPACE);
            iterator.setGenomeAssembly(GENOME_ASSEMBLY);
            iterator.setSpecies(SPECIES);
            iterator.setUri(URI);
            return iterator;
        };
    }

    /**
     * Use reference filename to create URI to go into header if URI was not passed on cmd line.
     */
    protected String[] customCommandLineValidation() {
        if (URI == null) {
            URI = referenceSequence.getHtsPath().getURIString();
        }
        if (OUTPUT == null) {
            final Path outputPath = ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(referenceSequence.getReferencePath());
            OUTPUT = new PicardHtsPath(outputPath.toString());
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
        public PicardHtsPath REFERENCE;

        @Override
        public PicardHtsPath getHtsPath() {
            return REFERENCE;
        }

        @Override
        public File getReferenceFile() {
            return ReferenceArgumentCollection.getFileSafe(REFERENCE, logger);
        }
    }

    /**
     * Do the work after command line has been parsed.
     * RuntimeException may be thrown by this method, and are reported appropriately.
     *
     * @return program exit status.
     */
    protected int doWork() {
        int sequencesWritten = 0;

        if (Files.exists(OUTPUT.toPath())) {
            throw new PicardException(OUTPUT.getURIString() +
                    " already exists.  Delete this file and try again, or specify a different output file.");
        }

        // We can check for writability provided the file is in a local filesystem and not in gcloud
        if (OUTPUT.getScheme().equals(PicardBucketUtils.FILE_SCHEME)){
            IOUtil.assertFileIsWritable(OUTPUT.toPath().toFile());
        }

        // map for aliases mapping a contig to its aliases
        final Map<String, Set<String>> aliasesByContig = loadContigAliasesMap();

        try (BufferedWriter writer = makeWriter()) {
            final Iterable<SAMSequenceRecord> samSequenceRecordIterable = getSamSequenceRecordsIterable();
            SAMSequenceDictionaryCodec samDictCodec = new SAMSequenceDictionaryCodec(writer);

            samDictCodec.encodeHeaderLine(false);
            // read reference sequence one by one and write its metadata
            for (SAMSequenceRecord samSequenceRecord : samSequenceRecordIterable) {
                // retrieve aliases, if any
                final Set<String> aliases = aliasesByContig.get(samSequenceRecord.getSequenceName());
                if (aliases != null) {
                    // "Alternative names is a comma separated list of alternative names"
                    samSequenceRecord.setAttribute(AN_ATTRIBUTE, String.join(",", aliases));
                }
                samDictCodec.encodeSequenceRecord(samSequenceRecord);

                if (++sequencesWritten >= NUM_SEQUENCES) {
                    break;
                }
            }
        } catch (IOException e) {
            throw new PicardException("Can't write to or close output file " + OUTPUT.getURIString(), e);
        } catch (IllegalArgumentException e) {
            // in case of an unexpected error delete the file so that there isn't a
            // truncated result which might be valid yet wrong.
            if (Files.exists(OUTPUT.toPath())){
                try {
                    Files.delete(OUTPUT.toPath());
                } catch (IOException e2) {
                    throw new PicardException("Unknown problem encountered, and failed to delete the incomplete output. " + e2.getMessage(), e);
                }
            }
            throw new PicardException("Unknown problem. Partial dictionary file was deleted.", e);
        }

        return 0;
    }

    private BufferedWriter makeWriter() throws IOException {
        return new BufferedWriter(
                new AsciiWriter(this.CREATE_MD5_FILE ?
                        new Md5CalculatingOutputStream(
                                Files.newOutputStream(OUTPUT.toPath(), StandardOpenOption.CREATE_NEW),
                                new PicardHtsPath(OUTPUT.getURIString() + ".md5").toPath()
                        )
                        : Files.newOutputStream(OUTPUT.toPath(), StandardOpenOption.CREATE_NEW)
                )
        );
    }


    /**
     * Load the file ALT_NAMES containing the alternative contig names
     *
     * @return a <code>Map&lt;src_contig,Set&lt;new_names&gt;&gt;</code>. Never null. May be empty if ALT_NAMES is null.
     * @throws PicardException if there is any error in the file ALT_NAMES
     * @author Pierre Lindenbaum
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
                if (aliasesByContig.keySet().stream()
                        // not an error if defined twice for same contig
                        .filter(K -> !K.equals(contigName))
                        .anyMatch(K -> aliasesByContig.get(K).contains(contigName))) {
                    throw new IOException("contig  " + contigName +
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
}

