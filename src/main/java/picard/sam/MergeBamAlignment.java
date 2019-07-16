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

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamPairUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.ReadDataManipulationProgramGroup;
import picard.sam.util.PGTagArgumentCollection;

import java.io.File;
import java.util.*;

/**
 * <h3>Summary</h3>
 * A command-line tool for merging BAM/SAM alignment info from a third-party aligner with the data in an
 * unmapped BAM file, producing a third BAM file that has alignment data (from the aligner) and all the remaining
 * data from the unmapped BAM.
 *
 * Quick note: this is <b>not</b> a tool for taking multiple sam files and creating a bigger file by merging them. For
 * that use-case, see {@link MergeSamFiles}.
 *
 * <h3>Details</h3>
 * Many alignment tools (still!) require fastq format input. The unmapped bam may contain useful information that will
 * be lost in the conversion to fastq (meta-data like sample alias, library, barcodes, etc., and read-level tags.)
 *
 * This tool takes an unaligned bam with meta-data, and the aligned bam produced by calling {@link SamToFastq} and
 * then passing the result to an aligner/mapper. It produces a new SAM file that includes all aligned and unaligned reads
 * and also carries forward additional read attributes from the unmapped BAM (attributes that are otherwise lost in the
 * process of converting to fastq). The resulting file will be valid for use by Picard and GATK tools.
 *
 * The output may be coordinate-sorted, in which case the tags, NM, MD, and UQ will be calculated and populated, or
 * query-name sorted, in which case the tags will not be calculated or populated.
 *
 * <h3>Usage example:</h3>
 * <pre>
 * java -jar picard.jar MergeBamAlignment \\
 *      ALIGNED=aligned.bam \\
 *      UNMAPPED=unmapped.bam \\
 *      O=merge_alignments.bam \\
 *      R=reference_sequence.fasta
 * </pre>
 *
 *
 * <h3>Caveats</h3>
 * This tool has been developing for a while and many arguments have been added to it over the years.
 * You may be particularly interested in the following (partial) list:
 * <ul>
 *     <li>CLIP_ADAPTERS -- Whether to (soft-)clip the ends of the reads that are identified as belonging to adapters</li>
 *     <li>IS_BISULFITE_SEQUENCE -- Whether the sequencing originated from bisulfite sequencing, in which case NM will be
 *     calculated differently</li>
 *     <li>ALIGNER_PROPER_PAIR_FLAGS -- Use if the aligner that was used cannot be trusted to set the "Proper pair" flag
 *     and then the tool will set this flag based on orientation and distance between pairs.</li>
 *     <li>ADD_MATE_CIGAR -- Whether to use this opportunity to add the MC tag to each read.</li>
 *     <li>UNMAP_CONTAMINANT_READS (and MIN_UNCLIPPED_BASES) -- Whether to identify extremely short alignments (with
 *     clipping on both sides) as cross-species contamination and unmap the reads.</li>
 * </ul>
 *
 * @author ktibbett@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = MergeBamAlignment.USAGE_SUMMARY + MergeBamAlignment.USAGE_DETAILS,
        oneLineSummary = MergeBamAlignment.USAGE_SUMMARY,
        programGroup = ReadDataManipulationProgramGroup.class
)
@DocumentedFeature
public class MergeBamAlignment extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Merge alignment data from a SAM or BAM with data in an unmapped BAM file.  ";
    static final String USAGE_DETAILS = "" +
            "<h3>Summary</h3>\n" +
            "A command-line tool for merging BAM/SAM alignment info from a third-party aligner with the data in an " +
            "unmapped BAM file, producing a third BAM file that has alignment data (from the aligner) and all the remaining " +
            "data from the unmapped BAM.\n" +
            "\n" +
            "Quick note: this is <b>not</b> a tool for taking multiple sam files and creating a bigger file by merging them. For " +
            "that use-case, see {@link MergeSamFiles}.\n" +
            "\n" +
            "<h3>Details</h3>\n" +
            "Many alignment tools (still!) require fastq format input. The unmapped bam may contain useful information that will " +
            "be lost in the conversion to fastq (meta-data like sample alias, library, barcodes, etc., and read-level tags.)\n" +
            "\n" +
            "This tool takes an unaligned bam with meta-data, and the aligned bam produced by calling {@link SamToFastq} and " +
            "then passing the result to an aligner/mapper. It produces a new SAM file that includes all aligned and unaligned reads " +
            "and also carries forward additional read attributes from the unmapped BAM (attributes that are otherwise lost in the " +
            "process of converting to fastq). The resulting file will be valid for use by Picard and GATK tools.\n" +
            "\n" +
            "The output may be coordinate-sorted, in which case the tags, NM, MD, and UQ will be calculated and populated, or " +
            "query-name sorted, in which case the tags will not be calculated or populated.\n" +
            "\n" +
            "<h3>Usage example:</h3>\n" +
            "\n" +
            "java -jar picard.jar MergeBamAlignment \\\n" +
            "     ALIGNED=aligned.bam \\\n" +
            "     UNMAPPED=unmapped.bam \\\n" +
            "     O=merge_alignments.bam \\\n" +
            "     R=reference_sequence.fasta\n" +
            "\n" +
            "<h3>Caveats</h3>\n" +
            "This tool has been developing for a while and many arguments have been added to it over the years. " +
            "You may be particularly interested in the following (partial) list:\n" +
            "<ul>\n" +
            "<li>CLIP_ADAPTERS -- Whether to (soft-)clip the ends of the reads that are identified as belonging to adapters</li>\n" +
            "<li>IS_BISULFITE_SEQUENCE -- Whether the sequencing originated from bisulfite sequencing, in which case NM will be " +
            "calculated differently</li>\n" +
            "<li>ALIGNER_PROPER_PAIR_FLAGS -- Use if the aligner that was used cannot be trusted to set the \"Proper pair\" flag " +
            "and then the tool will set this flag based on orientation and distance between pairs.</li>\n" +
            "<li>ADD_MATE_CIGAR -- Whether to use this opportunity to add the MC tag to each read.</li>\n" +
            "<li>UNMAP_CONTAMINANT_READS (and MIN_UNCLIPPED_BASES) -- Whether to identify extremely short alignments (with " +
            "clipping on both sides) as cross-species contamination and unmap the reads.</li>\n" +
            "</ul>\n" +
            "";

    @ArgumentCollection
    public final PGTagArgumentCollection pgTagArgumentCollection = new PGTagArgumentCollection();

    @Argument(shortName = "UNMAPPED",
            doc = "Original SAM or BAM file of unmapped reads, which must be in queryname order.")
    public File UNMAPPED_BAM;

    @Argument(shortName = "ALIGNED",
            doc = "SAM or BAM file(s) with alignment data.",
            mutex = {"READ1_ALIGNED_BAM", "READ2_ALIGNED_BAM"},
            optional = true)
    public List<File> ALIGNED_BAM;

    @Argument(shortName = "R1_ALIGNED",
            doc = "SAM or BAM file(s) with alignment data from the first read of a pair.",
            mutex = {"ALIGNED_BAM"},
            optional = true)
    public List<File> READ1_ALIGNED_BAM;

    @Argument(shortName = "R2_ALIGNED",
            doc = "SAM or BAM file(s) with alignment data from the second read of a pair.",
            mutex = {"ALIGNED_BAM"},
            optional = true)
    public List<File> READ2_ALIGNED_BAM;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Merged SAM or BAM file to write to.")
    public File OUTPUT;

    @Argument(shortName = StandardOptionDefinitions.PROGRAM_RECORD_ID_SHORT_NAME,
            doc = "The program group ID of the aligner (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_RECORD_ID;

    @Argument(shortName = "PG_VERSION",
            doc = "The version of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_VERSION;

    @Argument(shortName = "PG_COMMAND",
            doc = "The command line of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_COMMAND_LINE;

    @Argument(shortName = "PG_NAME",
            doc = "The name of the program group (if not supplied by the aligned file).",
            optional = true)
    public String PROGRAM_GROUP_NAME;

    @Deprecated
    @Argument(doc = "DEPRECATED. This argument is ignored and will be removed.", shortName = "PE", optional=true)
    public Boolean PAIRED_RUN = true;

    @Deprecated
    @Argument(doc = "The expected jump size (required if this is a jumping library). Deprecated. Use EXPECTED_ORIENTATIONS instead",
            shortName = "JUMP",
            mutex = "EXPECTED_ORIENTATIONS",
            optional = true)
    public Integer JUMP_SIZE;

    @Argument(doc = "Whether to clip adapters where identified.")
    public boolean CLIP_ADAPTERS = true;

    @Argument(doc = "Whether the lane is bisulfite sequence (used when calculating the NM tag).")
    public boolean IS_BISULFITE_SEQUENCE = false;

    @Argument(doc = "Whether to output only aligned reads.  ")
    public boolean ALIGNED_READS_ONLY = false;

    @Argument(doc = "The maximum number of insertions or deletions permitted for an alignment to be " +
            "included. Alignments with more than this many insertions or deletions will be ignored. " +
            "Set to -1 to allow any number of insertions or deletions.",
            shortName = "MAX_GAPS")
    public int MAX_INSERTIONS_OR_DELETIONS = 1;

    @Argument(doc = "Reserved alignment attributes (tags starting with X, Y, or Z) that should be " +
            "brought over from the alignment data when merging.", optional = true)
    public List<String> ATTRIBUTES_TO_RETAIN = new ArrayList<>();

    @Argument(doc = "Attributes from the alignment record that should be removed when merging." +
            "  This overrides ATTRIBUTES_TO_RETAIN if they share common tags.", optional = true)
    public List<String> ATTRIBUTES_TO_REMOVE = new ArrayList<>();

    @Argument(shortName="RV", doc="Attributes on negative strand reads that need to be reversed.", optional = true)
    public Set<String> ATTRIBUTES_TO_REVERSE = new TreeSet<>(SAMRecord.TAGS_TO_REVERSE);

    @Argument(shortName="RC", doc="Attributes on negative strand reads that need to be reverse complemented.", optional = true)
    public Set<String> ATTRIBUTES_TO_REVERSE_COMPLEMENT = new TreeSet<>(SAMRecord.TAGS_TO_REVERSE_COMPLEMENT);

    @Argument(shortName = "R1_TRIM",
            doc = "The number of bases trimmed from the beginning of read 1 prior to alignment")
    public int READ1_TRIM = 0;

    @Argument(shortName = "R2_TRIM",
            doc = "The number of bases trimmed from the beginning of read 2 prior to alignment")
    public int READ2_TRIM = 0;

    @Argument(shortName = "ORIENTATIONS",
            doc = "The expected orientation of proper read pairs. Replaces JUMP_SIZE",
            mutex = "JUMP_SIZE",
            optional = true)
    public List<SamPairUtil.PairOrientation> EXPECTED_ORIENTATIONS;

    @Argument(doc = "Use the aligner's idea of what a proper pair is rather than computing in this program.")
    public boolean ALIGNER_PROPER_PAIR_FLAGS = false;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME,
            doc = "The order in which the merged reads should be output.")
    public SortOrder SORT_ORDER = SortOrder.coordinate;

    @Argument(doc = "Strategy for selecting primary alignment when the aligner has provided more than one alignment " +
            "for a pair or fragment, and none are marked as primary, more than one is marked as primary, or the primary " +
            "alignment is filtered out for some reason. For all strategies, ties are resolved arbitrarily.")
    public PrimaryAlignmentStrategy PRIMARY_ALIGNMENT_STRATEGY = PrimaryAlignmentStrategy.BestMapq;

    @Argument(doc = "For paired reads, soft clip the 3' end of each read if necessary so that it does not extend past the 5' end of its mate.")
    public boolean CLIP_OVERLAPPING_READS = true;

    @Argument(doc = "If false, do not write secondary alignments to output.")
    public boolean INCLUDE_SECONDARY_ALIGNMENTS = true;

    @Argument(shortName = "MC", optional = true, doc = "Adds the mate CIGAR tag (MC) if true, does not if false.")
    public Boolean ADD_MATE_CIGAR = true;

    @Argument(shortName = "UNMAP_CONTAM", optional = true, doc = "Detect reads originating from foreign organisms (e.g. bacterial DNA in a non-bacterial sample)," +
            "and unmap + label those reads accordingly.")
    public boolean UNMAP_CONTAMINANT_READS = false;

    @Argument(doc = "If UNMAP_CONTAMINANT_READS is set, require this many unclipped bases or else the read will be marked as contaminant.")
    public int MIN_UNCLIPPED_BASES = 32;

    @Argument(doc = "List of Sequence Records tags that must be equal (if present) in the reference dictionary and in the aligned file. Mismatching tags will cause an error if in this list, and a warning otherwise.")
    public List<String> MATCHING_DICTIONARY_TAGS = SAMSequenceDictionary.DEFAULT_DICTIONARY_EQUAL_TAG;

    @Argument(doc = "How to deal with alignment information in reads that are being unmapped (e.g. due to cross-species contamination.) " +
            "Currently ignored unless UNMAP_CONTAMINANT_READS = true. Note that the DO_NOT_CHANGE strategy will actually reset the cigar and set the mapping quality on unmapped reads since otherwise" +
            "the result will be an invalid record. To force no change use the DO_NOT_CHANGE_INVALID strategy.", optional = true)
    public AbstractAlignmentMerger.UnmappingReadStrategy UNMAPPED_READ_STRATEGY = AbstractAlignmentMerger.UnmappingReadStrategy.DO_NOT_CHANGE;

    @Override
    protected boolean requiresReference() {
        return true;
    }

    /**
     * Mechanism to bridge between command line option and PrimaryAlignmentSelectionStrategy implementation.
     */
    enum PrimaryAlignmentStrategy implements CommandLineParser.ClpEnum{
        BestMapq(BestMapqPrimaryAlignmentSelectionStrategy.class,
                "Expects that multiple alignments will be correlated with HI tag, and prefers the pair of " +
                "alignments with the largest MAPQ, in the absence of a primary selected by the aligner."),
        EarliestFragment(EarliestFragmentPrimaryAlignmentSelectionStrategy.class,
                "Prefers the alignment which maps the earliest base in the read. Note that EarliestFragment " +
                "may not be used for paired reads."),
        BestEndMapq(BestEndMapqPrimaryAlignmentStrategy.class,"Appropriate for cases in which the aligner is not pair-aware, " +
                "and does not output the HI tag. It simply picks the alignment for each end with the highest MAPQ, and makes " +
                "those alignments primary, regardless of whether the two alignments make sense together."),
        MostDistant(MostDistantPrimaryAlignmentSelectionStrategy.class, "Appropriate for a non-pair-aware aligner. Picks " +
                "the alignment pair with the largest insert size. If all alignments would be chimeric, it picks the " +
                "alignments for each end with the best MAPQ. ");

        private final Class<PrimaryAlignmentSelectionStrategy> clazz;

        private final String description;

        public String getHelpDoc() {
            return description;
        }

        PrimaryAlignmentStrategy(final Class<?> clazz, final String description) {
            this.clazz = (Class<PrimaryAlignmentSelectionStrategy>) clazz;
            this.description = description;
        }

        PrimaryAlignmentSelectionStrategy newInstance() {
            try {
                return clazz.newInstance();
            } catch (Exception e) {
                throw new PicardException("Trouble instantiating " + clazz.getName(), e);
            }
        }
    }

    @Override
    protected int doWork() {
        // Check the files are readable/writable
        SAMProgramRecord prod = null;
        if (PROGRAM_RECORD_ID != null) {
            prod = new SAMProgramRecord(PROGRAM_RECORD_ID);
            prod.setProgramVersion(PROGRAM_GROUP_VERSION);
            prod.setCommandLine(PROGRAM_GROUP_COMMAND_LINE);
            prod.setProgramName(PROGRAM_GROUP_NAME);
        }
        // TEMPORARY FIX until internal programs all specify EXPECTED_ORIENTATIONS
        if (JUMP_SIZE != null) {
            EXPECTED_ORIENTATIONS = Collections.singletonList(SamPairUtil.PairOrientation.RF);
        } else if (EXPECTED_ORIENTATIONS == null || EXPECTED_ORIENTATIONS.isEmpty()) {
            EXPECTED_ORIENTATIONS = Collections.singletonList(SamPairUtil.PairOrientation.FR);
        }

        final SamAlignmentMerger merger = new SamAlignmentMerger(UNMAPPED_BAM, OUTPUT,
                referenceSequence.getReferenceFile(), prod, CLIP_ADAPTERS, IS_BISULFITE_SEQUENCE,
                ALIGNED_READS_ONLY, ALIGNED_BAM, MAX_INSERTIONS_OR_DELETIONS,
                ATTRIBUTES_TO_RETAIN, ATTRIBUTES_TO_REMOVE, READ1_TRIM, READ2_TRIM,
                READ1_ALIGNED_BAM, READ2_ALIGNED_BAM, EXPECTED_ORIENTATIONS, SORT_ORDER,
                PRIMARY_ALIGNMENT_STRATEGY.newInstance(), ADD_MATE_CIGAR, UNMAP_CONTAMINANT_READS,
                MIN_UNCLIPPED_BASES, UNMAPPED_READ_STRATEGY, MATCHING_DICTIONARY_TAGS);
        merger.setClipOverlappingReads(CLIP_OVERLAPPING_READS);
        merger.setMaxRecordsInRam(MAX_RECORDS_IN_RAM);
        merger.setKeepAlignerProperPairFlags(ALIGNER_PROPER_PAIR_FLAGS);
        merger.setIncludeSecondaryAlignments(INCLUDE_SECONDARY_ALIGNMENTS);
        merger.setAttributesToReverse(ATTRIBUTES_TO_REVERSE);
        merger.setAttributesToReverseComplement(ATTRIBUTES_TO_REVERSE_COMPLEMENT);
        merger.setAddPGTagToReads(pgTagArgumentCollection.ADD_PG_TAG_TO_READS);
        merger.mergeAlignment(referenceSequence.getReferenceFile());
        merger.close();

        return 0;
    }

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns
     * an array of error messages to be written to the appropriate place.
     */
    protected String[] customCommandLineValidation() {

        if ((PROGRAM_RECORD_ID != null || PROGRAM_GROUP_VERSION != null ||
                PROGRAM_GROUP_COMMAND_LINE != null) &&
                (PROGRAM_RECORD_ID == null || PROGRAM_GROUP_VERSION == null ||
                        PROGRAM_GROUP_COMMAND_LINE == null)) {

            return new String[]{"PROGRAM_RECORD_ID, PROGRAM_GROUP_VERSION, and " +
                    "PROGRAM_GROUP_COMMAND_LINE must all be supplied or none should " +
                    "be included."};
        }

        final boolean r1sExist = READ1_ALIGNED_BAM != null && !READ1_ALIGNED_BAM.isEmpty();
        final boolean r2sExist = READ2_ALIGNED_BAM != null && !READ2_ALIGNED_BAM.isEmpty();
        if ((r1sExist && !r2sExist) || (r2sExist && !r1sExist)) {
            return new String[]{"READ1_ALIGNED_BAM and READ2_ALIGNED_BAM " +
                    "must both be supplied or neither should be included.  For " +
                    "single-end read use ALIGNED_BAM."};
        }
        if (ALIGNED_BAM == null || ALIGNED_BAM.isEmpty() && !(r1sExist && r2sExist)) {
            return new String[]{"Either ALIGNED_BAM or the combination of " +
                    "READ1_ALIGNED_BAM and READ2_ALIGNED_BAM must be supplied."};
        }

        return null;
    }
}
