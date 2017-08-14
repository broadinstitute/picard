package picard.vcf;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Tool for lifting over a VCF to another genome build and producing a properly header'd,
 * sorted and indexed VCF in one go.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = LiftoverVcf.USAGE_SUMMARY + LiftoverVcf.USAGE_DETAILS,
        oneLineSummary = LiftoverVcf.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class)
@DocumentedFeature
public class LiftoverVcf extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Lifts over a VCF file from one reference build to another.  ";
    static final String USAGE_DETAILS = "This tool adjusts the coordinates of variants within a VCF file to match a new reference. The " +
            "output file will be sorted and indexed using the target reference build. To be clear, REFERENCE_SEQUENCE should be the " +
            "<em>target</em> reference build. The tool is based on the UCSC liftOver tool (see: http://genome.ucsc.edu/cgi-bin/hgLiftOver) " +
            "and uses a UCSC chain file to guide its operation. <br /><br />" +
            "Note that records may be rejected because they cannot be lifted over or because of sequence incompatibilities between the " +
            "source and target reference genomes.  Rejected records will be emitted with filters to the REJECT file, using the source " +
            "genome coordinates.<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar LiftoverVcf \\<br />" +
            "     I=input.vcf \\<br />" +
            "     O=lifted_over.vcf \\<br />" +
            "     CHAIN=b37tohg19.chain \\<br />" +
            "     REJECT=rejected_variants.vcf \\<br />" +
            "     R=reference_sequence.fasta" +
            "</pre>" +
            "For additional information, please see: http://genome.ucsc.edu/cgi-bin/hgLiftOver" +
            "<hr />";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input VCF/BCF file to be lifted over.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output location to write the lifted over VCF/BCF to.")
    public File OUTPUT;

    @Argument(shortName = "C", doc = "The liftover chain file. See https://genome.ucsc.edu/goldenPath/help/chain.html for a description" +
            " of chain files.  See http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files.")
    public File CHAIN;

    @Argument(doc = "File to which to write rejected records.")
    public File REJECT;

    // Option on whether or not to provide a warning, or error message and exit if a missing contig is encountered
    @Argument(shortName = "WMC", doc = "Warn on missing contig.", optional = true)
    public boolean WARN_ON_MISSING_CONTIG = false;

    // Option on whether or not to write the original contig/position of the variant to the INFO field
    @Argument(doc = "Write the original contig/position for lifted variants to the INFO field.", optional = true)
    public boolean WRITE_ORIGINAL_POSITION = false;

    @Argument(doc = "The minimum percent match required for a variant to be lifted.", optional = true)
    public double LIFTOVER_MIN_MATCH = 1.0;

    @Argument(doc = "Allow INFO and FORMAT in the records that are not found in the header", optional = true)
    public boolean ALLOW_MISSING_FIELDS_IN_HEADER = false;

    // When a contig used in the chain is not in the reference, exit with this value instead of 0.
    protected static int EXIT_CODE_WHEN_CONTIG_NOT_IN_REFERENCE = 1;

    /**
     * Filter name to use when a target cannot be lifted over.
     */
    public static final String FILTER_CANNOT_LIFTOVER_INDEL = "ReverseComplementedIndel";

    /**
     * Filter name to use when a target cannot be lifted over.
     */
    public static final String FILTER_NO_TARGET = "NoTarget";

    /**
     * Filter name to use when a target is lifted over, but the reference allele doesn't match the new reference.
     */
    public static final String FILTER_MISMATCHING_REF_ALLELE = "MismatchedRefAllele";

    /**
     * Filter name to use when an indel cannot be lifted over since it straddles two intervals in a chain which means
     * that it is unclear what are the right alleles to be used.
     */
    public static final String FILTER_INDEL_STRADDLES_TWO_INTERVALS = "IndelStraddlesMultipleIntevals";

    /**
     * Filters to be added to the REJECT file.
     */
    private static final List<VCFFilterHeaderLine> FILTERS = CollectionUtil.makeList(
            new VCFFilterHeaderLine(FILTER_CANNOT_LIFTOVER_INDEL, "Indel falls into a reverse complemented region in the target genome."),
            new VCFFilterHeaderLine(FILTER_NO_TARGET, "Variant could not be lifted between genome builds."),
            new VCFFilterHeaderLine(FILTER_MISMATCHING_REF_ALLELE, "Reference allele does not match reference genome sequence after liftover."),
            new VCFFilterHeaderLine(FILTER_INDEL_STRADDLES_TWO_INTERVALS, "Indel is straddling multiple intervalss in the chain, and so the results are not well defined.")
    );

    /**
     * Attribute used to store the name of the source contig/chromosome prior to liftover.
     */
    public static final String ORIGINAL_CONTIG = "OriginalContig";

    /**
     * Attribute used to store the position of the variant on the source contig prior to liftover.
     */
    public static final String ORIGINAL_START = "OriginalStart";

    /**
     * Attribute used to store the position of the failed variant on the target contig prior to finding out that alleles do not match.
     */
    public static final String ATTEMPTED_LOCUS = "AttemptedLocus";

    /**
     * Metadata to be added to the Passing file.
     */
    private static final List<VCFInfoHeaderLine> ATTRS = CollectionUtil.makeList(
            new VCFInfoHeaderLine(ORIGINAL_CONTIG, 1, VCFHeaderLineType.String, "The name of the source contig/chromosome prior to liftover."),
            new VCFInfoHeaderLine(ORIGINAL_START, 1, VCFHeaderLineType.String, "The position of the variant on the source contig prior to liftover.")
    );

    private VariantContextWriter rejects;
    private final Log log = Log.getInstance(LiftoverVcf.class);
    private SortingCollection<VariantContext> sorter;

    private long failedLiftover = 0, failedAlleleCheck = 0;

    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new ReferenceArgumentCollection() {
            @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, common=false,
                    doc = "The reference sequence (fasta) for the TARGET genome build.  The fasta file must have an " +
                            "accompanying sequence dictionary (.dict file).")
            public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

            @Override
            public File getReferenceFile() {
                return REFERENCE_SEQUENCE;
            }
        };
    }

    // Stock main method
    public static void main(final String[] args) {
        new LiftoverVcf().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsReadable(CHAIN);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(REJECT);

        ////////////////////////////////////////////////////////////////////////
        // Setup the inputs
        ////////////////////////////////////////////////////////////////////////
        final LiftOver liftOver = new LiftOver(CHAIN);
        final VCFFileReader in = new VCFFileReader(INPUT, false);

        log.info("Loading up the target reference genome.");
        final ReferenceSequenceFileWalker walker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final Map<String, ReferenceSequence> refSeqs = new HashMap<>();
        for (final SAMSequenceRecord rec : walker.getSequenceDictionary().getSequences()) {
            refSeqs.put(rec.getSequenceName(), walker.get(rec.getSequenceIndex()));
        }
        CloserUtil.close(walker);

        ////////////////////////////////////////////////////////////////////////
        // Setup the outputs
        ////////////////////////////////////////////////////////////////////////
        final VCFHeader inHeader = in.getFileHeader();
        final VCFHeader outHeader = new VCFHeader(inHeader);
        outHeader.setSequenceDictionary(walker.getSequenceDictionary());
        if (WRITE_ORIGINAL_POSITION) {
            for (final VCFInfoHeaderLine line : ATTRS) outHeader.addMetaDataLine(line);
        }
        final VariantContextWriter out = new VariantContextWriterBuilder()
                .setOption(Options.INDEX_ON_THE_FLY)
                .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, ALLOW_MISSING_FIELDS_IN_HEADER)
                .setOutputFile(OUTPUT).setReferenceDictionary(walker.getSequenceDictionary()).build();
        out.writeHeader(outHeader);

        rejects = new VariantContextWriterBuilder().setOutputFile(REJECT)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .modifyOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER, ALLOW_MISSING_FIELDS_IN_HEADER)
                .build();
        final VCFHeader rejectHeader = new VCFHeader(in.getFileHeader());
        for (final VCFFilterHeaderLine line : FILTERS) rejectHeader.addMetaDataLine(line);

        rejectHeader.addMetaDataLine(new VCFInfoHeaderLine(ATTEMPTED_LOCUS,1, VCFHeaderLineType.String, "The locus of the variant in the TARGET prior to failing due to mismatching alleles."));

        rejects.writeHeader(rejectHeader);

        ////////////////////////////////////////////////////////////////////////
        // Read the input VCF, lift the records over and write to the sorting
        // collection.
        ////////////////////////////////////////////////////////////////////////
        long total = 0;
        log.info("Lifting variants over and sorting.");

        sorter = SortingCollection.newInstance(VariantContext.class,
                new VCFRecordCodec(outHeader, ALLOW_MISSING_FIELDS_IN_HEADER || VALIDATION_STRINGENCY != ValidationStringency.STRICT),
                outHeader.getVCFRecordComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);

        ProgressLogger progress = new ProgressLogger(log, 1000000, "read");
        // a mapping from original allele to reverse complemented allele
        final Map<Allele, Allele> reverseComplementAlleleMap = new HashMap<>(10);

        for (final VariantContext ctx : in) {
            ++total;
            final Interval source = new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd(), false, ctx.getContig() + ":" + ctx.getStart() + "-" + ctx.getEnd());
            final Interval target = liftOver.liftOver(source, LIFTOVER_MIN_MATCH);

            // target is null when there is no good liftover for the context. This happens either when it fall in a gap
            // where there isn't a chain, or if a large enough proportion of it is diminished by the "deletion" at the
            // end of each interval in a chain.
            if (target == null) {
                rejectVariant(ctx, FILTER_NO_TARGET);
                continue;
            }

            // the target is the lifted-over interval comprised of the start/stop of the variant context,
            // if the sizes of target and ctx do not match, it means that the interval grew or shrank during
            // liftover which must be due to straddling multiple intervals in the liftover chain.
            // This would invalidate the indel as it isn't clear what the resulting alleles should be.
            if (ctx.getReference().length() != target.length()){
                rejectVariant(ctx, FILTER_INDEL_STRADDLES_TWO_INTERVALS);
                continue;
            }

            final ReferenceSequence refSeq;

            // if the target is null OR (the target is reverse complemented AND the variant is a non-biallelic indel or mixed), then we cannot lift it over
            if (target.isNegativeStrand() && (ctx.isMixed() || ctx.isIndel() && !ctx.isBiallelic())) {
                rejectVariant(ctx, FILTER_CANNOT_LIFTOVER_INDEL);

            } else if (!refSeqs.containsKey(target.getContig())) {
                rejectVariant(ctx, FILTER_NO_TARGET);

                final String missingContigMessage = "Encountered a contig, " + target.getContig() + " that is not part of the target reference.";
                if (WARN_ON_MISSING_CONTIG) {
                    log.warn(missingContigMessage);
                } else {
                    log.error(missingContigMessage);
                    return EXIT_CODE_WHEN_CONTIG_NOT_IN_REFERENCE;
                }
            } else if (target.isNegativeStrand() && ctx.isIndel() && ctx.isBiallelic()) {
                refSeq = refSeqs.get(target.getContig());
                //flipping indels:

                final VariantContext flippedIndel = flipIndel(ctx, liftOver, refSeq);
                if (flippedIndel == null) {
                    throw new IllegalArgumentException("Unexpectedly found null VC. This should have not happened.");
                } else {
                    tryToAddVariant(flippedIndel, refSeq, reverseComplementAlleleMap, ctx);
                }
            } else {
                refSeq = refSeqs.get(target.getContig());
                final VariantContext liftedVariant = liftSimpleVariant(ctx, target);

                tryToAddVariant(liftedVariant, refSeq, reverseComplementAlleleMap, ctx);
            }
            progress.record(ctx.getContig(), ctx.getStart());
        }

        final NumberFormat pfmt = new DecimalFormat("0.0000%");
        final String pct = pfmt.format((failedLiftover + failedAlleleCheck) / (double) total);
        log.info("Processed ", total, " variants.");
        log.info(failedLiftover, " variants failed to liftover.");
        log.info(failedAlleleCheck, " variants lifted over but had mismatching reference alleles after lift over.");
        log.info(pct, " of variants were not successfully lifted over and written to the output.");

        rejects.close();
        in.close();

        ////////////////////////////////////////////////////////////////////////
        // Write the sorted outputs to the final output file
        ////////////////////////////////////////////////////////////////////////
        sorter.doneAdding();
        progress = new ProgressLogger(log, 1000000, "written");
        log.info("Writing out sorted records to final VCF.");

        for (final VariantContext ctx : sorter) {
            out.add(ctx);
            progress.record(ctx.getContig(), ctx.getStart());
        }
        out.close();
        sorter.cleanup();

        return 0;
    }

    private void rejectVariant(final VariantContext ctx, final String reason) {
        rejects.add(new VariantContextBuilder(ctx).filter(reason).make());
        failedLiftover++;
    }

    /**
     *  utility function to attempt to add a variant. Checks that the reference allele still matches the reference (which may have changed)
     *
     * @param vc new {@link VariantContext}
     * @param refSeq {@link ReferenceSequence} of new reference
     * @param alleleMap a {@link Map} mapping the old alleles to the new alleles (for fixing the genotypes)
     * @param source the original {@link VariantContext} to use for putting the original location information into vc
     * @return true if successful, false if failed due to mismatching reference allele.
     */
    private void tryToAddVariant(final VariantContext vc, final ReferenceSequence refSeq, final Map<Allele, Allele> alleleMap, final VariantContext source) {

        final VariantContextBuilder builder = new VariantContextBuilder(vc);

        builder.genotypes(fixGenotypes(source.getGenotypes(), alleleMap));
        builder.filters(source.getFilters());
        builder.log10PError(source.getLog10PError());
        builder.attributes(source.getAttributes());
        builder.id(source.getID());

        if (WRITE_ORIGINAL_POSITION) {
            builder.attribute(ORIGINAL_CONTIG, source.getContig());
            builder.attribute(ORIGINAL_START, source.getStart());
        }

        // Check that the reference allele still agrees with the reference sequence
        boolean mismatchesReference = false;
        for (final Allele allele : builder.getAlleles()) {
            if (allele.isReference()) {
                final byte[] ref = refSeq.getBases();
                final String refString = StringUtil.bytesToString(ref, vc.getStart() - 1, vc.getEnd() - vc.getStart() + 1);

                if (!refString.equalsIgnoreCase(allele.getBaseString())) {
                    mismatchesReference = true;
                }
                break;
            }
        }

        if (mismatchesReference) {
            rejects.add(new VariantContextBuilder(source)
                    .filter(FILTER_MISMATCHING_REF_ALLELE)
                    .attribute(ATTEMPTED_LOCUS, String.format("%s:%d-%d",vc.getContig(),vc.getStart(),vc.getEnd()))
                    .make());
            failedAlleleCheck++;
        } else {
            sorter.add(builder.make());
        }
    }

    /** liftsOver snps on either positive or negative strand and biallelic indels on positive strand only
     *
     * @param source
     * @param target
     * @return
     */
    protected static VariantContext liftSimpleVariant(final VariantContext source, final Interval target) {
        // Fix the alleles if we went from positive to negative strand

        if (target == null) {
            return null;
        }

        if (source.getReference().length() != target.length()) {
            return null;
        }
        final List<Allele> alleles = new ArrayList<>();

        for (final Allele oldAllele : source.getAlleles()) {
            if (target.isPositiveStrand() || oldAllele.isSymbolic()) {
                alleles.add(oldAllele);
            } else {
                final Allele fixedAllele = Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference());
                alleles.add(fixedAllele);
            }
        }

        // Build the new variant context
        final VariantContextBuilder builder = new VariantContextBuilder(
                source.getSource(),
                target.getContig(),
                target.getStart(),
                target.getEnd(),
                alleles);

        builder.id(source.getID());

        return builder.make();
    }

    /**
     * @param source            original variant context
     * @param liftOver          the LiftOver object to use for flipping
     * @param referenceSequence the reference sequence of the target
     * @return a flipped variant-context.
     */
    protected static VariantContext flipIndel(final VariantContext source, final LiftOver liftOver, final ReferenceSequence referenceSequence) {
        if (!source.isBiallelic()) return null;  //only supporting biallelic indels, for now.

        final Interval originalLocus = new Interval(source.getContig(), source.getStart(), source.getEnd());
        final Interval target = liftOver.liftOver(originalLocus);

        if (target == null) return null;
        if (!target.isNegativeStrand()) {
            throw new IllegalArgumentException("Expecting a variant the is lifted over with an inversion. Got " +
                    source + " maps to " + target.toString());
        }

        // a boolean to protect against trying to access the -1 position in the reference array
        final boolean addToStart = target.getStart() > 1;

        final Map<Allele, Allele> reverseComplementAlleleMap = new HashMap<>(2);

        reverseComplementAlleleMap.clear();
        final List<Allele> alleles = new ArrayList<>();

        for (final Allele oldAllele : source.getAlleles()) {
            // target.getStart is 1-based, reference bases are 0-based
            final StringBuilder alleleBuilder = new StringBuilder(target.getEnd() - target.getStart() + 1);

            if (addToStart) alleleBuilder.append((char) referenceSequence.getBases()[target.getStart() - 2]);
            alleleBuilder.append(SequenceUtil.reverseComplement(oldAllele.getBaseString().substring(1, oldAllele.length())));
            if (!addToStart) alleleBuilder.append((char) referenceSequence.getBases()[target.getEnd() - 1]);

            final Allele fixedAllele = Allele.create(alleleBuilder.toString(), oldAllele.isReference());
            alleles.add(fixedAllele);
            reverseComplementAlleleMap.put(oldAllele, fixedAllele);
        }

        final VariantContextBuilder builder = new VariantContextBuilder(source.getSource(),
                target.getContig(),
                target.getStart() - (addToStart ? 1 : 0),
                target.getEnd() - (addToStart ? 1 : 0),
                alleles);

        builder.id(source.getID());
        builder.attributes(source.getAttributes());

        builder.genotypes(fixGenotypes(source.getGenotypes(), reverseComplementAlleleMap));
        builder.filters(source.getFilters());
        builder.log10PError(source.getLog10PError());

        return leftAlignVariant(builder.make(), referenceSequence);
    }

    protected static GenotypesContext fixGenotypes(final GenotypesContext originals, final Map<Allele, Allele> alleleMap) {
        // optimization: if nothing needs to be fixed then don't bother
        if (alleleMap.isEmpty()) {
            return originals;
        }

        final GenotypesContext fixedGenotypes = GenotypesContext.create(originals.size());
        for (final Genotype genotype : originals) {
            final List<Allele> fixedAlleles = new ArrayList<>();
            for (final Allele allele : genotype.getAlleles()) {
                final Allele fixedAllele = alleleMap.getOrDefault(allele, allele);
                fixedAlleles.add(fixedAllele);
            }
            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
        }
        return fixedGenotypes;
    }

    /**
     *    Normalizes and left aligns a {@link VariantContext}.
     *
     *    Based on Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. (2015)
     *    Unified Representation of Genetic Variants. Bioinformatics.
     *
     * @param vc the {@link VariantContext} to be normalized
     * @param referenceSequence the {@link ReferenceSequence} of the same contig as vc
     * @return a new {@link VariantContext} which represents the same variation as vc but has
     * been normalized and left-aligned
     */
    protected static VariantContext leftAlignVariant(final VariantContext vc, final ReferenceSequence referenceSequence) {

        boolean changesInAlleles = true;

        int start = vc.getStart();
        int end = vc.getEnd();

        if (!vc.getContig().equals(referenceSequence.getName())) {
            throw new IllegalArgumentException("vc contig doesn't match that of supplied reference: " + vc.getContig() + " != " + referenceSequence.getName());
        }

        final Map<Allele, byte[]> alleleBasesMap = new HashMap<>();
        vc.getAlleles().forEach(a -> alleleBasesMap.put(a, a.getBases()));

        // 1. while changes in alleles do
        while (changesInAlleles) {

            changesInAlleles = false;
            // 2. if alleles end with the same nucleotide then
            if (alleleBasesMap.values().stream()
                    .collect(Collectors.groupingBy(a -> a[a.length - 1], Collectors.toSet()))
                    .size() == 1 && end > 1) {
                // 3. truncate rightmost nucleotide of each allele
                for (final Allele allele : alleleBasesMap.keySet()) {
                    alleleBasesMap.put(allele, truncateBase(alleleBasesMap.get(allele), true));
                }
                changesInAlleles = true;
                end--;
                // 4. end if
            }

            // 5. if there exists an empty allele then
            if (alleleBasesMap.values().stream()
                    .map(a -> a.length)
                    .anyMatch(l -> l == 0)) {
                // 6. extend alleles 1 nucleotide to the left
                for (final Allele allele : alleleBasesMap.keySet()) {
                    // the first -1 for zero-base (getBases) versus 1-based (variant position)
                    // another   -1 to get the base prior to the location of the start of the allele
                    final byte extraBase =  (start > 1) ?
                            referenceSequence.getBases()[start - 2] :
                            referenceSequence.getBases()[end];

                    alleleBasesMap.put(allele, extendOneBase(alleleBasesMap.get(allele), extraBase));
                }
                changesInAlleles = true;
                start--;

                // 7. end if
            }
        }

        // 8. while leftmost nucleotide of each allele are the same and all alleles have length 2 or more do
        while (alleleBasesMap.values().stream()
                .allMatch(a -> a.length >= 2) &&

                alleleBasesMap.values().stream()
                        .collect(Collectors.groupingBy(a -> a[0], Collectors.toSet()))
                        .size() == 1
                ) {

            //9. truncate the leftmost base of the alleles
            for (final Allele allele : alleleBasesMap.keySet()) {
                alleleBasesMap.put(allele, truncateBase(alleleBasesMap.get(allele), false));
            }
            start++;
        }

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.start(start);
        builder.stop(end);

        final Map<Allele, Allele> fixedAlleleMap = alleleBasesMap.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, me -> Allele.create(me.getValue(), me.getKey().isReference())));
        builder.alleles(fixedAlleleMap.values());
        builder.genotypes(fixGenotypes(vc.getGenotypes(), fixedAlleleMap));

        return builder.make();
    }

    private static byte[] truncateBase(final byte[] allele, final boolean truncateRightmost) {
        return Arrays.copyOfRange(allele, truncateRightmost ? 0 : 1, truncateRightmost ?
                allele.length - 1 :
                allele.length);
    }

    //creates a new byte array with the base added at the begining
    private static byte[] extendOneBase(final byte[] bases, final byte base) {

        final byte[] newBases = new byte[bases.length + 1];

        System.arraycopy(bases, 0, newBases, 1, bases.length);
        newBases[0] = base;

        return newBases;
    }
}