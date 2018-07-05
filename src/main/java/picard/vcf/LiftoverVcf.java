package picard.vcf;

/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.argumentcollections.ReferenceArgumentCollection;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import picard.util.LiftoverUtils;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

/**
 * <h3>Summary</h3>
 * Tool for "lifting over" a VCF from one genome build to another, producing a properly headered,
 * sorted and indexed VCF in one go.
 *
 * <h3>Details</h3>
 * This tool adjusts the coordinates of variants within a VCF file to match a new reference. The
 * output file will be sorted and indexed using the target reference build. To be clear, {@link #REFERENCE_SEQUENCE} should be the
 * <em>target</em> reference build (that is, the "new" one). The tool is based on the <a href="http://genome.ucsc.edu/cgi-bin/hgLiftOver">UCSC LiftOver tool</a>
 * and uses a UCSC chain file to guide its operation. <br />
 *
 * For each variant, the tool will look for the target coordinate, reverse-complement and left-align the variant if needed,
 * and, in the case that the reference and alternate alleles of a SNP have been swapped in the new genome build, it will
 * adjust the SNP, and correct AF-like INFO fields and the relevant genotypes.
 * <br />
 *
 * <h3>Example</h3>
 * <pre>
 * java -jar picard.jar LiftoverVcf \\
 *     I=input.vcf \\
 *     O=lifted_over.vcf \\
 *     CHAIN=b37tohg38.chain \\
 *     REJECT=rejected_variants.vcf \\
 *     R=reference_sequence.fasta
 * </pre>
 * <h3>Caveats</h3>
 * <h4>Rejected Records</h4>
 * Records may be rejected because they cannot be lifted over or because of sequence incompatibilities between the
 * source and target reference genomes.  Rejected records will be emitted to the {@link #REJECT} file using the source
 * genome build coordinates. The reason for the rejection will be stated in the FILTER field, and more detail may be placed
 * in the INFO field.
 * <h4>Memory Use</h4>
 * LiftOverVcf sorts the output using a {@link htsjdk.samtools.util.SortingCollection} which relies on {@link #MAX_RECORDS_IN_RAM}
 * to specify how many (vcf) records to hold in memory before "spilling" to disk. The default value is reasonable when sorting SAM files,
 * but not for VCFs as there is no good default due to the dependence on the number of samples and amount of information in the INFO and FORMAT
 * fields. Consider lowering to 100,000 or even less if you have many genotypes.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = LiftoverVcf.USAGE_SUMMARY + LiftoverVcf.USAGE_DETAILS,
        oneLineSummary = LiftoverVcf.USAGE_SUMMARY,
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
public class LiftoverVcf extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Lifts over a VCF file from one reference build to another.  ";
    static final String USAGE_DETAILS = "<h3>Summary</h3>\n" +
            "Tool for \"lifting over\" a VCF from one genome build to another, producing a properly headered, " +
            "sorted and indexed VCF in one go.\n" +
            "\n" +
            "<h3>Details</h3>\n" +
            "This tool adjusts the coordinates of variants within a VCF file to match a new reference. The " +
            "output file will be sorted and indexed using the target reference build. To be clear, REFERENCE_SEQUENCE should be the " +
            "<em>target</em> reference build (that is, the \"new\" one). The tool is based on the UCSC LiftOver tool (see http://genome.ucsc.edu/cgi-bin/hgLiftOver) " +
            "and uses a UCSC chain file to guide its operation.\n" +
            "\n" +
            "For each variant, the tool will look for the target coordinate, reverse-complement and left-align the variant if needed, " +
            "and, in the case that the reference and alternate alleles of a SNP have been swapped in the new genome build, it will " +
            "adjust the SNP, and correct AF-like INFO fields and the relevant genotypes." +
            "\n" +
            "\n" +
            "<h3>Example</h3>\n" +
            "java -jar picard.jar LiftoverVcf \\\n" +
            "    I=input.vcf \\\n" +
            "    O=lifted_over.vcf \\\n" +
            "    CHAIN=b37tohg38.chain \\\n" +
            "    REJECT=rejected_variants.vcf \\\n" +
            "    R=reference_sequence.fasta\n" +
            "\n" +
            "<h3>Caveats</h3>\n" +
            "<h4>Rejected Records</h4>\n" +
            "Records may be rejected because they cannot be lifted over or because of sequence incompatibilities between the " +
            "source and target reference genomes.  Rejected records will be emitted to the REJECT file using the source " +
            "genome build coordinates. The reason for the rejection will be stated in the FILTER field, and more detail may be placed " +
            "in the INFO field.\n" +
            "<h4>Memory Use</h4>\n" +
            "LiftOverVcf sorts the output using a \"SortingCollection\" which relies on MAX_RECORDS_IN_RAM " +
            "to specify how many (vcf) records to hold in memory before \"spilling\" to disk. The default value is reasonable when sorting SAM files, " +
            "but not for VCFs as there is no good default due to the dependence on the number of samples and amount of information in the INFO and FORMAT " +
            "fields. Consider lowering to 100,000 or even less if you have many genotypes.\n";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input VCF/BCF file to be lifted over.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output location for the lifted over VCF/BCF.")
    public File OUTPUT;

    @Argument(shortName = "C", doc = "The liftover chain file. See https://genome.ucsc.edu/goldenPath/help/chain.html for a description" +
            " of chain files.  See http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files.")
    public File CHAIN;

    @Argument(doc = "File to which to write rejected records.")
    public File REJECT;

    // Option on whether or not to provide a warning, or error message and exit if a missing contig is encountered
    @Argument(shortName = "WMC", doc = "Warn on missing contig.", optional = true)
    public boolean WARN_ON_MISSING_CONTIG = false;

    @Argument(shortName = "LFI", doc = "If true, intervals failing due to match below LIFTOVER_MIN_MATCH will be logged as a warning to the console.", optional = true)
    public boolean LOG_FAILED_INTERVALS = true;

    // Option on whether or not to write the original contig/position of the variant to the INFO field
    @Argument(doc = "Write the original contig/position for lifted variants to the INFO field.", optional = true)
    public boolean WRITE_ORIGINAL_POSITION = false;

    // Option on whether or not to write the original alleles of the variant to the INFO field
    @Argument(doc = "Write the original alleles for lifted variants to the INFO field.  If the alleles are identical, this attribute will be omitted.", optional = true)
    public boolean WRITE_ORIGINAL_ALLELES = false;

    @Argument(doc = "The minimum percent match required for a variant to be lifted.", optional = true)
    public double LIFTOVER_MIN_MATCH = 1.0;

    @Argument(doc = "Allow INFO and FORMAT in the records that are not found in the header", optional = true)
    public boolean ALLOW_MISSING_FIELDS_IN_HEADER = false;


    @Argument(doc = "If the REF allele of the lifted site does not match the target genome, that variant is normally rejected. " +
            "For bi-allelic SNPs, if this is set to true and the ALT allele equals the new REF allele, the REF and ALT alleles will be swapped.  This can rescue " +
            "some variants; however, do this carefully as some annotations may become invalid, such as any that are alelle-specifc.  See also TAGS_TO_REVERSE and TAGS_TO_DROP.", optional = true)
    public boolean RECOVER_SWAPPED_REF_ALT = false;

    @Argument(doc = "INFO field annotations that behave like an Allele Frequency and should be transformed with x->1-x " +
            "when swapping reference with variant alleles.", optional = true)
    public Collection<String> TAGS_TO_REVERSE = new ArrayList<>(LiftoverUtils.DEFAULT_TAGS_TO_REVERSE);

    @Argument(doc = "INFO field annotations that should be deleted when swapping reference with variant alleles.", optional = true)
    public Collection<String> TAGS_TO_DROP = new ArrayList<>(LiftoverUtils.DEFAULT_TAGS_TO_DROP);

    // When a contig used in the chain is not in the reference, exit with this value instead of 0.
    public static int EXIT_CODE_WHEN_CONTIG_NOT_IN_REFERENCE = 1;

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
            new VCFFilterHeaderLine(FILTER_INDEL_STRADDLES_TWO_INTERVALS, "Reference allele in Indel is straddling multiple intervals in the chain, and so the results are not well defined.")
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
     * Attribute used to store the list of original alleles (including REF), in the original order prior to liftover.
     */
    public static final String ORIGINAL_ALLELES = "OriginalAlleles";

    /**
     * Attribute used to store the position of the failed variant on the target contig prior to finding out that alleles do not match.
     */
    public static final String ATTEMPTED_LOCUS = "AttemptedLocus";

    /**
     * Metadata to be added to the Passing file.
     */
    private static final List<VCFInfoHeaderLine> ATTRS = CollectionUtil.makeList(
            new VCFInfoHeaderLine(ORIGINAL_CONTIG, 1, VCFHeaderLineType.String, "The name of the source contig/chromosome prior to liftover."),
            new VCFInfoHeaderLine(ORIGINAL_START, 1, VCFHeaderLineType.String, "The position of the variant on the source contig prior to liftover."),
            new VCFInfoHeaderLine(ORIGINAL_ALLELES, VCFHeaderLineCount.R, VCFHeaderLineType.String, "A list of the original alleles (including REF) of the variant prior to liftover.  If the alleles were not changed during liftover, this attribute will be omitted.")
            );

    private VariantContextWriter rejects;
    private final Log log = Log.getInstance(LiftoverVcf.class);
    private SortingCollection<VariantContext> sorter;

    private long failedLiftover = 0, failedAlleleCheck = 0, totalTrackedAsSwapRefAlt = 0;
    private Map<String, Long> rejectsByContig = new TreeMap<>();
    private Map<String, Long> liftedByDestContig = new TreeMap<>();
    private Map<String, Long> liftedBySourceContig = new TreeMap<>();

    @Override
    protected ReferenceArgumentCollection makeReferenceArgumentCollection() {
        return new ReferenceArgumentCollection() {
            @Argument(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, common=false,
                    doc = "The reference sequence (fasta) for the TARGET genome build (i.e., the new one.  The fasta file must have an " +
                                "accompanying sequence dictionary (.dict file).")
            public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

            @Override
            public File getReferenceFile() {
                return REFERENCE_SEQUENCE;
            }
        };
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
        liftOver.setShouldLogFailedIntervalsBelowThreshold(LOG_FAILED_INTERVALS);

        final VCFFileReader in = new VCFFileReader(INPUT, false);

        log.info("Loading up the target reference genome.");
        final ReferenceSequenceFileWalker walker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final Map<String, ReferenceSequence> refSeqs = new HashMap<>();
        // check if sequence dictionary exists
        if (walker.getSequenceDictionary() == null) {
            log.error("Reference " + REFERENCE_SEQUENCE.getAbsolutePath() + " must have an associated Dictionary .dict file in the same directory.");
            return 1;
        }
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

        outHeader.addMetaDataLine(new VCFInfoHeaderLine(LiftoverUtils.SWAPPED_ALLELES, 0, VCFHeaderLineType.Flag,
                "The REF and the ALT alleles have been swapped in liftover due to changes in the reference. " +
                "It is possible that not all INFO annotations reflect this swap, and in the genotypes, " +
                "only the GT, PL, and AD fields have been modified. You should check the TAGS_TO_REVERSE parameter that was used " +
                        "during the LiftOver to be sure."));
        outHeader.addMetaDataLine(new VCFInfoHeaderLine(LiftoverUtils.REV_COMPED_ALLELES, 0, VCFHeaderLineType.Flag,
                "The REF and the ALT alleles have been reverse complemented in liftover since the mapping from the " +
                        "previous reference to the current one was on the negative strand."));

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
        log.info("Lifting variants over and sorting (not yet writing the output file.)");

        sorter = SortingCollection.newInstance(VariantContext.class,
                new VCFRecordCodec(outHeader, ALLOW_MISSING_FIELDS_IN_HEADER || VALIDATION_STRINGENCY != ValidationStringency.STRICT),
                outHeader.getVCFRecordComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);

        ProgressLogger progress = new ProgressLogger(log, 1000000, "read");

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
            if (ctx.getReference().length() != target.length()) {
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
            } else {
                refSeq = refSeqs.get(target.getContig());

                final VariantContext liftedVC = LiftoverUtils.liftVariant(ctx, target, refSeq, WRITE_ORIGINAL_POSITION, WRITE_ORIGINAL_ALLELES);
                // the liftedVC can be null if the liftover fails because of a problem with reverse complementing
                if (liftedVC == null) {
                    rejectVariant(ctx, FILTER_CANNOT_LIFTOVER_INDEL);
                } else {
                    tryToAddVariant(liftedVC, refSeq, ctx);
                }
            }
            progress.record(ctx.getContig(), ctx.getStart());
        }

        final NumberFormat pfmt = new DecimalFormat("0.0000%");
        final String pct = pfmt.format((failedLiftover + failedAlleleCheck) / (double) total);
        log.info("Processed ", total, " variants.");
        log.info(failedLiftover, " variants failed to liftover.");
        log.info(failedAlleleCheck, " variants lifted over but had mismatching reference alleles after lift over.");
        log.info(pct, " of variants were not successfully lifted over and written to the output.");

        final Set<String> contigUnion = new TreeSet<>();
        contigUnion.addAll(liftedBySourceContig.keySet());
        contigUnion.addAll(rejectsByContig.keySet());

        log.info("liftover success by source contig:");
        for (String contig : contigUnion) {
            final long success = liftedBySourceContig.getOrDefault(contig, 0L);
            final long fail = rejectsByContig.getOrDefault(contig, 0L);
            final String liftPct = pfmt.format((double)success / (double)(success + fail));

            log.info(contig, ": ", success, " / ", (success + fail), " (", liftPct, ")");
        }

        log.info("lifted variants by target contig:");
        for (String contig : liftedByDestContig.keySet()) {
            log.info(contig, ": ", liftedByDestContig.get(contig));
        }
        if (liftedByDestContig.isEmpty()) {
            log.info("no successfully lifted variants");
        }

        if (RECOVER_SWAPPED_REF_ALT) {
            log.info(totalTrackedAsSwapRefAlt, " variants were lifted by swapping REF/ALT alleles.");
        } else {
            log.warn(totalTrackedAsSwapRefAlt, " variants with a swapped REF/ALT were identified, but were not recovered.  See RECOVER_SWAPPED_REF_ALT and associated caveats.");
        }

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
        trackLiftedVariantContig(rejectsByContig, ctx.getContig());
    }

    private void trackLiftedVariantContig(Map<String, Long> map, String contig) {
        Long val = map.get(contig);
        if (val == null) {
            val = 0L;
        }

        map.put(contig, ++val);
    }

    private void addAndTrack(final VariantContext toAdd, final VariantContext source) {
        trackLiftedVariantContig(liftedBySourceContig, source.getContig());
        trackLiftedVariantContig(liftedByDestContig, toAdd.getContig());
        sorter.add(toAdd);
    }

    /**
     *  utility function to attempt to add a variant. Checks that the reference allele still matches the reference (which may have changed)
     *
     * @param vc new {@link VariantContext}
     * @param refSeq {@link ReferenceSequence} of new reference
     * @param source the original {@link VariantContext} to use for putting the original location information into vc
     * @return true if successful, false if failed due to mismatching reference allele.
     */
    private void tryToAddVariant(final VariantContext vc, final ReferenceSequence refSeq, final VariantContext source) {
        if (!refSeq.getName().equals(vc.getContig())) {
            throw new IllegalStateException("The contig of the VariantContext, " + vc.getContig() + ", doesnt match the ReferenceSequence: " + refSeq.getName());
        }

        // Check that the reference allele still agrees with the reference sequence
        boolean mismatchesReference = false;
        for (final Allele allele : vc.getAlleles()) {
            if (allele.isReference()) {
                final byte[] ref = refSeq.getBases();
                final String refString = StringUtil.bytesToString(ref, vc.getStart() - 1, vc.getEnd() - vc.getStart() + 1);

                if (!refString.equalsIgnoreCase(allele.getBaseString())) {
                    // consider that the ref and the alt may have been swapped in a simple biallelic SNP
                    if (vc.isBiallelic() && vc.isSNP() && refString.equalsIgnoreCase(vc.getAlternateAllele(0).getBaseString())) {
                        if (RECOVER_SWAPPED_REF_ALT) {
                            totalTrackedAsSwapRefAlt++;
                            addAndTrack(LiftoverUtils.swapRefAlt(vc, TAGS_TO_REVERSE, TAGS_TO_DROP), source);
                            return;
                        } else {
                            totalTrackedAsSwapRefAlt++;
                        }
                    }
                    mismatchesReference = true;
                }
                break;
            }
        }

        if (mismatchesReference) {
            rejects.add(new VariantContextBuilder(source)
                    .filter(FILTER_MISMATCHING_REF_ALLELE)
                    .attribute(ATTEMPTED_LOCUS, String.format("%s:%d-%d", vc.getContig(), vc.getStart(), vc.getEnd()))
                    .make());
            failedAlleleCheck++;
            trackLiftedVariantContig(rejectsByContig, source.getContig());
        } else {
            addAndTrack(vc, source);
        }
    }
}