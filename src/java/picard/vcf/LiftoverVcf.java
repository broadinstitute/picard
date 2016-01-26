package picard.vcf;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFRecordCodec;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tool for lifting over a VCF to another genome build and producing a properly header'd,
 * sorted and indexed VCF in one go.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = LiftoverVcf.USAGE_SUMMARY + LiftoverVcf.USAGE_DETAILS,
        usageShort = LiftoverVcf.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class
)
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
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input VCF/BCF file to be lifted over.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output location to write the lifted over VCF/BCF to.")
    public File OUTPUT;

    @Option(shortName="C", doc="The liftover chain file. See https://genome.ucsc.edu/goldenPath/help/chain.html for a description" +
            " of chain files.  See http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files.")
    public File CHAIN;

    @Option(doc="File to which to write rejected records.")
    public File REJECT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, common=false,
            doc = "The reference sequence (fasta) for the TARGET genome build.  The fasta file must have an " +
                    "accompanying sequence dictionary (.dict file).")
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    // Option on whether or not to provide a warning, or error message and exit if a missing contig is encountered
    @Option(shortName = "WMC", doc = "Warn on missing contig.", optional = true)
    public boolean WARN_ON_MISSING_CONTIG = false;

    // When a contig used in the chain is not in the reference, exit with this value instead of 0.
    protected static int EXIT_CODE_WHEN_CONTIG_NOT_IN_REFERENCE = 1;

    /** Filter name to use when a target cannot be lifted over. */
    public static final String FILTER_CANNOT_LIFTOVER_INDEL = "ReverseComplementedIndel";

    /** Filter name to use when a target cannot be lifted over. */
    public static final String FILTER_NO_TARGET = "NoTarget";

    /** Filter name to use when a target is lifted over, but the reference allele doens't match the new reference. */
    public static final String FILTER_MISMATCHING_REF_ALLELE = "MismatchedRefAllele";

    /** Filters to be added to the REJECT file. */
    private static final List<VCFFilterHeaderLine> FILTERS = CollectionUtil.makeList(
            new VCFFilterHeaderLine(FILTER_CANNOT_LIFTOVER_INDEL, "Indel falls into a reverse complemented region in the target genome."),
            new VCFFilterHeaderLine(FILTER_NO_TARGET, "Variant could not be lifted between genome builds."),
            new VCFFilterHeaderLine(FILTER_MISMATCHING_REF_ALLELE, "Reference allele does not match reference genome sequence after liftover.")
    );

    private final Log log = Log.getInstance(LiftoverVcf.class);

    // Stock main method
    public static void main(final String[] args) {
        new LiftoverVcf().instanceMainWithExit(args);
    }

    @Override protected int doWork() {
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
        final Map<String,byte[]> refSeqs = new HashMap<String,byte[]>();
        for (final SAMSequenceRecord rec: walker.getSequenceDictionary().getSequences()) {
            refSeqs.put(rec.getSequenceName(), walker.get(rec.getSequenceIndex()).getBases());
        }
        CloserUtil.close(walker);


        ////////////////////////////////////////////////////////////////////////
        // Setup the outputs
        ////////////////////////////////////////////////////////////////////////
        final VCFHeader inHeader = in.getFileHeader();
        final VCFHeader outHeader = new VCFHeader(inHeader);
        outHeader.setSequenceDictionary(walker.getSequenceDictionary());
        final VariantContextWriter out = new VariantContextWriterBuilder().setOption(Options.INDEX_ON_THE_FLY)
                .setOutputFile(OUTPUT).setReferenceDictionary(walker.getSequenceDictionary()).build();
        out.writeHeader(outHeader);

        final VariantContextWriter rejects = new VariantContextWriterBuilder().setOutputFile(REJECT).unsetOption(Options.INDEX_ON_THE_FLY).build();
        final VCFHeader rejectHeader = new VCFHeader(in.getFileHeader());
        for (final VCFFilterHeaderLine line : FILTERS) rejectHeader.addMetaDataLine(line);
        rejects.writeHeader(rejectHeader);


        ////////////////////////////////////////////////////////////////////////
        // Read the input VCF, lift the records over and write to the sorting
        // collection.
        ////////////////////////////////////////////////////////////////////////
        long failedLiftover = 0, failedAlleleCheck = 0, total = 0;
        log.info("Lifting variants over and sorting.");

        final SortingCollection<VariantContext> sorter = SortingCollection.newInstance(VariantContext.class,
                new VCFRecordCodec(outHeader, VALIDATION_STRINGENCY != ValidationStringency.STRICT),
                outHeader.getVCFRecordComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);

        ProgressLogger progress = new ProgressLogger(log, 1000000, "read");
        // a mapping from original allele to reverse complemented allele
        final Map<Allele, Allele> reverseComplementAlleleMap = new HashMap<Allele, Allele>(10);

        for (final VariantContext ctx : in) {
            ++total;
            final Interval source = new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd(), false, ctx.getContig() + ":" + ctx.getStart() + "-" + ctx.getEnd());
            final Interval target = liftOver.liftOver(source, 1.0);

            // if the target is null OR (the target is reverse complemented AND the variant is an indel or mixed), then we cannot lift it over
            if (target == null || (target.isNegativeStrand() && (ctx.isMixed() || ctx.isIndel()))) {
                final String reason = (target == null) ? FILTER_NO_TARGET : FILTER_CANNOT_LIFTOVER_INDEL;
                rejects.add(new VariantContextBuilder(ctx).filter(reason).make());
                failedLiftover++;
            } else if (!refSeqs.containsValue(target.getContig())) {
                rejects.add(new VariantContextBuilder(ctx).filter(FILTER_NO_TARGET).make());
                failedLiftover++;

                String missingContigMessage = "Encountered a contig, " + target.getContig() + " that is not part of the target reference.";
                if(WARN_ON_MISSING_CONTIG) {
                    log.warn(missingContigMessage);
                } else {
                    log.error(missingContigMessage);
                    return EXIT_CODE_WHEN_CONTIG_NOT_IN_REFERENCE;
                }
            } else {
                // Fix the alleles if we went from positive to negative strand
                reverseComplementAlleleMap.clear();
                final List<Allele> alleles = new ArrayList<Allele>();

                for (final Allele oldAllele : ctx.getAlleles()) {
                    if (target.isPositiveStrand() || oldAllele.isSymbolic()) {
                        alleles.add(oldAllele);
                    }
                    else {
                        final Allele fixedAllele = Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference());
                        alleles.add(fixedAllele);
                        reverseComplementAlleleMap.put(oldAllele, fixedAllele);
                    }
                }

                // Build the new variant context
                final VariantContextBuilder builder = new VariantContextBuilder(
                        ctx.getSource(),
                        target.getContig(),
                        target.getStart(),
                        target.getEnd(),
                        alleles);

                builder.id(ctx.getID());
                builder.attributes(ctx.getAttributes());
                builder.genotypes(fixGenotypes(ctx.getGenotypes(), reverseComplementAlleleMap));
                builder.filters(ctx.getFilters());
                builder.log10PError(ctx.getLog10PError());

                // Check that the reference allele still agrees with the reference sequence
                boolean mismatchesReference = false;
                for (final Allele allele : builder.getAlleles()) {
                    if (allele.isReference()) {
                        final byte[] ref = refSeqs.get(target.getContig());
                        final String refString = StringUtil.bytesToString(ref, target.getStart()-1, target.length());

                        if (!refString.equalsIgnoreCase(allele.getBaseString())) {
                            mismatchesReference = true;
                        }

                        break;
                    }
                }

                if (mismatchesReference) {
                    rejects.add(new VariantContextBuilder(ctx).filter(FILTER_MISMATCHING_REF_ALLELE).make());
                    failedAlleleCheck++;
                }
                else {
                    sorter.add(builder.make());
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

    protected static GenotypesContext fixGenotypes(final GenotypesContext originals, final Map<Allele, Allele> reverseComplementAlleleMap) {
        // optimization: if nothing needs to be fixed then don't bother
        if ( reverseComplementAlleleMap.isEmpty() ) {
            return originals;
        }

        final GenotypesContext fixedGenotypes = GenotypesContext.create(originals.size());
        for ( final Genotype genotype : originals ) {
            final List<Allele> fixedAlleles = new ArrayList<Allele>();
            for ( final Allele allele : genotype.getAlleles() ) {
                final Allele fixedAllele = reverseComplementAlleleMap.containsKey(allele) ? reverseComplementAlleleMap.get(allele) : allele;
                fixedAlleles.add(fixedAllele);
            }
            fixedGenotypes.add(new GenotypeBuilder(genotype).alleles(fixedAlleles).make());
        }
        return fixedGenotypes;
    }
}