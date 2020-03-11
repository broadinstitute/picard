package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.util.Set;
import java.util.TreeSet;

/**
 * Creates a VCF that contains all the site-level information for all records in the input VCF but no genotype information.
 *
 * <h3> Summary </h3>
 * This tool reads a VCF/VCF.gz/BCF and removes all genotype information from it while retaining all site level information,
 * including annotations based on genotypes (e.g. AN, AF). Output can be any supported variant format including .vcf,
 * .vcf.gz or .bcf.
 *
 * <h3> Inputs</h3>
 * <ul>
 *     <li> Input VCF or BCF file containing genotype and site-level information. </li>
 *     <li> Output VCF or BCF file containing only site-level information. </li>
 *     <li> [Optional] Names of one or more samples to include in the output VCF. </li>
 * </ul>
 *
 * <h3>Usage example:</h3>
 * <pre>
 *     java -jar picard.jar MakeSitesOnlyVcf \
 *      INPUT=input_variants.vcf \
 *      OUTPUT=output_variants.vcf
 * </pre>
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = MakeSitesOnlyVcf.USAGE_DETAILS,
        oneLineSummary = MakeSitesOnlyVcf.USAGE_SUMMARY,
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
public class MakeSitesOnlyVcf extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Creates a VCF that contains all the site-level information for all records in the input VCF but no genotype information.";
    static final String USAGE_DETAILS = "This tool reads a VCF/VCF.gz/BCF and removes all genotype information from it while retaining" +
            "all site level information, including annotations based on genotypes (e.g. AN, AF). Output can be" +
            "any supported variant format including .vcf, .vcf.gz or .bcf. <br /><br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar MakeSitesOnlyVcf \\ <br />" +
            "      INPUT=input_variants.vcf \\ <br />" +
            "      OUTPUT=output_variants.vcf" +
            "</pre>";

    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF or BCF containing genotype and site-level information.")
    public File INPUT;

    @Argument(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF or BCF file containing only site-level information.")
    public File OUTPUT;

    @Argument(shortName="S", doc="Names of one or more samples to include in the output VCF.", optional=true)
    public Set<String> SAMPLE = new TreeSet<String>();

    public MakeSitesOnlyVcf() {
        CREATE_INDEX = true;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final VCFFileReader reader = new VCFFileReader(INPUT, false);
        final VCFHeader inputVcfHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder());
        final SAMSequenceDictionary sequenceDictionary = inputVcfHeader.getSequenceDictionary();

        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
        }

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(MakeSitesOnlyVcf.class), 10000);

        // Setup the site-only file writer
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);
        if (CREATE_INDEX)
            builder.setOption(Options.INDEX_ON_THE_FLY);
        else
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        final VariantContextWriter writer = builder.build();

        final VCFHeader header = new VCFHeader(inputVcfHeader.getMetaDataInInputOrder(), SAMPLE);
        writer.writeHeader(header);

        // Go through the input, strip the records and write them to the output
        final CloseableIterator<VariantContext> iterator = reader.iterator();
        while (iterator.hasNext()) {
            final VariantContext full = iterator.next();
            final VariantContext site = subsetToSamplesWithOriginalAnnotations(full, SAMPLE);
            writer.add(site);
            progress.record(site.getContig(), site.getStart());
        }

        CloserUtil.close(iterator);
        CloserUtil.close(reader);
        writer.close();

        return 0;
    }

    /** Makes a new VariantContext with only the desired samples. */
    private static VariantContext subsetToSamplesWithOriginalAnnotations(final VariantContext ctx, final Set<String> samples) {
        final VariantContextBuilder builder = new VariantContextBuilder(ctx);
        final GenotypesContext newGenotypes = ctx.getGenotypes().subsetToSamples(samples);
        builder.alleles(ctx.getAlleles());
        return builder.genotypes(newGenotypes).make();
    }
}
