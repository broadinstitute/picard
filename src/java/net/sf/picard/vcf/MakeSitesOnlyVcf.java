package net.sf.picard.vcf;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterBuilder;
import org.broadinstitute.variant.vcf.*;

import java.io.File;
import java.util.*;

/**
 * Writes out a VCF that contains all the site-level information for all records in the input VCF and no per-sample information.
 *
 * @author Tim Fennell
 */
public class MakeSitesOnlyVcf extends CommandLineProgram {
    @Usage
    public final String usage = "Reads a VCF/VCF.gz/BCF and removes all genotype information from it while retaining " +
            "all site level information, including annotations based on genotypes (e.g. AN, AF). Output an be " +
            "any support variant format including .vcf, .vcf.gz or .bcf.";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF or BCF")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF or BCF to emit without per-sample info.")
    public File OUTPUT;

    @Option(shortName="S", doc="Optionally one or more samples to retain when building the 'sites-only' VCF.", optional=true)
    public Set<String> SAMPLE = new TreeSet<String>();

    // Stock main method
    public static void main(final String[] args) {
        new MakeSitesOnlyVcf().instanceMainWithExit(args);
    }

	public MakeSitesOnlyVcf() {
		CREATE_INDEX = true;
	}

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

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
            progress.record(site.getChr(), site.getStart());
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
