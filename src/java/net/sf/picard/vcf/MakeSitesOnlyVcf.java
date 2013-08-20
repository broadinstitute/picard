package net.sf.picard.vcf;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.File;
import java.util.Collections;
import java.util.Set;

/**
 * Writes out a VCF that contains all the site-level information for all records in the input VCF and no per-sample information.
 *
 * @author Tim Fennell
 */
public class MakeSitesOnlyVcf extends CommandLineProgram {
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF or BCF")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF or BCF to emit without per-sample info.")
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc="Sequence dictionary to use when indexing the VCF.")
    public File SEQUENCE_DICTIONARY;

    private static final Set<String> NO_SAMPLES = Collections.emptySet();
    
    // Stock main method
    public static void main(final String[] args) {
        new MakeSitesOnlyVcf().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IoUtil.assertFileIsWritable(OUTPUT);

        final VariantContextIterator in = VariantContextIteratorFactory.create(INPUT);
        final VariantContextWriter out = VariantContextWriterFactory.create(OUTPUT, VariantContextUtils.getSequenceDictionary(SEQUENCE_DICTIONARY));

        final VCFHeader header = new VCFHeader(in.getHeader());
        out.writeHeader(header);

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(MakeSitesOnlyVcf.class), 10000);

        while (in.hasNext()) {
            final VariantContext ctx = in.next();
            out.add(ctx.subContextFromSamples(
                    NO_SAMPLES, 
                    false // Do not re-derive the alleles from the new, subsetted genotypes: our site-only VCF should retain these values.
            ));
            progress.record(ctx.getChr(), ctx.getStart());
        }

        out.close();
        in.close();

        return 0;
    }
}
