package net.sf.picard.vcf;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.CloserUtil;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
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

    @Option(shortName=StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc="Sequence dictionary to use when indexing the VCF.", optional = true)
    public File SEQUENCE_DICTIONARY;

    private static final Set<String> NO_SAMPLES = Collections.emptySet();
    
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
        if (SEQUENCE_DICTIONARY != null) IoUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
        IoUtil.assertFileIsWritable(OUTPUT);

	    final VCFFileReader reader = new VCFFileReader(INPUT);
	    final VCFHeader header = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder());
	    final SAMSequenceDictionary sequenceDictionary =
			    SEQUENCE_DICTIONARY != null
			            ? SAMFileReader.getSequenceDictionary(SEQUENCE_DICTIONARY)
					    : header.getSequenceDictionary();
	    if (CREATE_INDEX && sequenceDictionary == null) {
		    throw new PicardException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
	    }
	    final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
	    final VariantContextWriter writer = VariantContextWriterFactory.create(OUTPUT, sequenceDictionary, options);

	    writer.writeHeader(header);

        final ProgressLogger progress = new ProgressLogger(Log.getInstance(MakeSitesOnlyVcf.class), 10000);

	    final CloseableIterator<VariantContext> iterator = reader.iterator();
	    while (iterator.hasNext()) {
		    final VariantContext context = iterator.next();
		    writer.add(context.subContextFromSamples(
                    NO_SAMPLES, 
                    false // Do not re-derive the alleles from the new, subsetted genotypes: our site-only VCF should retain these values.
            ));
            progress.record(context.getChr(), context.getStart());
        }

	    CloserUtil.close(iterator);
	    CloserUtil.close(reader);
	    writer.close();

        return 0;
    }
}
