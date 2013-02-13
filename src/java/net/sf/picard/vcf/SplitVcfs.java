package net.sf.picard.vcf;

import net.sf.picard.cmdline.CommandLineParser;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.io.File;
import java.util.EnumSet;

/**
 * Splits the input VCF file into two, one for indels and one for SNPs. The headers of the two output
 * files will be identical.
 *
 * An index file is created for the output file by default. Using an output file name with a ".gz"
 * extension will create gzip-compressed output.
 */
public class SplitVcfs extends CommandLineProgram {

	@Usage
	public final String USAGE =
			CommandLineParser.getStandardUsagePreamble(getClass()) +
			"Splits an input VCF file into two files, one for indel records and one for SNPs.\n";

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc="VCF input file")
	public File INPUT;

	@Option(doc="VCF file to which SNP records should be written")
	public File SNP_OUTPUT;

	@Option(doc="VCF file to which indel records should be written")
	public File INDEL_OUTPUT;

	@Option(shortName="D", doc="The index sequence dictionary (required if CREATE_INDEX=true)", optional = true)
	public File SEQUENCE_DICTIONARY;

	private final Log log = Log.getInstance(SplitVcfs.class);

	public static void main(final String[] argv) {
		new SplitVcfs().instanceMainWithExit(argv);
	}

	public SplitVcfs() {
		this.CREATE_INDEX = true;
	}

	@Override
	protected int doWork() {
		IoUtil.assertFileIsReadable(INPUT);
		final ProgressLogger progress = new ProgressLogger(log, 10000);
		final VariantContextIterator variantIterator = new VariantContextIterator(INPUT);
		final VCFHeader header = variantIterator.getHeader();

		final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
		final SAMSequenceDictionary sequenceDictionary =
				SEQUENCE_DICTIONARY != null ? VariantContextUtils.getSequenceDictionary(SEQUENCE_DICTIONARY) : null;

		final VariantContextWriter snpOutput = VariantContextUtils.getConditionallyCompressingWriter(SNP_OUTPUT, sequenceDictionary, options);
		final VariantContextWriter indelOutput = VariantContextUtils.getConditionallyCompressingWriter(INDEL_OUTPUT, sequenceDictionary, options);
		snpOutput.writeHeader(header);
		indelOutput.writeHeader(header);

		while (variantIterator.hasNext()) {
			final VariantContext context = variantIterator.next();
			if (context.isIndel()) indelOutput.add(context);
			else if (context.isSNP()) snpOutput.add(context);
			else throw new IllegalStateException("Found a record with type " + context.getType().name());
			progress.record(context.getChr(), context.getStart());
		}

		snpOutput.close();
		indelOutput.close();

		return 0;
	}

	protected String[] customCommandLineValidation() {
		if (this.CREATE_INDEX && (this.SEQUENCE_DICTIONARY == null)) {
			return new String[] { "If CREATE_INDEX is set a sequence dictionary must be specified." };
		}
		return null;
	}
}
