package net.sf.picard.vcf;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineParser;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
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
			"Splits an input VCF or BCF file into two VCF files, one for indel records and one for SNPs. The" +
			"headers of the two output files will be identical. An index file is created and a" +
			"sequence dictionary is required by default.";

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The VCF or BCF input file")
	public File INPUT;

	@Option(doc="The VCF or BCF file to which SNP records should be written. The file format is determined by file extension.")
	public File SNP_OUTPUT;

	@Option(doc="The VCF or BCF file to which indel records should be written. The file format is determined by file extension.")
	public File INDEL_OUTPUT;

	@Option(shortName="D", doc="The index sequence dictionary to use instead of the sequence dictionaries in the input files", optional = true)
	public File SEQUENCE_DICTIONARY;

    @Option(doc="If true an exception will be thrown if an event type other than SNP or indel is encountered")
    public Boolean STRICT = true;

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

		final VCFFileReader fileReader = new VCFFileReader(INPUT);
		final VCFHeader fileHeader = fileReader.getFileHeader();

		final SAMSequenceDictionary sequenceDictionary =
				SEQUENCE_DICTIONARY != null
						? SAMFileReader.getSequenceDictionary(SEQUENCE_DICTIONARY)
						: fileHeader.getSequenceDictionary();
		if (CREATE_INDEX && sequenceDictionary == null) {
			throw new PicardException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
		}

		final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);

		final VariantContextWriter snpWriter = VariantContextWriterFactory.create(SNP_OUTPUT, sequenceDictionary, options);
		final VariantContextWriter indelWriter = VariantContextWriterFactory.create(INDEL_OUTPUT, sequenceDictionary, options);
		snpWriter.writeHeader(fileHeader);
		indelWriter.writeHeader(fileHeader);

        int incorrectVariantCount = 0;

		final CloseableIterator<VariantContext> iterator = fileReader.iterator();
		while (iterator.hasNext()) {
			final VariantContext context = iterator.next();
			if (context.isIndel()) indelWriter.add(context);
			else if (context.isSNP()) snpWriter.add(context);
			else {
                if (STRICT) throw new IllegalStateException("Found a record with type " + context.getType().name());
                else incorrectVariantCount++;
            }

            progress.record(context.getChr(), context.getStart());
		}

        if (incorrectVariantCount > 0) {
            log.debug("Found " + incorrectVariantCount + " records that didn't match SNP or INDEL");
        }

		CloserUtil.close(iterator);
		CloserUtil.close(fileReader);
		snpWriter.close();
		indelWriter.close();

		return 0;
	}
}
