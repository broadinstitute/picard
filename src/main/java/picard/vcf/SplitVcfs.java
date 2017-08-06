package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;

/**
 * Splits the input VCF file into two, one for indels and one for SNPs. The headers of the two output
 * files will be identical.
 * <p/>
 * An index file is created for the output file by default. Using an output file name with a ".gz"
 * extension will create gzip-compressed output.
 */
@CommandLineProgramProperties(
        summary = SplitVcfs.USAGE_SUMMARY + SplitVcfs.USAGE_DETAILS,
        oneLineSummary = SplitVcfs.USAGE_SUMMARY,
        programGroup = VcfOrBcf.class
)
public class SplitVcfs extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Splits SNPs and INDELs into separate files.  ";
    static final String USAGE_DETAILS = "This tool reads in a VCF or BCF file and writes out the SNPs and INDELs it contains to separate " +
            "files. The headers of the two output files will be identical and index files will be created for both outputs. If records " +
            "other than SNPs or INDELs are present, set the STRICT option to \"false\", otherwise the tool will raise an exception and " +
            "quit. <br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar SplitVcfs \\<br />" +
            "      I=input.vcf \\<br />" +
            "      SNP_OUTPUT=snp.vcf \\<br />" +
            "      INDEL_OUTPUT=indel.vcf \\<br />" +
            "      STRICT=false" +
            "</pre>" +
            "<hr />" ;
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The VCF or BCF input file")
    public File INPUT;

    @Argument(doc = "The VCF or BCF file to which SNP records should be written. The file format is determined by file extension.")
    public File SNP_OUTPUT;

    @Argument(doc = "The VCF or BCF file to which indel records should be written. The file format is determined by file extension.")
    public File INDEL_OUTPUT;

    @Argument(shortName = "D", doc = "The index sequence dictionary to use instead of the sequence dictionaries in the input files", optional = true)
    public File SEQUENCE_DICTIONARY;

    @Argument(doc = "If true an exception will be thrown if an event type other than SNP or indel is encountered")
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
        IOUtil.assertFileIsReadable(INPUT);
        final ProgressLogger progress = new ProgressLogger(log, 10000);

        final VCFFileReader fileReader = new VCFFileReader(INPUT);
        final VCFHeader fileHeader = fileReader.getFileHeader();

        final SAMSequenceDictionary sequenceDictionary =
                SEQUENCE_DICTIONARY != null
                        ? SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(SEQUENCE_DICTIONARY).getSequenceDictionary()
                        : fileHeader.getSequenceDictionary();
        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
        }

        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setReferenceDictionary(sequenceDictionary)
                .clearOptions();
        if (CREATE_INDEX)
            builder.setOption(Options.INDEX_ON_THE_FLY);

        final VariantContextWriter snpWriter = builder.setOutputFile(SNP_OUTPUT).build();
        final VariantContextWriter indelWriter = builder.setOutputFile(INDEL_OUTPUT).build();
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

            progress.record(context.getContig(), context.getStart());
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
