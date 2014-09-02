package picard.vcf;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFRecordCodec;
import htsjdk.variant.vcf.VCFUtils;
import picard.PicardException;
import picard.cmdline.CommandLineParser;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

/**
 * Sorts one or more VCF files according to the order of the contigs in the header/sequence dictionary and then
 * by coordinate.  Can accept an external dictionary. If no external dictionary is supplied, multiple inputs' headers must have
 * the same sequence dictionaries
 *
 */
public class SortVcf extends CommandLineProgram {
    @Usage
    public String USAGE =
            CommandLineParser.getStandardUsagePreamble(getClass()) +
            "Sorts one or more VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate. " +
            "Can accept an external sequence dictionary. If no external dictionary is supplied, multiple inputs' headers must have " +
            "the same sequence dictionaries\n";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF(s) to be sorted. Multiple inputs must have the same sample names (in order)")
    public List<File> INPUT;

    @Option(shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF to be written.")
    public File OUTPUT;

    @Option(shortName=StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, optional=true)
    public File SEQUENCE_DICTIONARY;

    private final Log log = Log.getInstance(SortVcf.class);

    public static void main(final String[] args) {
        new SortVcf().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        final List<String> sampleList = new ArrayList<String>();

        for (final File input : INPUT) IOUtil.assertFileIsReadable(input);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (SEQUENCE_DICTIONARY != null) IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);

        SAMSequenceDictionary samSequenceDictionary = null;
        if (SEQUENCE_DICTIONARY != null) samSequenceDictionary = SamReaderFactory.makeDefault().open(SEQUENCE_DICTIONARY).getFileHeader().getSequenceDictionary();

        // Gather up a file reader and file header for each input file. Check for sequence dictionary compatibility along the way.
        final List<VCFFileReader> inputReaders = new ArrayList<VCFFileReader>();
        final List<VCFHeader> inputHeaders = new ArrayList<VCFHeader>();
        collectFileReadersAndHeaders(sampleList, samSequenceDictionary, inputReaders, inputHeaders);

        // Create the merged output header from the input headers
        final VCFHeader outputHeader = new VCFHeader(VCFUtils.smartMergeHeaders(inputHeaders, false), sampleList);

        // Write to the sorting collection
        final SortingCollection<VariantContext> sortedOutput = sortInputs(inputReaders, outputHeader);

        // Output to the final file
        writeSortedOutput(outputHeader, sortedOutput);

        return 0;
    }

    private void collectFileReadersAndHeaders(final List<String> sampleList, SAMSequenceDictionary samSequenceDictionary, final List<VCFFileReader> inputReaders, final List<VCFHeader> inputHeaders) {
        for (final File input : INPUT) {
            final VCFFileReader in = new VCFFileReader(input, false);
            final VCFHeader header = in.getFileHeader();
            final SAMSequenceDictionary dict = in.getFileHeader().getSequenceDictionary();
            if (dict.isEmpty()) {
                if (null == samSequenceDictionary) {
                    throw new PicardException("Please specify SEQUENCE_DICTIONARY. Sequence dictionary was empty for the VCF: " + input.getAbsolutePath());
                }
                header.setSequenceDictionary(samSequenceDictionary);
            } else {
                if (null == samSequenceDictionary) samSequenceDictionary = dict;
                else samSequenceDictionary.assertSameDictionary(dict);
            }
            if (sampleList.isEmpty()) {
                sampleList.addAll(header.getSampleNamesInOrder());
            } else {
                if ( ! sampleList.equals(header.getSampleNamesInOrder())) {
                    throw new IllegalArgumentException("Input file " + input.getAbsolutePath() + " has sample entries that don't match the other files.");
                }
            }
            inputReaders.add(in);
            inputHeaders.add(header);
        }
    }

    private SortingCollection<VariantContext> sortInputs(final List<VCFFileReader> readers, final VCFHeader outputHeader) {ProgressLogger progress = new ProgressLogger(log, 25000, "read", "records");
        final SortingCollection<VariantContext> sorter =
                SortingCollection.newInstance(
                        VariantContext.class,
                        new VCFRecordCodec(outputHeader),
                        new VariantContextComparator(outputHeader.getSequenceDictionary()),
                        MAX_RECORDS_IN_RAM,
                        TMP_DIR);
        for (final VCFFileReader reader : readers) {
            for (final VariantContext variantContext : reader) {
                sorter.add(variantContext);
                progress.record(variantContext.getChr(), variantContext.getStart());
            }
            reader.close();
        }
        return sorter;
    }

    private void writeSortedOutput(final VCFHeader outputHeader, final SortingCollection<VariantContext> sortedOutput) {
        final ProgressLogger progress;
        final EnumSet<Options> options = CREATE_INDEX ? EnumSet.of(Options.INDEX_ON_THE_FLY) : EnumSet.noneOf(Options.class);
        progress = new ProgressLogger(log, 25000, "wrote", "records");
        final VariantContextWriter out = new VariantContextWriterBuilder().
                setReferenceDictionary(outputHeader.getSequenceDictionary()).
                setOptions(options).
                setOutputFile(OUTPUT).build();
        out.writeHeader(outputHeader);
        for (final VariantContext variantContext : sortedOutput) {
            out.add(variantContext);
            progress.record(variantContext.getChr(), variantContext.getStart());
        }
        out.close();
    }
}
