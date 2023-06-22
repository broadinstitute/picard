package picard.vcf;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalListWriter;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

/**
 * Converts a VCF or BCF file to a Picard Interval List.
 *
 * <p>This tool creates a Picard Interval List from a VCF or BCF. It is important
 * that the file extension is included as the file format is determined by the file
 * extension. Variants that were filtered can be included in the output interval list by
 * setting INCLUDE_FILTERED to true. </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li> A BCF or VCF input file  </li>
 *     <li> Boolean if variants that were filtered should be included in the output interval list </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 * A Picard Interval List
 * </ul>
 *
 * <p>
 * <h4>Usage example:</h4>
 * <pre>
 *     java -jar picard.jar VcfToIntervalList \
 *          I=input_variants.vcf \
 *          O=output.interval_list
 * </pre>
 * </p>
 *
 * @author ggrant@broadinstitute.org
 */
@CommandLineProgramProperties(
        summary = VcfToIntervalList.USAGE_DETAILS,
        oneLineSummary = VcfToIntervalList.USAGE_SUMMARY,
        programGroup = VariantManipulationProgramGroup.class)
@DocumentedFeature
public class VcfToIntervalList extends CommandLineProgram {

    static final String USAGE_SUMMARY = "Converts a VCF or BCF file to a Picard Interval List";
    static final String USAGE_DETAILS = "This tool creates a Picard Interval List from a VCF or BCF. " +
            "It is important that the file extension is included as the file format is determined by the file" +
            "extension. Variants that were filtered can be included in the output interval list by setting" +
            "INCLUDE_FILTERED to true." +
            "<p>"+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar VcfToIntervalList <br />" +
            "      I=sample.vcf <br />" +
            "      O=sample.interval_list <br />"+
            "</pre>";

    public static final String INCLUDE_FILTERED_SHORT_NAME = "IF";

    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfToIntervalList.class);

    @Argument(doc="The BCF or VCF input file. The file format is determined by file extension.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List.")
    public File OUTPUT;


    public enum VARIANT_ID_TYPES {
        CONCAT_ALL,
        USE_FIRST
    }

    @Argument(doc="Controls the naming of the resulting intervals. When set to CONCAT_ALL (the default), each resulting " +
            "interval will be named the concatenation of the variant ID fields (if present), or 'interval-<number>' " +
            "(if not) with a pipe '|' separator. If set to USE_FIRST, only the first name will be used.")
    public static VARIANT_ID_TYPES VARIANT_ID_METHOD = VARIANT_ID_TYPES.CONCAT_ALL;

    @Argument(shortName = INCLUDE_FILTERED_SHORT_NAME,
            doc = "Include variants that were filtered in the output interval list.",
            optional = true)
    public boolean INCLUDE_FILTERED = false;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        final boolean concatenate_ids = (VARIANT_ID_METHOD == VARIANT_ID_TYPES.CONCAT_ALL);

        try (VCFFileReader vcfReader = new VCFFileReader(INPUT.toPath(), false)) {
            final Iterator<Interval> samFileIterator = VCFFileReader.toIntervals(vcfReader, INCLUDE_FILTERED);
            try (IntervalListWriter writer = new IntervalListWriter(OUTPUT.toPath(), new SAMFileHeader(vcfReader.getFileHeader().getSequenceDictionary()))) {
                final IntervalList.IntervalMergerIterator mergingIterator =
                        new IntervalList.IntervalMergerIterator(samFileIterator, true, false, concatenate_ids);
                for (final Interval interval : new IterableAdapter<>(mergingIterator)){
                    writer.write(interval);
                }
            } catch (IOException e) {
                if(!OUTPUT.renameTo(new File(OUTPUT.getAbsolutePath() + ".incomplete"))){
                    OUTPUT.delete();
                }
                throw new PicardException("Trouble writing IntervalList. Renamed <OUTPUT> to <OUTPUT>.incomplete to avoid misuse.", e);
            }
        }
        return 0;
    }
}

