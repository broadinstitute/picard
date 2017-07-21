package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;

/**
 * Creates an interval list from a VCF
 *
 * @author ggrant@broadinstitute.org
 */

@CommandLineProgramProperties(
        summary = "Converts a VCF or BCF file to a Picard Interval List.",
        oneLineSummary = "Converts a VCF or BCF file to a Picard Interval List.",
        programGroup = VcfOrBcf.class)
@DocumentedFeature
public class VcfToIntervalList extends CommandLineProgram {
    public static final String INCLUDE_FILTERED_SHORT_NAME = "IF";

    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfToIntervalList.class);

    @Argument(doc="The BCF or VCF input file. The file format is determined by file extension.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
    public File OUTPUT;

    @Argument(shortName = INCLUDE_FILTERED_SHORT_NAME,
            doc = "Include variants that were filtered in the output interval list.",
            optional = true)
    public boolean INCLUDE_FILTERED = false;

    public static void main(final String[] argv) {
        new VcfToIntervalList().instanceMainWithExit(argv);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervalList = VCFFileReader.fromVcf(INPUT, INCLUDE_FILTERED);

        // Sort and write the output
        intervalList.uniqued().write(OUTPUT);
        return 0;
    }
}