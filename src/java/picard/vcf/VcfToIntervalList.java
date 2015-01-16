package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;

import java.io.File;

/**
 * Creates an interval list from a VCF
 *
 * @author ggrant@broadinstitute.org
 */

@CommandLineProgramProperties(
        usage = "Converts a VCF or BCF file to a Picard Interval List.",
        usageShort = "Converts a VCF or BCF file to a Picard Interval List.",
        programGroup = VcfOrBcf.class
)
public class VcfToIntervalList extends CommandLineProgram {
    // The following attributes define the command-line arguments
    public static final Log LOG = Log.getInstance(VcfToIntervalList.class);

    @Option(doc="The BCF or VCF input file. The file format is determined by file extension.", shortName= StandardOptionDefinitions.INPUT_SHORT_NAME)
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output Picard Interval List")
    public File OUTPUT;

    public static void main(final String[] argv) {
        new VcfToIntervalList().instanceMainWithExit(argv);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervalList = VCFFileReader.fromVcf(INPUT);
        // Sort and write the output
        intervalList.uniqued().write(OUTPUT);
        return 0;
    }
}