package picard.sam;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

@CommandLineProgramProperties(
        summary = CalculateReadGroupChecksum.USAGE_SUMMARY + CalculateReadGroupChecksum.USAGE_DETAILS,
        oneLineSummary = CalculateReadGroupChecksum.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class CalculateReadGroupChecksum extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Creates a hash code based on the read groups (RG).  ";
    static final String USAGE_DETAILS = "This tool creates a hash code based on identifying information in the read groups " +
            "(RG) of a \".BAM\" or \"SAM\" file header.  Addition or removal of RGs changes the hash code, enabling the user to " +
            "quickly determine if changes have been made to the read group information. " +
            "<br />" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar CalculateReadGroupChecksum \\<br />" +
            "      I=input.bam" +
            "</pre>"  +
            "Please see the AddOrReplaceReadGroups tool documentation for information regarding the addition, subtraction, or merging of read groups."     +
            "<hr />";
    private static final String OUTPUT_FILE_EXTENSION = ".read_group_md5";

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file. ")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The file to which the hash code should be written.", optional = true)
    public File OUTPUT;

    public static void main(final String[] args) {
        new CalculateReadGroupChecksum().instanceMainWithExit(args);
    }

    /**
     * Creates a file name (not including the path) for an RG MD5 file based on the name of the input file.
     */
    public static String getOutputFileName(final File inputFile) {
        return inputFile.getName() + OUTPUT_FILE_EXTENSION;
    }

    @Override
    protected int doWork() {
        final File output =
                OUTPUT == null
                        ? new File(INPUT.getParentFile(), getOutputFileName(INPUT))
                        : OUTPUT;

        IOUtil.assertFileIsWritable(output);
        final String hashText = SAMUtils.calculateReadGroupRecordChecksum(INPUT, REFERENCE_SEQUENCE);

        try {
            final FileWriter outputWriter = new FileWriter(output);
            outputWriter.write(hashText);
            outputWriter.close();
        } catch (final IOException ioe) {
            throw new PicardException(
                    "Could not write the computed hash (" + hashText + ") to the output file: " + ioe.getMessage(), ioe);
        }
        return 0;
    }
}