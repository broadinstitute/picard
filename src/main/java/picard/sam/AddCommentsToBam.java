package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.util.List;

/**
 * A tool to add comments to a BAM file header. Effectively copies the BAM file except for the addition of the @CO records
 * in the header. This tool does not support SAM files. A block copying method is used to ensure efficient transfer to the output file.
 *
 * @author jgentry
 */
@CommandLineProgramProperties(
        summary =  AddCommentsToBam.USAGE_SUMMARY + AddCommentsToBam.USAGE_DETAILS,
        oneLineSummary = AddCommentsToBam.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class AddCommentsToBam extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Adds comments to the header of a BAM file.";
    static final String USAGE_DETAILS = "This tool makes a copy of the input bam file, with a modified header that includes the comments " +
            "specified at the command line (prefixed by @CO). Use double quotes to wrap comments that include whitespace or special characters. <br /><br />" +
            "Note that this tool cannot be run on SAM files.<br />"          +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar AddCommentsToBam \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=modified_bam.bam \\<br />" +
            "      C=comment_1 \\<br />" +
            "      C=\"comment 2\"" +
            "</pre>" +
            "" +
            "<hr />";
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input BAM file to add a comment to the header")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output BAM file to write results")
    public File OUTPUT;

    @Argument(shortName = "C", doc = "Comments to add to the BAM file")
    public List<String> COMMENT;

    public static void main(final String[] args) { new AddCommentsToBam().instanceMainWithExit(args); }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (INPUT.getAbsolutePath().endsWith(".sam")) {
            throw new PicardException("SAM files are not supported");
        }

        final SAMFileHeader samFileHeader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).getFileHeader(INPUT);
        for (final String comment : COMMENT) {
            if (comment.contains("\n")) {
                throw new PicardException("Comments can not contain a new line");
            }
            samFileHeader.addComment(comment);
        }

        BamFileIoUtils.reheaderBamFile(samFileHeader, INPUT, OUTPUT, CREATE_MD5_FILE, CREATE_INDEX);

        return 0;
    }
}
