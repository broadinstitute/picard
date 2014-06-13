package picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.OneLineUsage;
import picard.cmdline.Option;
import picard.cmdline.ProviderFor;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.Usage;

import java.io.File;
import java.util.List;

/**
 * A tool to add comments to a BAM file header. Effectively copies the BAM file except for the addition of the @CO records
 * in the header. This tool does not support SAM files.
 *
 * @author jgentry
 */
@ProviderFor(CommandLineProgram.class)
public class AddCommentsToBam extends CommandLineProgram {
    @Usage public final String USAGE = "Adds one or more comments to the header of a specified BAM file. Copies the file with the " +
            "modified header to a specified output file. Note that a block copying method is used to ensure efficient transfer to the " +
            "output file. SAM files are not supported";

    @OneLineUsage
    public String ONE_LINE_USAGE = "Adds comments to the header of a BAM file";

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file to add a comment to the header")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file to write results")
    public File OUTPUT;

    @Option(shortName="C", doc="Comments to add to the BAM file")
    public List<String> COMMENT;

    public static void main(final String[] args) { new AddCommentsToBam().instanceMainWithExit(args); }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (INPUT.getAbsolutePath().endsWith(".sam")) {
            throw new PicardException("SAM files are not supported");
        }

        final SAMFileHeader samFileHeader = new SAMFileReader(INPUT).getFileHeader();
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
