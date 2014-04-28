package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.BamFileIoUtils;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

import java.io.File;
import java.util.List;

/**
 * A tool to add comments to a BAM file header. Effectively copies the BAM file except for the addition of the @CO records
 * in the header. This tool does not support SAM files.
 *
 * @author jgentry
 */
public class AddCommentsToBam extends CommandLineProgram {
    @Usage public final String USAGE = "Adds one or more comments to the header of a specified BAM file. Copies the file with the " +
            "modified header to a specified output file. Note that a block copying method is used to ensure efficient transfer to the " +
            "output file. SAM files are not supported";
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input BAM file to add a comment to the header")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file to write results")
    public File OUTPUT;

    @Option(shortName="C", doc="Comments to add to the BAM file")
    public List<String> COMMENT;

    public static void main(final String[] args) { new AddCommentsToBam().instanceMainWithExit(args); }

    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

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
