package picard.sam;

import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedInputStream.FileTermination;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.io.IOException;

/**
 * Simple class to check the terminator block of a SAM file.
 */
@CommandLineProgramProperties(
        usage = "Checks a BAM file or other block-compressed GZIP'd file and determines " +
                "if the file is terminated with an empty block, a healthy non-empty block or if the ending of the " +
                "file is defective (implying a truncated or non-block compressed file).",
        usageShort = "Checks a BAM file for the terminating empty-block",
        programGroup = SamOrBam.class
)
public class CheckTerminatorBlock extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The block compressed file to check.")
    public File INPUT;

    public static void main(final String[] args) {
        new CheckTerminatorBlock().instanceMainWithExit(args);
    }

    @Override protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        try {
            final FileTermination term = BlockCompressedInputStream.checkTermination(INPUT);
            System.err.println(term.name());
            if (term == FileTermination.DEFECTIVE) {
                return 100;
            }
            else {
                return 0;
            }
        }
        catch (IOException ioe) {
            throw new PicardException("Exception reading terminator block of file: " + INPUT.getAbsolutePath());
        }
    }
}