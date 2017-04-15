package picard.sam;

import org.broadinstitute.barclay.argparser.Argument;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
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
        summary = CheckTerminatorBlock.USAGE,
        oneLineSummary = CheckTerminatorBlock.USAGE,
        programGroup = SamOrBam.class
)
public class CheckTerminatorBlock extends CommandLineProgram {
    static final String USAGE = "Asserts the provided gzip file's (e.g., BAM) last block is well-formed; RC 100 otherwise";
    
    @Argument(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The block compressed file to check.")
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