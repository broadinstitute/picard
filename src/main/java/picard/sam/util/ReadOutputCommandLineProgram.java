package picard.sam.util;

import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.CommandLineProgram;

/**
 * Abstract class that holds parameters and methods common to classes that want to add PG tags to reads in SAM/BAM files
 *
 */
public abstract class ReadOutputCommandLineProgram extends CommandLineProgram {

    @Argument(doc = "Add PG tag to each read in a SAM or BAM", common = true)
    public Boolean ADD_PG_TAG_TO_READS = true;
}
