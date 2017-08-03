package picard.sam.util;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;

/**
 * Abstract class that holds parameters and methods common to classes that want to add PG tags to reads in SAM/BAM files
 *
 */
public abstract class ReadOutputCommandLineProgram extends CommandLineProgram {

    @Option(doc = "Add PG tag to each read in a SAM or BAM", common = true)
    public Boolean ADD_PG_TAG_TO_READS = false;
}
