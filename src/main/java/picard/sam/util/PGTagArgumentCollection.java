package picard.sam.util;

import org.broadinstitute.barclay.argparser.Argument;

/**
 * Argument Collection which holds parameters common to classes that want to add PG tags to reads in SAM/BAM files
 */
public class PGTagArgumentCollection {

    @Argument(doc = "Add PG tag to each read in a SAM or BAM", common = true)
    public Boolean ADD_PG_TAG_TO_READS = true;

}