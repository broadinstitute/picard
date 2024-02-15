package picard.sam.util;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;

/**
 * Argument collection for SAM comparison
 */
public class SAMComparisonArgumentCollection {
    @Argument(doc = "Perform lenient checking of header.  In this mode, species, assembly, ur, m5, fields of sequence records, and pg fields in the header may all differ. " +
        "Sequence record length must also be the same.")
    public boolean LENIENT_HEADER;

    @Argument(doc = "Perform lenient checking of duplicate marks.  In this mode, will reduce the number of mismatches by allowing the choice of the representative read in each duplicate set" +
            " to differ between the input files, as long as the duplicate sets agree.")
    public boolean LENIENT_DUP;

    @Argument(doc = "Count reads which have mapping quality below LOW_MQ_THRESHOLD in both files but are mapped to different locations as matches.  By default we count such reads as mismatching.")
    public boolean LENIENT_LOW_MQ_ALIGNMENT;

    @Argument(doc = "Count reads for which no mapping quality is available (mapping quality value 255) in both files but are mapped to different locations as matches.  By default we count such reads as mismatching.")
    public boolean LENIENT_UNKNOWN_MQ_ALIGNMENT;

    @Argument(doc = "When running in LENIENT_LOW_MQ_ALIGNMENT mode, reads which have mapping quality below this value will be counted as matches. " +
            "if LENIENT_LOW_MQ_ALIGNMENT is false (default), then this argument has no effect.")
    public int LOW_MQ_THRESHOLD = 3;

    @Argument(doc = "If set to true, generate a histogram for mapping quality concordance between the two SAM files and write it to the output metrics file. " +
            "In this histogram, an entry of 10 at bin \"20,30\" means that 10 reads in the left file have mapping quality 20 while those same reads have mapping quality 30 in the right file. " +
            "The reads are associated based solely on read names and not the mapped position. Only primary alignments are included, but all duplicate reads are counted individually.",
            optional = true)
    public boolean COMPARE_MQ = false;
}
