package picard.sam.markduplicates;

import org.broadinstitute.barclay.argparser.Argument;

public class MarkDuplicatesForFlowArgumentCollection {

    @Argument(doc = "enable parameters and behavior specific to flow based reads.", optional = true)
    public boolean FLOW_MODE = false;

    @Argument(doc = "Use specific quality summing strategy for flow based reads. The strategy ensures that the same " +
            "(and correct) quality value is used for all bases of the same homopolymer.", optional = true)
    public boolean FLOW_QUALITY_SUM_STRATEGY = false;

    @Argument(doc = "Make the end location of single end read be significant when considering duplicates, " +
            "in addition to the start location, which is always significant (i.e. require single-ended reads to start and" +
            "end on the same position to be considered duplicate) " +
            "(for this argument, \"read end\" means 3' end).", optional = true)
    public boolean USE_END_IN_UNPAIRED_READS = false;

    @Argument(doc = "Use position of the clipping as the end position, when considering duplicates (or use the unclipped end position) " +
            "(for this argument, \"read end\" means 3' end).", optional = true)
    public boolean USE_UNPAIRED_CLIPPED_END = false;

    @Argument(doc = "Maximal difference of the read end position that counted as equal. Useful for flow based " +
            "reads where the end position might vary due to sequencing errors. " +
            "(for this argument, \"read end\" means 3' end)", optional = true)
    public int UNPAIRED_END_UNCERTAINTY = 0;

    @Argument(doc = "Skip first N flows, starting from the read's start, when considering duplicates. Useful for flow based reads where sometimes there " +
            "is noise in the first flows " +
            "(for this argument, \"read start\" means 5' end).", optional = true)
    public int FLOW_SKIP_FIRST_N_FLOWS = 0;

    @Argument(doc = "Treat position of read trimming based on quality as the known end (relevant for flow based reads). Default false - if the read " +
            "is trimmed on quality its end is not defined and the read is duplicate of any read starting at the same place.", optional = true)
    public boolean FLOW_Q_IS_KNOWN_END = false;

    @Argument(doc = "Threshold for considering a quality value high enough to be included when calculating FLOW_QUALITY_SUM_STRATEGY calculation.", optional = true)
    public int FLOW_EFFECTIVE_QUALITY_THRESHOLD = 15;

}
