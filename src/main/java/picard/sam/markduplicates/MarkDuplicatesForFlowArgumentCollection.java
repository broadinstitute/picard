package picard.sam.markduplicates;

import org.broadinstitute.barclay.argparser.Argument;

public class MarkDuplicatesForFlowArgumentCollection {

    public enum FLOW_DUPLICATE_SELECTION_STRATEGY {
        FLOW_QUALITY_SUM_STRATEGY,
        FLOW_END_QUALITY_STRATEGY
    }
    @Argument(doc = "enable parameters and behavior specific to flow based reads.", optional = true)
    public boolean FLOW_MODE = false;

    @Argument(doc = "Use specific quality summing strategy for flow based reads. Two strategies are available: " + 
          "FLOW_QUALITY_SUM_STRATEG: The selects the read with the best total homopolymer quality." + 
            " FLOW_END_QUALITY_STRATEGY: The strategy selects the read with the best homopolymer quality close to the end (10 bases) of the read. " +
            " The latter strategy is recommended for samples with high duplication rate ", optional = true)
    public FLOW_DUPLICATE_SELECTION_STRATEGY FLOW_DUP_STRATEGY = FLOW_DUPLICATE_SELECTION_STRATEGY.FLOW_QUALITY_SUM_STRATEGY;

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

    @Argument(doc = "Maximal difference of the read start position that counted as equal. Useful for flow based " +
            "reads where the end position might vary due to sequencing errors. " +
            "(for this argument, \"read start\" means 5' end in the direction of sequencing)", optional = true)
    public int UNPAIRED_START_UNCERTAINTY = 0;

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
