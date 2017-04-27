package picard.analysis.directed;

import picard.metrics.MultilevelMetrics;

/** Metrics class for the analysis of reads obtained from targeted pcr experiments e.g. the TruSeq Custom Amplicon
 * (TSCA) kit (Illumina).  */
public class TargetedPcrMetrics extends MultilevelMetrics {

    /**  The name of the amplicon set used in this metrics collection run */
    public String CUSTOM_AMPLICON_SET;

    /** The number of bases in the reference genome used for alignment */
    public long GENOME_SIZE;

    /** The number of unique bases covered by the intervals of all amplicons in the amplicon set */
    public long AMPLICON_TERRITORY;

    /** The number of unique bases covered by the intervals of all targets that should be covered */
    public long TARGET_TERRITORY;

    /** The total number of reads in the SAM or BAM file examined */
    public long TOTAL_READS;

    /** The total number of reads passing filter (PF), where the filter(s) can be platform/vendor quality controls*/
    public long PF_READS;

    /** The total number of bases within the PF_READS of the SAM or BAM file to be examined */
    public long PF_BASES;

    /** The number of PF_READS that were not marked as sample or optical duplicates. */
    public long PF_UNIQUE_READS;

    /** The fraction of reads passing filter, PF_READS/TOTAL_READS. */
    public double PCT_PF_READS;

    /** The fraction of TOTAL_READS that are unique, PF, and are not duplicates, PF_UNIQUE_READS/TOTAL_READS */
    public double PCT_PF_UQ_READS;

    /** The total number of PF_UNIQUE_READS that align to the reference genome with mapping scores > 0 */
    public long PF_UQ_READS_ALIGNED;

    /** Tracks the number of PF read pairs (used to calculate library size) */
    public long PF_SELECTED_PAIRS;

    /** Tracks the number of unique, PF, read pairs, observed (used to calculate library size) */
    public long PF_SELECTED_UNIQUE_PAIRS;

    /** Fraction of PF_READS that are unique and align to the reference genome, PF_UQ_READS_ALIGNED/PF_READS */
    public double PCT_PF_UQ_READS_ALIGNED;

    /** The number of bases from PF_READS that align to the reference genome with mapping score > 0 */
    public long PF_BASES_ALIGNED;

    /** The number of bases from PF_UNIQUE_READS that align to the reference genome and have a mapping score > 0 */
    public long PF_UQ_BASES_ALIGNED;

    /** The number of PF_BASES_ALIGNED that mapped to an amplified region of the genome. */
    public long ON_AMPLICON_BASES;

    /** The number of PF_BASES_ALIGNED that mapped to within a fixed interval of an amplified region, but not on a
     * baited region. */
    public long NEAR_AMPLICON_BASES;

    /** The number of PF_BASES_ALIGNED that mapped neither on or near an amplicon. */
    public long OFF_AMPLICON_BASES;

    /** The number of PF_BASES_ALIGNED that mapped to a targeted region of the genome. */
    public long ON_TARGET_BASES;

    /** The number of bases from PF_SELECTED_UNIQUE_PAIRS that mapped to a targeted region of the genome. */
    public long ON_TARGET_FROM_PAIR_BASES;

    /** The fraction of PF_BASES_ALIGNED that mapped to or near an amplicon, (ON_AMPLICON_BASES +
     * NEAR_AMPLICON_BASES)/PF_BASES_ALIGNED. */
    public double PCT_AMPLIFIED_BASES;

    /** The fraction of PF_BASES_ALIGNED that mapped neither onto or near an amplicon,
     * OFF_AMPLICON_BASES/PF_BASES_ALIGNED */
    public double PCT_OFF_AMPLICON;

    /**
     * The fraction of bases mapping to regions on or near amplicons, which mapped directly to but not near
     * amplicons, ON_AMPLICON_BASES/(NEAR_AMPLICON_BASES + ON_AMPLICON_BASES)
     * */
    public double ON_AMPLICON_VS_SELECTED;

    /** The mean read coverage of all amplicon regions in the experiment. */
    public double MEAN_AMPLICON_COVERAGE;

    /** The mean read coverage of all target regions in an experiment. */
    public double MEAN_TARGET_COVERAGE;

    /** The median coverage of reads that mapped to target regions of an experiment. */
    public double MEDIAN_TARGET_COVERAGE;

    /** The maximum coverage of reads that mapped to target regions of an experiment. */
    public long MAX_TARGET_COVERAGE;

    /** The fold by which the amplicon region has been amplified above genomic background. */
    public double FOLD_ENRICHMENT;

    /** The fraction of targets that did not reach coverage=1 over any base. */
    public double ZERO_CVG_TARGETS_PCT;

    /** The fraction of aligned bases that were filtered out because they were in reads marked as duplicates. */
    public double PCT_EXC_DUPE;

    /** The fraction of aligned bases that were filtered out because they were in reads with low mapping quality. */
    public double PCT_EXC_MAPQ;

    /** The fraction of aligned bases that were filtered out because they were of low base quality. */
    public double PCT_EXC_BASEQ;

    /** The fraction of aligned bases that were filtered out because they were the second observation from an
     * insert with overlapping reads. */
    public double PCT_EXC_OVERLAP;

    /** The fraction of bases that were filtered out because they did not map to a base within a target region. */
    public double PCT_EXC_OFF_TARGET;
    /**
     * The fold over-coverage necessary to raise 80% of bases in "non-zero-cvg" targets to
     * the mean coverage level in those targets.
     */
    public double FOLD_80_BASE_PENALTY;

    /** The fraction of all target bases achieving 1X or greater coverage. */
    public double PCT_TARGET_BASES_1X;
    /** The fraction of all target bases achieving 2X or greater coverage depth. */
    public double PCT_TARGET_BASES_2X;
    /** The fraction of all target bases achieving 10X or greater coverage depth. */
    public double PCT_TARGET_BASES_10X;
    /** The fraction of all target bases achieving 20X or greater coverage depth. */
    public double PCT_TARGET_BASES_20X;
	/** The fraction of all target bases achieving 30X or greater coverage depth. */
	public double PCT_TARGET_BASES_30X;

    /**
     * A measure of how regions with low GC content (<= 50%), are undercovered relative to mean coverage.
     * After binning the GC content [0..50], we calculate a = fraction of target territory, and b = fraction of
     * aligned reads aligned to these targets for each bin.  AT DROPOUT is then abs(sum(a-b when a-b < 0)).
     * For example, if the AT_DROPOUT value is 5% this implies that 5% of total reads that
     * should have mapped to GC<=50% regions, mapped elsewhere.
     * */
    public double AT_DROPOUT;

    /**
     * A measure of how regions of high GC content (>= 50% GC) are undercovered relative to the mean coverage
     * value. For each GC bin [50..100], we calculate a = % of target territory, and b = % of aligned reads aligned
     * to these targets.  GC DROPOUT is then abs(sum(a-b when a-b < 0)).  For example, if the value is 5%, this
     * implies that 5% of total reads that should have mapped to GC>=50% regions, mapped elsewhere.
     * */
    public double GC_DROPOUT;

    /** The theoretical HET SNP sensitivity. */
    public double HET_SNP_SENSITIVITY;

    /** The Q Score of the theoretical HET SNP sensitivity. */
    public double HET_SNP_Q;
}
