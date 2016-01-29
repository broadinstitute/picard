package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.IntervalList;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * Collect metric information for target pcr metrics runs.  See CollectTargetedMetrics and TargetPcrMetricsCollector for
 * more information
 */
@CommandLineProgramProperties(
        usage = CollectTargetedPcrMetrics.USAGE_SUMMARY + CollectTargetedPcrMetrics.USAGE_DETAILS,
        usageShort = CollectTargetedPcrMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectTargetedPcrMetrics extends CollectTargetedMetrics<TargetedPcrMetrics, TargetedPcrMetricsCollector> {
    static final String USAGE_SUMMARY = "Calculate PCR-related metrics from targeted sequencing data. ";
    static final String USAGE_DETAILS = "This tool calculates a set of PCR-related metrics from an aligned SAM or " +
            "BAM file containing targeted sequencing data. It is appropriate for data produced with multiple small-target technologies " +
            "including exome sequencing an custom amplicon panels such as the Illumina " +
            "<a href='http://www.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_truseq_custom_amplicon.pdf'>TruSeq Custom Amplicon (TSCA)</a> kit. <br /><br />" +
            "If a reference sequence is provided, AT/GC dropout metrics will be calculated and the PER_TARGET_COVERAGE  option can be " +
            "used to output GC content and mean coverage information for each target. The AT/GC dropout metrics indicate the degree of " +
            "inadequate coverage of a particular region based on its AT or GC content. The PER_TARGET_COVERAGE option can be used to " +
            "output GC content and mean sequence depth information for every target interval. <br /><br />" +
            "Please note that coverage depth at each locus should not exceed a limit of java MAX_SHORT ~32K.  This is because " +
            "CollectTargetedPcrMetrics tool uses a short array to calculate coverage metrics." +
            "<h4>Usage Example</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectTargetedPcrMetrics \\<br /> " +
            "      I=input.bam \\<br /> " +
            "      O=pcr_metrics.txt \\<br /> " +
            "      R=reference_sequence.fasta \\<br /> " +
            "      AMPLICON_INTERVALS=amplicon.interval_list \\<br /> " +
            "      TARGET_INTERVALS=targets.interval_list " +
            "</pre>" +
            "For explanations of the output metrics, see " +
            "http://broadinstitute.github.io/picard/picard-metric-definitions.html#TargetedPcrMetrics" +
            "<hr />";
    @Option(shortName = "AI", doc = "An interval list file that contains the locations of the baits used.")
    public File AMPLICON_INTERVALS;

    @Option(shortName = "N", doc = "Custom amplicon set name. If not provided it is inferred from the filename of the AMPLICON_INTERVALS intervals.", optional = true)
    public String CUSTOM_AMPLICON_SET_NAME;

    /**
     * @return AMPLICON_INTERVALS
     */
    @Override
    protected IntervalList getProbeIntervals() {
        return IntervalList.fromFile(AMPLICON_INTERVALS);
    }

    /**
     * @return CUSTOM_AMPLICON_SET_NAME
     */
    @Override
    protected String getProbeSetName() {
        return CUSTOM_AMPLICON_SET_NAME != null ? CUSTOM_AMPLICON_SET_NAME : CollectTargetedMetrics.renderProbeNameFromFile(AMPLICON_INTERVALS);
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CollectTargetedPcrMetrics().instanceMain(argv));
    }

    @Override
    protected TargetedPcrMetricsCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                                        final List<SAMReadGroupRecord> samRgRecords,
                                                        final ReferenceSequenceFile refFile,
                                                        final File perTargetCoverage,
                                                        final File perBaseCoverage,
                                                        final IntervalList targetIntervals,
                                                        final IntervalList probeIntervals,
                                                        final String probeSetName,
                                                        final int nearProbeDistance) {
        return new TargetedPcrMetricsCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, perBaseCoverage, targetIntervals, probeIntervals, probeSetName, nearProbeDistance,
                MINIMUM_MAPPING_QUALITY, MINIMUM_BASE_QUALITY, CLIP_OVERLAPPING_READS, true, COVERAGE_CAP, SAMPLE_SIZE);
    }
}
