package net.sf.picard.analysis.directed;

import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMReadGroupRecord;

import java.io.File;
import java.util.List;
import java.util.Set;

/**
 * Collect metric information for target pcr metrics runs.  See CollectTargetedMetrics and TargetPcrMetricsCollector for
 * more information
 */
public class CollectTargetedPcrMetrics extends CollectTargetedMetrics {

    @Usage
    public final String USAGE =
            "Calculates a set of metrics to Illumina Truseq Custom Amplicon sequencing from an aligned SAM" +
                    "or BAM file. If a reference sequence is provided, AT/GC dropout metrics will " +
                    "be calculated, and the PER_TARGET_COVERAGE option can be used to output GC and " +
                    "mean coverage information for every target.";
    @Option(shortName="AI", doc="An interval list file that contains the locations of the baits used.")
    public File AMPLICON_INTERVALS;

    @Option(shortName="N",  doc="Custom amplicon set name. If not provided it is inferred from the filename of the AMPLICON_INTERVALS intervals.", optional=true)
    public String CUSTOM_AMPLICON_SET_NAME;

    /**
     * @return AMPLICON_INTERVALS
     */
    @Override
    protected File getProbeIntervals() {
        return AMPLICON_INTERVALS;
    }

    /**
     * @return CUSTOM_AMPLICON_SET_NAME
     */
    @Override
    protected String getProbeSetName() {
        return CUSTOM_AMPLICON_SET_NAME;
    }

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CollectTargetedPcrMetrics().instanceMain(argv));
    }

    @Override
    protected TargetMetricsCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                                   final List<SAMReadGroupRecord> samRgRecords,
                                                   final ReferenceSequenceFile refFile,
                                                   final File perTargetCoverage,
                                                   final File targetIntervals,
                                                   final File probeIntervals,
                                                   final String probeSetName) {
        return new TargetedPcrMetricsCollector(accumulationLevels, samRgRecords, refFile, perTargetCoverage, targetIntervals, probeIntervals, probeSetName);
    }
}
