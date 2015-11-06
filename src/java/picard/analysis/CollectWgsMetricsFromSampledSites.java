package picard.analysis;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.Metrics;
import picard.filter.CountingFilter;

import java.io.File;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of whole genome sequencing experiments,
 * but only at a set of sampled positions.
 * It is important that the sampled positions be chosen so that they are spread out at least further than a read's length apart;
 * otherwise, you run the risk of double-counting reads in the metrics.
 *
 * @author ebanks
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "whole genome sequencing experiments, but only at a set of sampled positions.  " +
                "It is important that the sampled positions be chosen so that they are spread out " +
                "at least further than a read's length apart; otherwise, you run the risk of double-counting " +
                "reads in the metrics.",
        usageShort = "Writes whole genome sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectWgsMetricsFromSampledSites extends CollectWgsMetrics {

    @Option(shortName = "INTERVALS", doc = "An interval list file that contains the locations of the positions to assess.", optional = false)
    public File INTERVALS = null;

    public static void main(final String[] args) {
        new CollectWgsMetricsFromSampledSites().instanceMainWithExit(args);
    }

    @Override
    protected SamLocusIterator getLocusIterator(final SamReader in) {
        IOUtil.assertFileIsReadable(INTERVALS);
        return new SamLocusIterator(in, IntervalList.fromFile(INTERVALS));
    }

    /**
     * By design we want to count just those bases at the positions we care about, not across the entire read.
     * Therefore, we call filter.getFilteredRecords() so that only the bases in the pileup at a given position
     * are included in the calculations (with filter.getFilteredBases() we would be including other bases in
     * the read too).
     */
    @Override
    protected long getBasesExcludedBy(final CountingFilter filter) {
        return filter.getFilteredRecords();
    }

    // rename the class so that in the metric file it is annotated differently.
    public static class SampledWgsMetrics extends WgsMetrics {}

    @Override
    protected WgsMetrics generateWgsMetrics() {
        return new SampledWgsMetrics();
    }
}

