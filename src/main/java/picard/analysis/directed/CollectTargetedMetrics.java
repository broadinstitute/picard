package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import picard.analysis.MetricAccumulationLevel;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.metrics.MultilevelMetrics;

import java.io.File;
import java.util.List;
import java.util.Set;

import static picard.cmdline.StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME;

/**
 * <p>Both CollectTargetedPCRMetrics and CollectHsSelection share virtually identical program structures except
 * for the name of their targeting mechanisms (e.g. bait set or amplicon set).  The shared behavior of these programs
 * is encapsulated in CollectTargetedMetrics which is then subclassed by CollectHsMetrics and CollectTargetedPcrMetrics.
 * <p/>
 * <p>This program verifies the input parameters to TargetMetricsCollector and converts all files to
 * the format desired by TargetMetricsCollector.  Then it instantiates a TargetMetricsCollector and
 * collects metric information for all reads in the INPUT sam file.</p>
 */
public abstract class CollectTargetedMetrics<METRIC extends MultilevelMetrics, COLLECTOR extends TargetMetricsCollector<METRIC>> extends CommandLineProgram {

    private final Log log = Log.getInstance(getClass());

    protected abstract IntervalList getProbeIntervals();

    protected abstract String getProbeSetName();

    /**
     * A factory method for the TargetMetricsCollector to use this time.  Examples of TargetMetricsCollector:
     * (TargetedPcrMetricsCollector, HsMetricsCalculator)
     *
     * @return A TargetMetricsCollector to which we will pass SAMRecords
     */
    protected abstract COLLECTOR makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                               final List<SAMReadGroupRecord> samRgRecords,
                                               final ReferenceSequenceFile refFile,
                                               final File perTargetCoverage,
                                               final File perBaseCoverage,
                                               final IntervalList targetIntervals,
                                               final IntervalList probeIntervals,
                                               final String probeSetName,
                                               final int nearProbeDistance);


    @Argument(shortName = "TI", doc = "An interval list file that contains the locations of the targets.", minElements=1)
    public List<File> TARGET_INTERVALS;

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "An aligned SAM or BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output file to write the metrics to.")
    public File OUTPUT;

    @Argument(shortName = "LEVEL", doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(optional = true, doc = "An optional file to output per target coverage information to.")
    public File PER_TARGET_COVERAGE;

    @Argument(optional = true, doc = "An optional file to output per base coverage information to. The per-base file contains " +
            "one line per target base and can grow very large. It is not recommended for use with large target sets.")
    public File PER_BASE_COVERAGE;

    @Argument(optional = true, doc= "The maximum distance between a read and the nearest probe/bait/amplicon for the read to be " +
            "considered 'near probe' and included in percent selected.")
    public int NEAR_DISTANCE = TargetedPcrMetricsCollector.NEAR_PROBE_DISTANCE_DEFAULT;

    @Argument(shortName = MINIMUM_MAPPING_QUALITY_SHORT_NAME, doc = "Minimum mapping quality for a read to contribute coverage.")
    public int MINIMUM_MAPPING_QUALITY = 1;

    @Argument(shortName = "Q", doc = "Minimum base quality for a base to contribute coverage.")
    public int MINIMUM_BASE_QUALITY = 0;

    @Argument(doc = "True if we are to clip overlapping reads, false otherwise.", optional=true)
    public boolean CLIP_OVERLAPPING_READS = false;

    @Argument(shortName = "covMax", doc = "Parameter to set a max coverage limit for Theoretical Sensitivity calculations. Default is 200.", optional = true)
    public int COVERAGE_CAP = 200;

    @Argument(doc="Sample Size used for Theoretical Het Sensitivity sampling. Default is 10000.", optional = true)
    public int SAMPLE_SIZE=10000;

    /**
     * Asserts that files are readable and writable and then fires off an
     * HsMetricsCalculator instance to do the real work.
     */
    protected int doWork() {
        for (final File targetInterval : TARGET_INTERVALS) IOUtil.assertFileIsReadable(targetInterval);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (PER_TARGET_COVERAGE != null) IOUtil.assertFileIsWritable(PER_TARGET_COVERAGE);

        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final IntervalList targetIntervals = IntervalList.fromFiles(TARGET_INTERVALS);

        // Validate that the targets and baits have the same references as the reads file
        SequenceUtil.assertSequenceDictionariesEqual(
                reader.getFileHeader().getSequenceDictionary(),
                targetIntervals.getHeader().getSequenceDictionary());
        SequenceUtil.assertSequenceDictionariesEqual(
                reader.getFileHeader().getSequenceDictionary(),
                getProbeIntervals().getHeader().getSequenceDictionary()
        );

        ReferenceSequenceFile ref = null;
        if (REFERENCE_SEQUENCE != null) {
            IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
            SequenceUtil.assertSequenceDictionariesEqual(
                    reader.getFileHeader().getSequenceDictionary(), ref.getSequenceDictionary(),
                    INPUT, REFERENCE_SEQUENCE
            );
        }

        final COLLECTOR collector = makeCollector(
                METRIC_ACCUMULATION_LEVEL,
                reader.getFileHeader().getReadGroups(),
                ref,
                PER_TARGET_COVERAGE,
                PER_BASE_COVERAGE,
                targetIntervals,
                getProbeIntervals(),
                getProbeSetName(),
                NEAR_DISTANCE
        );

        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord record : reader) {
            collector.acceptRecord(record, null);
            progress.record(record);
        }

        // Write the output file
        final MetricsFile<METRIC, Integer> metrics = getMetricsFile();
        collector.finish();

        collector.addAllLevelsToFile(metrics);

        metrics.write(OUTPUT);

        CloserUtil.close(reader);
        return 0;
    }

    /** Renders a probe name from the provided file, returning {@link java.io.File#getName()} with all extensions stripped. */
    static String renderProbeNameFromFile(final File probeIntervalFile) {
        final String name = probeIntervalFile.getName();
        final int firstPeriodIndex = name.indexOf('.');
        if (firstPeriodIndex == -1) {
            return name;
        } else {
            return name.substring(0, firstPeriodIndex);
        }
    }

    protected String[] customCommandLineValidation() {
        if (PER_TARGET_COVERAGE != null && (METRIC_ACCUMULATION_LEVEL.size() != 1 ||
                METRIC_ACCUMULATION_LEVEL.iterator().next() != MetricAccumulationLevel.ALL_READS)) {
            return new String[]{"PER_TARGET_COVERAGE can be specified only when METRIC_ACCUMULATION_LEVEL is set " +
                    "to ALL_READS."};
        }

        if (PER_TARGET_COVERAGE != null && REFERENCE_SEQUENCE == null) {
            return new String[]{"Must supply REFERENCE_SEQUENCE when supplying PER_TARGET_COVERAGE"};
        }

        return super.customCommandLineValidation();
    }
}
