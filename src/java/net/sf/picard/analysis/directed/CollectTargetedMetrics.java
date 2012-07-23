package net.sf.picard.analysis.directed;

import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.util.CollectionUtil;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.util.StopWatch;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * Both CollectTargetedPCRMetrics and CalculateHybridSelection metrics share virtually identical program structures except
 * for the name of their targeting mechanisms (e.g. bait set or amplicon set).  The shared behavior of these programs
 * is encapsulated in CollectTargetedMetrics which is then subclassed by CalculateHsMetrics and CollectTargetedPcrMetrics.
 *
 * This program verifies the input parameters to TargetMetricsCollector and converts all files to
 * the format desired by TargetMetricsCollector.  Then it instantiates a TargetMetricsCollector and
 * collects metric information for all reads in the INPUT sam file.
 */
public abstract class CollectTargetedMetrics extends CommandLineProgram {

    private static final Log log = Log.getInstance(CalculateHsMetrics.class);

    /**
     * The interval file to be fed to TargetMetricsCollector
     * @return An interval file that denotes the intervals of the regions targeted by the probes for this run that is
     *         passed to the TargetMetricsCollector produced by makeCollector
     */
    protected abstract File getProbeIntervals();

    /**
     * @return The name of the probe set used in this run, getProbeIntervals().getName() is
     */
    protected abstract String getProbeSetName();

    /**
     *  A factory method for the TargetMetricsCollector to use this time.  Examples of TargetMetricsCollector:
     *  (TargetedPcrMetricsCollector, HsMetricsCalculator)
     *  @return A TargetMetricsCollector to which we will pass SAMRecords
     */
    protected abstract TargetMetricsCollector makeCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                                             final List<SAMReadGroupRecord> samRgRecords,
                                                             final ReferenceSequenceFile refFile,
                                                             final File perTargetCoverage,
                                                             final File targetIntervals,
                                                             final File probeIntervals,
                                                             final String probeSetName);


    @Option(shortName="TI", doc="An interval list file that contains the locations of the targets.")
    public File TARGET_INTERVALS;

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="An aligned SAM or BAM file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file to write the metrics to.")
    public File OUTPUT;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, optional=true, doc="The reference sequence aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(optional=true, doc="An optional file to output per target coverage information to.")
    public File PER_TARGET_COVERAGE;

    /**
     * Asserts that files are readable and writable and then fires off an
     * HsMetricsCalculator instance to do the real work.
     */
    protected int doWork() {
        IoUtil.assertFileIsReadable(getProbeIntervals());
        IoUtil.assertFileIsReadable(TARGET_INTERVALS);
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        if (PER_TARGET_COVERAGE != null) IoUtil.assertFileIsWritable(PER_TARGET_COVERAGE);

        final SAMFileReader samReader = new SAMFileReader(INPUT);

        final File probeIntervals = getProbeIntervals();

        // Validate that the targets and baits have the same references as the reads file
        SequenceUtil.assertSequenceDictionariesEqual(samReader.getFileHeader().getSequenceDictionary(),
                IntervalList.fromFile(TARGET_INTERVALS).getHeader().getSequenceDictionary(),
                INPUT, TARGET_INTERVALS);
        SequenceUtil.assertSequenceDictionariesEqual(samReader.getFileHeader().getSequenceDictionary(),
                IntervalList.fromFile(probeIntervals).getHeader().getSequenceDictionary(),
                INPUT, probeIntervals);

        ReferenceSequenceFile ref = null;
        if (REFERENCE_SEQUENCE != null) {
            IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
            SequenceUtil.assertSequenceDictionariesEqual(samReader.getFileHeader().getSequenceDictionary(), ref.getSequenceDictionary(),
                    INPUT, REFERENCE_SEQUENCE);
        }

        final TargetMetricsCollector collector = makeCollector(METRIC_ACCUMULATION_LEVEL, samReader.getFileHeader().getReadGroups(), ref,
                PER_TARGET_COVERAGE, TARGET_INTERVALS, probeIntervals, getProbeSetName());


        // Add each record to the requested collectors
        final Iterator<SAMRecord> records = samReader.iterator();
        final ProgressLogger progress = new ProgressLogger(log);

        while (records.hasNext()) {
            final SAMRecord sam = records.next();
            collector.acceptRecord(sam, null);
            progress.record(sam);
        }

        // Write the output file
        final MetricsFile<HsMetrics, Integer> metrics = getMetricsFile();
        collector.finish();

        collector.addAllLevelsToFile(metrics);

        metrics.write(OUTPUT);

        return 0;
    }

    protected String[] customCommandLineValidation() {
        if (PER_TARGET_COVERAGE != null && (METRIC_ACCUMULATION_LEVEL.size() != 1 ||
                METRIC_ACCUMULATION_LEVEL.iterator().next() != MetricAccumulationLevel.ALL_READS)) {
            return new String[] {"PER_TARGET_COVERAGE can be specified only when METRIC_ACCUMULATION_LEVEL is set " +
                    "to ALL_READS."};
        }

        if(PER_TARGET_COVERAGE != null && REFERENCE_SEQUENCE == null) {
            return new String[] {"Must supply REFERENCE_SEQUENCE when supplying PER_TARGET_COVERAGE"};
        }

        return super.customCommandLineValidation();
    }
}