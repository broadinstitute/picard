package picard.analysis;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.QuerySortedReadPairIteratorUtil;

import java.io.File;
import java.util.List;

/**
 * Computes a number of metrics that are useful for evaluating coverage and performance of sequencing experiments.
 *
 * @author ebanks
 */
@CommandLineProgramProperties(
        usage = "Computes a number of metrics that are useful for evaluating coverage and performance of " +
                "sequencing experiments.",
        usageShort = "Writes sequencing-related metrics for a SAM or BAM file",
        programGroup = Metrics.class
)
public class CollectWgsMetricsFromQuerySorted extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output metrics file.")
    public File OUTPUT;

    @Option(shortName = "USABLE_MQ", doc = "Minimum mapping quality for a read to contribute to usable coverage.", overridable = true)
    public int MINIMUM_USABLE_MAPPING_QUALITY = 20;

    @Option(shortName = "USABLE_Q", doc = "Minimum base quality for a base to contribute to usable coverage.", overridable = true)
    public int MINIMUM_USABLE_BASE_QUALITY = 20;

    @Option(shortName = "RAW_MQ", doc = "Minimum mapping quality for a read to contribute to raw coverage.", overridable = true)
    public int MINIMUM_RAW_MAPPING_QUALITY = 0;

    @Option(shortName = "RAW_Q", doc = "Minimum base quality for a base to contribute to raw coverage.", overridable = true)
    public int MINIMUM_RAW_BASE_QUALITY = 3;

    private final Log log = Log.getInstance(CollectWgsMetricsFromQuerySorted.class);

    /** the various metrics types */
    public enum FILTERING_STRINGENCY { RAW, USABLE }

    /** Metrics for evaluating the performance of whole genome sequencing experiments. */
    public static class QuerySortedSeqMetrics extends CollectWgsMetrics.WgsMetrics {
        /** Identifier for metrics type */
        public FILTERING_STRINGENCY TYPE;

        /** The total number of bases, before any filters are applied. */
        public long TOTAL_BASES = 0;
        /** The number of passing bases, after all filters are applied. */
        public long TOTAL_PASSING_BASES = 0;

        /** The number of read pairs, before all filters are applied. */
        public long TOTAL_READ_PAIRS = 0;
        /** The number of duplicate read pairs, before all filters are applied. */
        public long TOTAL_DUPE_PAIRS = 0;

        /** The number of read pairs with standard orientations from which to calculate mean insert size, after filters are applied. */
        public long TOTAL_ORIENTED_PAIRS = 0;
        /** The mean insert size, after filters are applied. */
        public double MEAN_INSERT_SIZE = 0.0;
    }

    /** A private class to track the intermediate values of certain metrics */
    private class IntermediateMetrics {
        final QuerySortedSeqMetrics metrics = new QuerySortedSeqMetrics();
        long basesExcludedByDupes = 0;
        long basesExcludedByMapq = 0;
        long basesExcludedByPairing = 0;
        long basesExcludedByBaseq = 0;
        long basesExcludedByOverlap = 0;
        double insertSizeSum = 0.0;
    }

    public static void main(final String[] args) {
        new CollectWgsMetricsFromQuerySorted().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        // progress tracker
        final ProgressLogger progress = new ProgressLogger(log, 50000000, "Processed", "read pairs");

        // the SAM reader
        final SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        final PeekableIterator<SAMRecord> iterator = new PeekableIterator<>(reader.iterator());

        // the metrics to keep track of
        final IntermediateMetrics usableMetrics = new IntermediateMetrics();
        usableMetrics.metrics.TYPE = FILTERING_STRINGENCY.USABLE;
        final IntermediateMetrics rawMetrics = new IntermediateMetrics();
        rawMetrics.metrics.TYPE = FILTERING_STRINGENCY.RAW;

        // Loop through all the loci by read pairs
        QuerySortedReadPairIteratorUtil.ReadPair pairToAnalyze = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        while (pairToAnalyze != null) {
            // calculate intermediate metrics
            calculateMetricsForRead(pairToAnalyze, usableMetrics, MINIMUM_USABLE_MAPPING_QUALITY, MINIMUM_USABLE_BASE_QUALITY);
            calculateMetricsForRead(pairToAnalyze, rawMetrics, MINIMUM_RAW_MAPPING_QUALITY, MINIMUM_RAW_BASE_QUALITY);

            // record progress
            progress.record(pairToAnalyze.read1);

            // iterate
            pairToAnalyze = QuerySortedReadPairIteratorUtil.getNextReadPair(iterator);
        }

        // finalize and write the metrics
        final long genomeTerritory = reader.getFileHeader().getSequenceDictionary().getReferenceLength();
        usableMetrics.metrics.GENOME_TERRITORY = genomeTerritory;
        finalizeMetrics(usableMetrics);
        rawMetrics.metrics.GENOME_TERRITORY = genomeTerritory;
        finalizeMetrics(rawMetrics);

        final MetricsFile<QuerySortedSeqMetrics, Integer> out = getMetricsFile();
        out.addMetric(usableMetrics.metrics);
        out.addMetric(rawMetrics.metrics);
        out.write(OUTPUT);
        return 0;
    }

    /**
     * Calculate the contribution to the intermediate metrics for a given read pair
     *
     * @param pairToAnalyze the read pair to grab metrics from
     * @param metrics the intermediate metrics with all the data we need
     * @param minimumMappingQuality the minimum mapping quality
     * @param minimumBaseQuality the minimum base quality
     */
    private void calculateMetricsForRead(final QuerySortedReadPairIteratorUtil.ReadPair pairToAnalyze,
                                         final IntermediateMetrics metrics,
                                         final int minimumMappingQuality,
                                         final int minimumBaseQuality) {

        final boolean isPaired = (pairToAnalyze.read2 != null);

        // how many bases do we have?
        final int read1bases = pairToAnalyze.read1.getReadLength();
        final int read2bases = isPaired ? pairToAnalyze.read2.getReadLength() : 0;
        final int totalReadBases = read1bases + read2bases;

        // now compute metrics...
        metrics.metrics.TOTAL_BASES += totalReadBases;
        if (isPaired) metrics.metrics.TOTAL_READ_PAIRS++;

        if (!isPaired || pairToAnalyze.read1.getMateUnmappedFlag() || pairToAnalyze.read2.getMateUnmappedFlag()) {
            metrics.basesExcludedByPairing += totalReadBases;
        } else if (pairToAnalyze.read1.getDuplicateReadFlag()) {
            metrics.metrics.TOTAL_DUPE_PAIRS++;
            metrics.basesExcludedByDupes += totalReadBases;
        } else {

            // determine the bad bases from the reads
            final BaseExclusionHelper read1exclusions = determineBaseExclusions(pairToAnalyze.read1, minimumMappingQuality, minimumBaseQuality);
            final BaseExclusionHelper read2exclusions = determineBaseExclusions(pairToAnalyze.read2, minimumMappingQuality, minimumBaseQuality);
            metrics.basesExcludedByMapq += read1exclusions.basesExcludedByMapq + read2exclusions.basesExcludedByMapq;
            metrics.basesExcludedByBaseq += read1exclusions.lowBQcount + read2exclusions.lowBQcount;

            // keep track of the total usable bases
            int usableBaseCount = totalReadBases;
            usableBaseCount -= (read1exclusions.basesExcludedByMapq + read1exclusions.lowBQcount);
            usableBaseCount -= (read2exclusions.basesExcludedByMapq + read2exclusions.lowBQcount);

            // subtract out bad bases from overlaps between the reads, but only if both reads pass mapping quality thresholds
            if (read1exclusions.basesExcludedByMapq == 0 && read2exclusions.basesExcludedByMapq == 0) {
                final int overlapCount = getOverlappingBaseCount(read1exclusions, read2exclusions, minimumBaseQuality);
                metrics.basesExcludedByOverlap += overlapCount;
                usableBaseCount -= overlapCount;
            }

            metrics.metrics.TOTAL_PASSING_BASES += usableBaseCount;

            final int insertSize = Math.abs(pairToAnalyze.read1.getInferredInsertSize());
            if (insertSize > 0 && pairToAnalyze.read1.getProperPairFlag()) {
                metrics.metrics.TOTAL_ORIENTED_PAIRS++;
                metrics.insertSizeSum += insertSize;
            }
        }
    }

    /**
     * Finalize the metrics by doing some fun but easy math
     *
     * @param metrics the intermediate metrics with all the data we need
     */
    private void finalizeMetrics(final IntermediateMetrics metrics) {
        setUnusedMetrics(metrics.metrics);
        metrics.metrics.MEAN_COVERAGE = metrics.metrics.TOTAL_PASSING_BASES / (double)metrics.metrics.GENOME_TERRITORY;
        metrics.metrics.PCT_EXC_DUPE = metrics.basesExcludedByDupes / (double)metrics.metrics.TOTAL_BASES;
        metrics.metrics.PCT_EXC_MAPQ = metrics.basesExcludedByMapq / (double)metrics.metrics.TOTAL_BASES;
        metrics.metrics.PCT_EXC_UNPAIRED = metrics.basesExcludedByPairing / (double)metrics.metrics.TOTAL_BASES;
        metrics.metrics.PCT_EXC_BASEQ = metrics.basesExcludedByBaseq / (double)metrics.metrics.TOTAL_BASES;
        metrics.metrics.PCT_EXC_OVERLAP = metrics.basesExcludedByOverlap / (double)metrics.metrics.TOTAL_BASES;
        final double totalExcludedBases = metrics.metrics.TOTAL_BASES - metrics.metrics.TOTAL_PASSING_BASES;
        metrics.metrics.PCT_EXC_TOTAL = totalExcludedBases / metrics.metrics.TOTAL_BASES;
        metrics.metrics.MEAN_INSERT_SIZE = metrics.insertSizeSum / metrics.metrics.TOTAL_ORIENTED_PAIRS;
    }

    /**
     * Get the count of low quality and/or softclip bases in the given read
     *
     * @param exclusions  the helper object
     * @param minimumBaseQuality the minimum base quality
     * @return non-negative int
     */
    private int getLowQualityOrSoftclipBaseCount(final BaseExclusionHelper exclusions, final int minimumBaseQuality) {
        final byte[] quals = exclusions.read.getBaseQualities();

        int badCount = exclusions.firstUnclippedBaseIndex + (quals.length - exclusions.firstTrailingClippedBaseIndex);
        for (int i = exclusions.firstUnclippedBaseIndex; i < exclusions.firstTrailingClippedBaseIndex; i++) {
            if (quals[i] < minimumBaseQuality)
                badCount++;
        }
        return badCount;
    }

    /**
     * set the values of the unused metrics to -1
     *
     * @param metrics the metrics object
     */
    private void setUnusedMetrics(final QuerySortedSeqMetrics metrics) {
        metrics.SD_COVERAGE = -1;
        metrics.MEDIAN_COVERAGE = -1;
        metrics.MAD_COVERAGE = -1;
        metrics.PCT_1X = -1;
        metrics.PCT_5X = -1;
        metrics.PCT_10X = -1;
        metrics.PCT_15X = -1;
        metrics.PCT_20X = -1;
        metrics.PCT_25X = -1;
        metrics.PCT_30X = -1;
        metrics.PCT_40X = -1;
        metrics.PCT_50X = -1;
        metrics.PCT_60X = -1;
        metrics.PCT_70X = -1;
        metrics.PCT_80X = -1;
        metrics.PCT_90X = -1;
        metrics.PCT_100X = -1;
        metrics.PCT_EXC_CAPPED = -1;
    }

    /**
     * Get the count of overlapping bases for the given reads
     *
     * @param read1exclusions  the 1st read exclusions
     * @param read2exclusions  the 2nd read exclusions
     * @param minimumBaseQuality the minimum base quality
     * @return non-negative int
     */
    private int getOverlappingBaseCount(final BaseExclusionHelper read1exclusions, final BaseExclusionHelper read2exclusions, final int minimumBaseQuality) {
        // make life easy by ensuring that reads come in coordinate order
        if ( read2exclusions.read.getAlignmentStart() < read1exclusions.read.getAlignmentStart() ) {
            return getOverlappingBaseCount(read2exclusions, read1exclusions, minimumBaseQuality);
        }

        // must be overlapping
        if ( read1exclusions.read.getAlignmentEnd() < read2exclusions.read.getAlignmentStart() ||
                !read1exclusions.read.getReferenceIndex().equals(read2exclusions.read.getReferenceIndex()) )
            return 0;

        final byte[] read1quals = read1exclusions.read.getBaseQualities();
        final byte[] read2quals = read2exclusions.read.getBaseQualities();
        final int indexOfOverlapInFirstRead = read1exclusions.read.getReadPositionAtReferencePosition(read2exclusions.read.getAlignmentStart(), true) - 1;
        final int maxPossibleOverlap = read1exclusions.firstTrailingClippedBaseIndex - indexOfOverlapInFirstRead;
        // the overlap cannot actually be larger than the usable bases in read2
        final int actualOverlap = Math.min(maxPossibleOverlap, read2exclusions.firstTrailingClippedBaseIndex - read2exclusions.firstUnclippedBaseIndex);
        int numHighQualityOverlappingBases = 0;

        for (int i = 0; i < actualOverlap; i++) {
            // we count back from the end of the aligned bases (i.e. not included soft-clips) in read1 and from the front of read2
            final int posInRead1 = read1exclusions.firstTrailingClippedBaseIndex - actualOverlap + i;
            final int posInRead2 = read2exclusions.firstUnclippedBaseIndex + i;

            // we only want to count it if they are both high quality (i.e. not already counted among bad bases)
            if (read1quals[posInRead1] >= minimumBaseQuality && read2quals[posInRead2] >= minimumBaseQuality) {
                numHighQualityOverlappingBases++;
            }
        }

        return numHighQualityOverlappingBases;
    }

    /**
     * Determine how many bases are excluded because of low mapping or base quality.
     *
     * @param read the read
     * @param minimumMappingQuality the minimum mapping quality
     * @param minimumBaseQuality the minimum base quality
     * @return non-null object
     */
    private BaseExclusionHelper determineBaseExclusions(final SAMRecord read, final int minimumMappingQuality, final int minimumBaseQuality) {
        final BaseExclusionHelper exclusions = new BaseExclusionHelper(read);

        if (read.getMappingQuality() < minimumMappingQuality) {
            exclusions.basesExcludedByMapq = read.getReadLength();
        } else {
            exclusions.lowBQcount = getLowQualityOrSoftclipBaseCount(exclusions, minimumBaseQuality);
        }

        return exclusions;
    }

    private static class BaseExclusionHelper {
        public SAMRecord read;
        public int firstUnclippedBaseIndex;
        public int firstTrailingClippedBaseIndex;
        public int basesExcludedByMapq = 0;
        public int lowBQcount = 0;

        public BaseExclusionHelper(final SAMRecord read) {
            this.read = read;

            final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
            firstUnclippedBaseIndex = 0;
            for (final CigarElement element : cigarElements) {
                final CigarOperator op = element.getOperator();

                if (op == CigarOperator.SOFT_CLIP) {
                    firstUnclippedBaseIndex = element.getLength();
                } else if (op != CigarOperator.HARD_CLIP) {
                    break;
                }
            }

            firstTrailingClippedBaseIndex = read.getReadLength();
            for (int i = cigarElements.size() - 1; i >= 0; --i) {
                final CigarElement element = cigarElements.get(i);
                final CigarOperator op = element.getOperator();

                if (op == CigarOperator.SOFT_CLIP) {
                    firstTrailingClippedBaseIndex -= element.getLength();
                } else if (op != CigarOperator.HARD_CLIP) {
                    break;
                }
            }
        }
    }
}

