/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package net.sf.picard.analysis;

import net.sf.picard.util.CollectionUtil;
import net.sf.picard.util.Histogram;
import net.sf.picard.sam.ReservedTagConstants;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.*;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.analysis.AlignmentSummaryMetrics.Category;
import net.sf.picard.util.IlluminaUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CoordMath;
import net.sf.samtools.util.SequenceUtil;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;

import java.io.File;
import java.util.*;

/**
 * A command line tool to read a BAM file and produce standard alignment metrics that would be applicable to any alignment.  
 * Metrics to include, but not limited to:
 * <ul>
 * <li>Total number of reads (total, period, no exclusions)</li>
 * <li>Total number of PF reads (PF == does not fail vendor check flag)</li>
 * <li>Number of PF noise reads (does not fail vendor check and has noise attr set)</li>
 * <li>Total aligned PF reads (any PF read that has a sequence and position)</li>
 * <li>High quality aligned PF reads (high quality == mapping quality >= 20)</li>
 * <li>High quality aligned PF bases (actual aligned bases, calculate off alignment blocks)</li>
 * <li>High quality aligned PF Q20 bases (subset of above where base quality >= 20)</li>
 * <li>Median mismatches in HQ aligned PF reads (how many aligned bases != ref on average)</li>
 * <li>Reads aligned in pairs (vs. reads aligned with mate unaligned/not present)</li>
 * <li>Read length (how to handle mixed lengths?)</li>
 * <li>Bad Cycles - how many machine cycles yielded combined no-call and mismatch rates of >= 80%</li>
 * <li>Strand balance - reads mapped to positive strand / total mapped reads</li>
 * </ul>
 * Metrics are written for the first read of a pair, the second read, and combined for the pair.
 * 
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class CollectAlignmentSummaryMetrics extends CommandLineProgram {
    private static final int MAPPING_QUALITY_THRESHOLD = 20;
    private static final int BASE_QUALITY_THRESHOLD = 20;

    private static final int ADAPTER_MATCH_LENGTH = 16;
    private static final int MAX_ADAPTER_ERRORS = 1;
    private byte[][] ADAPTER_SEQUENCES;
    private static final Log log = Log.getInstance(CollectAlignmentSummaryMetrics.class);



    // Usage and parameters
    @Usage
    public String USAGE = "Reads a SAM or BAM file and writes a file containing summary alignment metrics.\n";
    @Option(shortName="I", doc="SAM or BAM file") public File INPUT;
    @Option(shortName="O", doc="File to write alignment summary metrics to") public File OUTPUT;
    @Option(shortName="R", doc="Reference sequence file", optional=true) public File REFERENCE_SEQUENCE;
    @Option(doc="If true (default), \"unsorted\" SAM/BAM files will be considerd coordinate sorted",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public Boolean ASSUME_SORTED = Boolean.TRUE;
    @Option(doc="Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.")
    public int MAX_INSERT_SIZE = 100000;
    @Option() public List<String> ADAPTER_SEQUENCE = CollectionUtil.makeList(
        IlluminaUtil.AdapterPair.SINGLE_END.get5PrimeAdapter(),
        IlluminaUtil.AdapterPair.SINGLE_END.get3PrimeAdapter(),
        IlluminaUtil.AdapterPair.PAIRED_END.get5PrimeAdapter(),
        IlluminaUtil.AdapterPair.PAIRED_END.get3PrimeAdapter(),
        IlluminaUtil.AdapterPair.INDEXED.get5PrimeAdapter(),
        IlluminaUtil.AdapterPair.INDEXED.get3PrimeAdapter()
    );
    // also List options are only appended to, not replaced
    @Option(shortName="BS", doc="Whether the SAM or BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    private ReferenceSequenceFileWalker referenceSequenceWalker;
    private SAMFileHeader samFileHeader;
    private boolean doRefMetrics;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        System.exit(new CollectAlignmentSummaryMetrics().instanceMain(argv));
    }

    @Override
    protected int doWork() {
        prepareAdapterSequences();

        doRefMetrics = REFERENCE_SEQUENCE != null;

        // Check the files are readable/writable
        IoUtil.assertFileIsReadable(INPUT);
        if(doRefMetrics) IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IoUtil.assertFileIsWritable(OUTPUT);
        final SAMFileReader in = new SAMFileReader(INPUT);
        assertCoordinateSortOrder(in);

        if(doRefMetrics) {this.referenceSequenceWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);}
        this.samFileHeader = in.getFileHeader();

        if (!samFileHeader.getSequenceDictionary().isEmpty()) {
            if(doRefMetrics){
                SequenceUtil.assertSequenceDictionariesEqual(

                    samFileHeader.getSequenceDictionary(),
                    this.referenceSequenceWalker.getSequenceDictionary());
            }
        } else {
            log.warn(INPUT.getAbsoluteFile() + " has no sequence dictionary.  If any reads " +
                    "in the file are aligned then alignment summary metrics collection will fail.");
        }

        final MetricCollector<AlignmentSummaryMetrics, SAMRecord> unpairedCollector =  constructCollector(Category.UNPAIRED);
        final MetricCollector<AlignmentSummaryMetrics, SAMRecord> firstOfPairCollector = constructCollector(Category.FIRST_OF_PAIR);
        final MetricCollector<AlignmentSummaryMetrics, SAMRecord> secondOfPairCollector = constructCollector(Category.SECOND_OF_PAIR);
        final MetricCollector<AlignmentSummaryMetrics, SAMRecord> pairCollector =  constructCollector(Category.PAIR);

        // Loop over the reads applying them to the correct collectors
        for (final SAMRecord record : in) {
            if (record.getReadPairedFlag()) {
                if (record.getFirstOfPairFlag()) firstOfPairCollector.addRecord(record);
                else secondOfPairCollector.addRecord(record);

                pairCollector.addRecord(record);
            }
            else {
                unpairedCollector.addRecord(record);
            }

        }

        in.close();

        // Let the collectors do any summary computations etc.
        firstOfPairCollector.onComplete();
        secondOfPairCollector.onComplete();
        pairCollector.onComplete();
        unpairedCollector.onComplete();

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> file = getMetricsFile();

        if (firstOfPairCollector.getMetrics().TOTAL_READS > 0) {
            // override how bad cycle is determined for paired reads, it should be
            // the sum of first and second reads
            pairCollector.getMetrics().BAD_CYCLES = firstOfPairCollector.getMetrics().BAD_CYCLES +
                                                    secondOfPairCollector.getMetrics().BAD_CYCLES;

            file.addMetric(firstOfPairCollector.getMetrics());
            file.addMetric(secondOfPairCollector.getMetrics());
            file.addMetric(pairCollector.getMetrics());
        }

        //if there are no reads in any category then we will returned an unpaired alignment summary metric with all zero values
        if (unpairedCollector.getMetrics().TOTAL_READS > 0 ||
            (firstOfPairCollector.getMetrics().TOTAL_READS == 0 && unpairedCollector.getMetrics().TOTAL_READS == 0)) {
            file.addMetric(unpairedCollector.getMetrics());
        }

        file.write(OUTPUT);

        return 0;
    }

    /**
     * Checks that the SAM is either coordinate sorted according to it's header, or that the header
     * doesn't specify a sort and that the user has supplied the argument to assume that it is
     * sorted correctly.
     */
    private void assertCoordinateSortOrder(final SAMFileReader in) {
        switch (in.getFileHeader().getSortOrder()) {
        case coordinate:
            break;
        case unsorted:
            if (this.ASSUME_SORTED) {
                break;
            }
        default:
            log.warn("May not be able collect summary statistics in file " + INPUT.getAbsoluteFile() +
            " because it is not sorted in coordinate order.  If any of the reads are aligned this will blow up.");
        }
    }

    /** Converts the supplied adapter sequences to byte arrays in both fwd and rc. */
    protected void prepareAdapterSequences() {
        final int count = ADAPTER_SEQUENCE.size();
        ADAPTER_SEQUENCES = new byte[count * 2][];

        for (int i=0; i<count; ++i) {
            final String adapter = ADAPTER_SEQUENCE.get(i).toUpperCase();
            ADAPTER_SEQUENCES[i] = StringUtil.stringToBytes(adapter);
            ADAPTER_SEQUENCES[i + count] = StringUtil.stringToBytes(SequenceUtil.reverseComplement(adapter));

        }
    }

    /**
     * Checks the first ADAPTER_MATCH_LENGTH bases of the read against known adapter sequences and returns
     * true if the read matches an adapter sequence with MAX_ADAPTER_ERRORS mismsatches or fewer.
     *
     * @param read the basecalls for the read in the order and orientation the machine read them
     * @return true if the read matches an adapter and false otherwise 
     */
    protected boolean isAdapterSequence(final byte[] read) {
        StringUtil.toUpperCase(read);
        
        for (final byte[] adapter : ADAPTER_SEQUENCES) {
            final int lastKmerStart = adapter.length - ADAPTER_MATCH_LENGTH;

            for (int adapterStart=0; adapterStart<lastKmerStart; ++adapterStart) {
                int errors = 0;

                for (int i=0; i<ADAPTER_MATCH_LENGTH && errors <= MAX_ADAPTER_ERRORS; ++i) {
                    if (read[i] != adapter[i + adapterStart]) ++errors;
                }

                if (errors <= MAX_ADAPTER_ERRORS) return true;
            }
        }

        return false;
    }

    /** Constructs a metrics collector for the supplied category. */
    private MetricCollector<AlignmentSummaryMetrics, SAMRecord> constructCollector(final Category category) {
        final MetricCollector<AlignmentSummaryMetrics, SAMRecord> collector =
            new AggregateMetricCollector<AlignmentSummaryMetrics, SAMRecord>(new ReadCounter(), new QualityMappingCounter());
        collector.setMetrics(new AlignmentSummaryMetrics());
        collector.getMetrics().CATEGORY = category;
        return collector;
    }

    /**
     * Class that counts reads that match various conditions
     */
    private class ReadCounter implements MetricCollector<AlignmentSummaryMetrics, SAMRecord> {
        private long numPositiveStrand = 0;
        private final Histogram<Integer> readLengthHistogram = new Histogram<Integer>();
        private AlignmentSummaryMetrics metrics;
        private long chimeras;
        private long adapterReads;

        public void addRecord(final SAMRecord record) {
            if (record.getNotPrimaryAlignmentFlag()) {
                // only want 1 count per read so skip non primary alignments
                return;
            }

            metrics.TOTAL_READS++;
            readLengthHistogram.increment(record.getReadBases().length);

            if (!record.getReadFailsVendorQualityCheckFlag()) {
                metrics.PF_READS++;

                if (isNoiseRead(record)) {
                    metrics.PF_NOISE_READS++;
                }
                else if (record.getReadUnmappedFlag()) {
                    // If the read is unmapped see if it's adapter sequence
                    if (isAdapterSequence(record.getReadBases())) {
                        this.adapterReads++;
                    }
                }
                else {
                    if(doRefMetrics) {
                        metrics.PF_READS_ALIGNED++;

                        if (!record.getReadNegativeStrandFlag()) {
                            numPositiveStrand++;
                        }

                        if (record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {
                            metrics.READS_ALIGNED_IN_PAIRS++;

                            // With both reads mapped we can see if this pair is chimeric
                            if (Math.abs(record.getInferredInsertSize()) > MAX_INSERT_SIZE ||
                                 !record.getReferenceIndex().equals(record.getMateReferenceIndex())) {

                                // Check that both ends have mapq > minimum
                                final Integer mateMq = record.getIntegerAttribute("MQ");
                                if (mateMq == null || mateMq >= MAPPING_QUALITY_THRESHOLD &&
                                        record.getMappingQuality() >= MAPPING_QUALITY_THRESHOLD) {
                                    ++this.chimeras;
                                }
                            }
                        }
                    }
                }
            }
        }

        public void onComplete() {
            if (metrics.TOTAL_READS > 0)
            {
                metrics.PCT_PF_READS = (double) metrics.PF_READS / (double) metrics.TOTAL_READS;
                metrics.PCT_ADAPTER = this.adapterReads / (double) metrics.PF_READS;
                metrics.MEAN_READ_LENGTH = readLengthHistogram.getMean();

                if(doRefMetrics) {
                    metrics.PCT_PF_READS_ALIGNED = (double) metrics.PF_READS_ALIGNED / (double) metrics.PF_READS;
                    metrics.PCT_READS_ALIGNED_IN_PAIRS = (double) metrics.READS_ALIGNED_IN_PAIRS/ (double) metrics.PF_READS_ALIGNED;
                    metrics.STRAND_BALANCE = numPositiveStrand / (double) metrics.PF_READS_ALIGNED;
                    metrics.PCT_CHIMERAS = this.chimeras / (double) metrics.PF_HQ_ALIGNED_READS;
                }
            }
        }

        private boolean isNoiseRead(final SAMRecord record) {
            final Object noiseAttribute = record.getAttribute(ReservedTagConstants.XN);
            return (noiseAttribute != null && noiseAttribute.equals(1));
        }

        public void setMetrics(final AlignmentSummaryMetrics metrics) {
            this.metrics = metrics;
        }

        public AlignmentSummaryMetrics getMetrics() {
            return this.metrics;
        }
    }

    /**
     * Class that counts quality mappings & base calls that match various conditions
     */
    private class QualityMappingCounter implements MetricCollector<AlignmentSummaryMetrics, SAMRecord> {
        private final Histogram<Long> mismatchHistogram = new Histogram<Long>();
        private final Histogram<Integer> badCycleHistogram = new Histogram<Integer>();
        private AlignmentSummaryMetrics metrics;

        public void addRecord(final SAMRecord record) {
            if (record.getNotPrimaryAlignmentFlag()) {
                return;
            }
            if (record.getReadUnmappedFlag() || !doRefMetrics) {
                final byte[] readBases = record.getReadBases();
                for (int i = 0; i < readBases.length; i++) {
                    if (SequenceUtil.isNoCall(readBases[i])) {
                        badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                    }
                }
            } else {
                final boolean highQualityMapping = isHighQualityMapping(record);
                if (highQualityMapping) metrics.PF_HQ_ALIGNED_READS++;

                final byte[] readBases = record.getReadBases();
                final byte[] refBases = referenceSequenceWalker.get(record.getReferenceIndex()).getBases();
                final byte[] qualities  = record.getBaseQualities();
                final int refLength = refBases.length;
                long mismatchCount = 0;

                for (final AlignmentBlock alignmentBlock : record.getAlignmentBlocks()) {
                    final int readIndex = alignmentBlock.getReadStart() - 1;
                    final int refIndex  = alignmentBlock.getReferenceStart() - 1;
                    final int length    = alignmentBlock.getLength();

                    for (int i=0; i<length && refIndex+i<refLength; ++i) {
                        final int readBaseIndex = readIndex + i;
                        boolean mismatch = !SequenceUtil.basesEqual(readBases[readBaseIndex], refBases[refIndex+i]);
                        boolean bisulfiteBase = false;
                        if (mismatch && IS_BISULFITE_SEQUENCED) {
                            if ( (record.getReadNegativeStrandFlag() &&
                                   (refBases[refIndex+i] == 'G' || refBases[refIndex+i] =='g') &&
                                   (readBases[readBaseIndex] == 'A' || readBases[readBaseIndex] == 'a'))
                                || ((!record.getReadNegativeStrandFlag()) &&
                                    (refBases[refIndex+i] == 'C' || refBases[refIndex+i] == 'c') &&
                                    (readBases[readBaseIndex] == 'T') || readBases[readBaseIndex] == 't') ) {

                                bisulfiteBase = true;
                                mismatch = false;
                            }
                        }
                        if (highQualityMapping) {
                            metrics.PF_HQ_ALIGNED_BASES++;
                            if (!bisulfiteBase) {
                                metrics.incrementErrorRateDenominator();
                            }
                            if (qualities[readBaseIndex] >= BASE_QUALITY_THRESHOLD) {
                                metrics.PF_HQ_ALIGNED_Q20_BASES++;
                            }
                            if (mismatch) {
                                mismatchCount++;
                            }
                        }
                        if (mismatch || SequenceUtil.isNoCall(readBases[readBaseIndex])) {
                            badCycleHistogram.increment(CoordMath.getCycle(record.getReadNegativeStrandFlag(), readBases.length, i));
                        }
                    }
                }
                mismatchHistogram.increment(mismatchCount);
            }
        }

        private boolean isHighQualityMapping(final SAMRecord record) {
            return !record.getReadFailsVendorQualityCheckFlag() &&
            record.getMappingQuality() >= MAPPING_QUALITY_THRESHOLD;
        }

        public void onComplete() {
            if(metrics.TOTAL_READS > 0) {
                if(doRefMetrics) {
                metrics.PF_HQ_MEDIAN_MISMATCHES = mismatchHistogram.getMedian();
                metrics.PF_HQ_ERROR_RATE = mismatchHistogram.getSum() / (double)metrics.getErrorRateDenominator();
                }

                metrics.BAD_CYCLES = 0;

                for (final Histogram<Integer>.Bin cycleBin : badCycleHistogram.values()) {
                    final double badCyclePercentage = cycleBin.getValue() / metrics.TOTAL_READS;
                    if (badCyclePercentage >= .8) {
                        metrics.BAD_CYCLES++;
                    }
                }
            }
        }

        public void setMetrics(final AlignmentSummaryMetrics metrics) {
            this.metrics = metrics;
        }

        public AlignmentSummaryMetrics getMetrics() {
            return this.metrics;
        }
    }
}
