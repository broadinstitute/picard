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

package net.sf.picard.analysis.directed;

import net.sf.picard.analysis.MetricAccumulationLevel;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;

import java.io.File;
import java.util.*;

import net.sf.picard.util.*;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SequenceUtil;

/**
 * Calculates a set of HS metrics from a sam or bam file.
 *
 * @author Tim Fennell
 */
public class CalculateHsMetrics extends CommandLineProgram {

    private static final Log log = Log.getInstance(CalculateHsMetrics.class);
    private HsMetricsCalculator baseCalculator = null;

    @Usage public final String USAGE =
            "Calculates a set of Hybrid Selection specific metrics from an aligned SAM" +
            "or BAM file. If a reference sequence is provided, AT/GC dropout metrics will " +
            "be calculated, and the PER_TARGET_COVERAGE option can be used to output GC and " +
            "mean coverage information for every target.";
    @Option(shortName="BI", doc="An interval list file that contains the locations of the baits used.")
    public File BAIT_INTERVALS;

    @Option(shortName="TI", doc="An interval list file that contains the locations of the targets.")
    public File TARGET_INTERVALS;
    
    @Option(shortName="N", doc="Bait set name. If not provided it is inferred from the filename of the bait intervals.")
    public String BAIT_SET_NAME;

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="An aligned SAM or BAM file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file to write the metrics to.")
    public File OUTPUT;

    @Option(shortName="M", mutex="OUTPUT", doc="Legacy synonym for OUTPUT, should not be used.")
    public File METRICS_FILE;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, optional=true, doc="The reference sequence aligned to.")
    public File REFERENCE_SEQUENCE;

    @Option(optional=true, doc="An optional file to output per target coverage information to.")
    public File PER_TARGET_COVERAGE;

    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new CalculateHsMetrics().instanceMain(argv));
    }

    /**
     * Asserts that files are readable and writable and then fires off an
     * HsMetricsCalculator instance to do the real work.
     */
    protected int doWork() {
        if (OUTPUT==null) OUTPUT = METRICS_FILE;

        IoUtil.assertFileIsReadable(BAIT_INTERVALS);
        IoUtil.assertFileIsReadable(TARGET_INTERVALS);
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        if (PER_TARGET_COVERAGE != null) IoUtil.assertFileIsWritable(PER_TARGET_COVERAGE);

        final SAMFileReader samReader = new SAMFileReader(INPUT);

        if (REFERENCE_SEQUENCE == null && PER_TARGET_COVERAGE != null) {
            throw new IllegalArgumentException("Must supply REFERENCE_SEQUENCE when supplying PER_TARGET_COVERAGE");
        }

        final boolean calculateAll = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.ALL_READS);
        final boolean calculateSample = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.SAMPLE);
        final boolean calculateLibrary = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.LIBRARY);
        final boolean calculateReadGroup = METRIC_ACCUMULATION_LEVEL.contains(MetricAccumulationLevel.READ_GROUP);

        // Validate that the targets and baits have the same references as the reads file
        SequenceUtil.assertSequenceDictionariesEqual(samReader.getFileHeader().getSequenceDictionary(),
                IntervalList.fromFile(TARGET_INTERVALS).getHeader().getSequenceDictionary(),
                INPUT, TARGET_INTERVALS);
        SequenceUtil.assertSequenceDictionariesEqual(samReader.getFileHeader().getSequenceDictionary(),
                IntervalList.fromFile(BAIT_INTERVALS).getHeader().getSequenceDictionary(),
                INPUT, BAIT_INTERVALS);

        ReferenceSequenceFile ref = null;
        if (REFERENCE_SEQUENCE != null) {
            IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
            SequenceUtil.assertSequenceDictionariesEqual(samReader.getFileHeader().getSequenceDictionary(), ref.getSequenceDictionary(),
                INPUT, REFERENCE_SEQUENCE);
        }

        final HsMetricsCalculator allReadsCalculator = calculateAll
              ? createCalculator(null, null, null, ref)
              : null;
        if (calculateAll && PER_TARGET_COVERAGE != null) {
            allReadsCalculator.setPerTargetOutput(PER_TARGET_COVERAGE);
        }

        final Map<String,HsMetricsCalculator> sampleCalculators    = new HashMap<String,HsMetricsCalculator>();
        final Map<String,HsMetricsCalculator> libraryCalculators   = new HashMap<String,HsMetricsCalculator>();
        final Map<String,HsMetricsCalculator> readGroupCalculators = new HashMap<String,HsMetricsCalculator>();

        for (final SAMReadGroupRecord rg : samReader.getFileHeader().getReadGroups()) {
            if (calculateSample) {
                if (!sampleCalculators.containsKey(rg.getSample())) {
                    sampleCalculators.put(rg.getSample(), createCalculator(rg.getSample(), null, null, ref));
                }
            }
            if (calculateLibrary) {
                if (!libraryCalculators.containsKey(rg.getLibrary())) {
                    libraryCalculators.put(rg.getLibrary(), createCalculator(rg.getSample(), rg.getLibrary(), null, ref));
                }
            }
            if (calculateReadGroup) {
                if (!readGroupCalculators.containsKey(rg.getPlatformUnit())) {
                    readGroupCalculators.put(rg.getPlatformUnit(), createCalculator(rg.getSample(), rg.getLibrary(),
                            rg.getPlatformUnit(), ref));
                }
            }
        }


        // Add each record to the requested collectors
        final Iterator<SAMRecord> records = samReader.iterator();
        int i = 0;
        while (records.hasNext()) {
            final SAMRecord sam = records.next();

            if (calculateAll) {
                allReadsCalculator.analyze(sam);
            }

            if (calculateSample) {
                sampleCalculators.get(sam.getReadGroup().getSample()).analyze(sam);
            }

            if (calculateLibrary) {
                libraryCalculators.get(sam.getReadGroup().getLibrary()).analyze(sam);
            }

            if (calculateReadGroup) {
                readGroupCalculators.get(sam.getReadGroup().getPlatformUnit()).analyze(sam);
            }

            if (++i % 1000000 == 0) {
                log.info("Processed " + i + " records so far.");
            }
        }

        // Write the output file
        final MetricsFile<HsMetrics, Integer> metrics = getMetricsFile();

        if (calculateAll)  metrics.addMetric(allReadsCalculator.getMetrics());

        if (calculateSample) {
            for (final HsMetricsCalculator calculator : sampleCalculators.values()) {
                metrics.addMetric(calculator.getMetrics());
            }
        }

        if (calculateLibrary) {
            for (final HsMetricsCalculator calculator : libraryCalculators.values()) {
                metrics.addMetric(calculator.getMetrics());
            }
        }

        if (calculateReadGroup) {
            for (final HsMetricsCalculator calculator : readGroupCalculators.values()) {
                metrics.addMetric(calculator.getMetrics());
            }
        }

        metrics.write(OUTPUT);
        return 0;
    }

    private HsMetricsCalculator createCalculator(final String sample, final String library, final String readGroup,
                                                 final ReferenceSequenceFile ref) {
        if (baseCalculator == null) {
            baseCalculator = new HsMetricsCalculator(BAIT_INTERVALS, TARGET_INTERVALS, ref, sample, library, readGroup);
            if (this.BAIT_SET_NAME != null) baseCalculator.setBaitSetName(BAIT_SET_NAME);
            return baseCalculator;
        }
        else {
            return baseCalculator.cloneMetricsCalculator(sample, library, readGroup);
        }
    }

    protected String[] customCommandLineValidation() {
        if (PER_TARGET_COVERAGE != null && (METRIC_ACCUMULATION_LEVEL.size() != 1 ||
            METRIC_ACCUMULATION_LEVEL.iterator().next() != MetricAccumulationLevel.ALL_READS)) {
            return new String[] {"PER_TARGET_COVERAGE can be specified only when METRIC_ACCUMULATION_LEVEL is set " +
                    "to ALL_READS."};
        }
        return super.customCommandLineValidation();
    }

}
