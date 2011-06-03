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

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.MetricsFile;

import java.io.File;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.util.IntervalList;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.SequenceUtil;

/**
 * Calculates a set of HS metrics from a sam or bam file.
 *
 * @author Tim Fennell
 */
public class CalculateHsMetrics extends CommandLineProgram {
    @Usage public final String USAGE =
            "Calculates a set of Hybrid Selection specific metrics from an aligned SAM" +
            "or BAM file. If a reference sequence is provided, AT/GC dropout metrics will " +
            "be calculated, and the PER_TARGET_COVERAGE option can be used to output GC and " +
            "mean coverage information for every target.";
    @Option(shortName="BI", doc="An interval list file that contains the locations of the baits used.")
    public File BAIT_INTERVALS;

    @Option(shortName="TI", doc="An interval list file that contains the locations of the targets.")
    public File TARGET_INTERVALS;

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="An aligned SAM or BAM file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file to write the metrics to.")
    public File OUTPUT;

    @Option(shortName="M", mutex="OUTPUT", doc="Legacy synonym for OUTPUT, should not be used.")
    public File METRICS_FILE;

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

        if (REFERENCE_SEQUENCE == null && PER_TARGET_COVERAGE != null) {
            throw new IllegalArgumentException("Must supply REFERENCE_SEQUENCE when supplying PER_TARGET_COVERAGE");
        }

        final HsMetricsCalculator calculator = new HsMetricsCalculator(BAIT_INTERVALS, TARGET_INTERVALS);

        final SAMFileReader sam = new SAMFileReader(INPUT);

        // Validate that the targets and baits fore for the same references as the reads files
        SequenceUtil.assertSequenceDictionariesEqual(sam.getFileHeader().getSequenceDictionary(),
                IntervalList.fromFile(TARGET_INTERVALS).getHeader().getSequenceDictionary());
        SequenceUtil.assertSequenceDictionariesEqual(sam.getFileHeader().getSequenceDictionary(),
                IntervalList.fromFile(BAIT_INTERVALS).getHeader().getSequenceDictionary());
        
        final ReferenceSequenceFile ref;
        if (REFERENCE_SEQUENCE != null) {
            IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
            ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);
            SequenceUtil.assertSequenceDictionariesEqual(sam.getFileHeader().getSequenceDictionary(), ref.getSequenceDictionary());
            calculator.setReference(ref);
        }
        else {
            ref = null;
        }

        calculator.setPerTargetOutput(PER_TARGET_COVERAGE);
        calculator.analyze(sam.iterator());

        final MetricsFile<HsMetrics, Integer> metrics = getMetricsFile();
        metrics.addMetric(calculator.getMetrics());

        metrics.write(OUTPUT);
        return 0;
    }

}
