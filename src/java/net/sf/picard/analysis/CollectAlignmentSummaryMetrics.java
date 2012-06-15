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

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.CollectionUtil;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.*;
import net.sf.picard.util.IlluminaUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.*;

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
public class CollectAlignmentSummaryMetrics extends SinglePassSamProgram {

    private static final Log log = Log.getInstance(CollectAlignmentSummaryMetrics.class);

    // Usage and parameters
    @Usage
    public String USAGE = "Reads a SAM or BAM file and writes a file containing summary alignment metrics.\n";

    @Option(doc="Paired end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.")
    public int MAX_INSERT_SIZE = 100000;

    @Option() public List<String> ADAPTER_SEQUENCE = CollectionUtil.makeList(
        IlluminaUtil.IlluminaAdapterPair.SINGLE_END.get5PrimeAdapter(),
        IlluminaUtil.IlluminaAdapterPair.SINGLE_END.get3PrimeAdapter(),
        IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get5PrimeAdapter(),
        IlluminaUtil.IlluminaAdapterPair.PAIRED_END.get3PrimeAdapter(),
        IlluminaUtil.IlluminaAdapterPair.INDEXED.get5PrimeAdapter(),
        IlluminaUtil.IlluminaAdapterPair.INDEXED.get3PrimeAdapter()
    );

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Option(shortName="BS", doc="Whether the SAM or BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    private AlignmentSummaryMetricsCollector collector;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectAlignmentSummaryMetrics().instanceMainWithExit(argv);
    }

    /** Silly method that is necessary to give unit test access to call doWork() */
    protected final int testDoWork() { return doWork(); }

    @Override protected void setup(final SAMFileHeader header, final File samFile) {
        IoUtil.assertFileIsWritable(OUTPUT);

        if (header.getSequenceDictionary().isEmpty()) {
            log.warn(INPUT.getAbsoluteFile() + " has no sequence dictionary.  If any reads " +
                    "in the file are aligned then alignment summary metrics collection will fail.");
        }

        final boolean doRefMetrics = REFERENCE_SEQUENCE != null;
        collector = new AlignmentSummaryMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), doRefMetrics,
                ADAPTER_SEQUENCE, MAX_INSERT_SIZE,  IS_BISULFITE_SEQUENCED);
    }

    @Override protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        collector.acceptRecord(rec, ref);
    }

    @Override protected void finish() {
        collector.finish();

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> file = getMetricsFile();
        collector.addAllLevelsToFile(file);

        file.write(OUTPUT);
    }
}
