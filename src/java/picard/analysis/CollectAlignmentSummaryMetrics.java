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

package picard.analysis;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.util.EnumSet;
import java.util.List;
import java.util.Set;

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
 * Chimeras are identified if any of the of following criteria are met:
 * <ul>
 * <li>the insert size is larger than MAX_INSERT_SIZE</li>
 * <li>the ends of a pair map to different contigs</li>
 * <li>the paired end orientation is different that the expected orientation</li>
 * <li>the read contains an SA tag (chimeric alignment)</li>
 * </ul>
 * 
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
@CommandLineProgramProperties(
        usage = CollectAlignmentSummaryMetrics.USAGE_SUMMARY + CollectAlignmentSummaryMetrics.USAGE_DETAILS,
        usageShort = CollectAlignmentSummaryMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectAlignmentSummaryMetrics extends SinglePassSamProgram {
    static final String USAGE_SUMMARY = "Produce a summary of alignment metrics from a SAM or BAM file";
    static final String USAGE_DETAILS = "Using read outputs from high throughput sequencing (HTS) technologies, this tool provides " +
            "metrics regarding the quality of read alignments to a reference sequence, as well as the proportion of the reads " +
            "that passed machine signal-to-noise threshold quality filters (Illumina)."+
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "    java -jar picard.jar CollectAlignmentSummaryMetrics \\<br />" +
            "          R=reference_sequence.fasta \\<br />" +
            "          I=input.bam \\<br />" +
            "          O=output.txt" +
            "</pre>"+
            "Please see <a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics'>" +
            "the AlignmentSummaryMetrics documentation</a> for detailed explanations of each metric. <br /> <br />" +
            "Additional information about Illumina's quality filters can be found in the following documents on the Illumina website: " +
            "<ul><li>http://support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-x-percent-pf-technical-note-770-2014-043.pdf</li> " +
            "<li>http://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/hiseqx/hiseq-x-system-guide-15050091-d.pdf</li></ul>" +
            "<hr />";
    private static final Log log = Log.getInstance(CollectAlignmentSummaryMetrics.class);

    @Option(doc="Paired-end reads above this insert size will be considered chimeric along with inter-chromosomal pairs.")
    public int MAX_INSERT_SIZE = ChimeraUtil.DEFAULT_INSERT_SIZE_LIMIT;

    @Option(doc="Paired-end reads that do not have this expected orientation will be considered chimeric.")
    public Set<PairOrientation> EXPECTED_PAIR_ORIENTATIONS = EnumSet.copyOf(ChimeraUtil.DEFAULT_EXPECTED_ORIENTATIONS);

    @Option(doc="List of adapter sequences to use when processing the alignment metrics")
	public List<String> ADAPTER_SEQUENCE = AdapterUtility.DEFAULT_ADAPTER_SEQUENCE;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Option(shortName="BS", doc="Whether the SAM or BAM file consists of bisulfite sequenced reads.  ")
    public boolean IS_BISULFITE_SEQUENCED = false;

    //overridden to make it visible on the commandline and to change the doc.
    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc = "Reference sequence file. Note that while this argument isn't required, without it only a small subset of the metrics will be calculated. Note also that if a reference sequence is provided, it must be accompanied by a sequence dictionary.",  optional = true, overridable = true)
    public File REFERENCE_SEQUENCE = Defaults.REFERENCE_FASTA;

    private AlignmentSummaryMetricsCollector collector;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectAlignmentSummaryMetrics().instanceMainWithExit(argv);
    }

    /** Silly method that is necessary to give unit test access to call doWork() */
    protected final int testDoWork() { return doWork(); }

    @Override protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(OUTPUT);

        if (header.getSequenceDictionary().isEmpty()) {
            log.warn(INPUT.getAbsoluteFile() + " has no sequence dictionary.  If any reads " +
                    "in the file are aligned, then alignment summary metrics collection will fail.");
        }

        final boolean doRefMetrics = REFERENCE_SEQUENCE != null;
        collector = new AlignmentSummaryMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), doRefMetrics,
                ADAPTER_SEQUENCE, MAX_INSERT_SIZE, EXPECTED_PAIR_ORIENTATIONS, IS_BISULFITE_SEQUENCED);
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
