/*
 * The MIT License
 *
 * Copyright (c) 2011 The Broad Institute
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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.analysis.directed.RnaSeqMetricsCollector;
import picard.annotation.Gene;
import picard.annotation.GeneAnnotationReader;
import picard.cmdline.programgroups.Metrics;
import picard.util.RExecutor;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@CommandLineProgramProperties(
        summary = CollectRnaSeqMetrics.USAGE_SUMMARY + CollectRnaSeqMetrics.USAGE_DETAILS,
        oneLineSummary = CollectRnaSeqMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectRnaSeqMetrics extends SinglePassSamProgram {
static final String USAGE_SUMMARY = "Produces RNA alignment metrics for a SAM or BAM file.  ";
static final String USAGE_DETAILS = "<p>This tool takes a SAM/BAM file containing the aligned reads from an RNAseq experiment "+
"and produces metrics describing the distribution of the bases within the transcripts.  It calculates the total numbers and the "+
"fractions of nucleotides within specific genomic regions including untranslated regions (UTRs), introns, intergenic sequences "+
"(between discrete genes), and peptide-coding sequences (exons). This tool also determines the numbers of bases that pass quality filters "+
"that are specific to Illumina data (PF_BASES).  For more information please see the corresponding GATK "+
"<a href='https://www.broadinstitute.org/gatk/guide/article?id=6329'>Dictionary</a> entry.</p>" +

"<p>Other metrics include the median coverage (depth), the ratios of 5 prime /3 prime-biases, and the numbers of reads with the "+
"correct/incorrect strand designation. The 5 prime /3 prime-bias results from errors introduced by reverse transcriptase enzymes "+
"during library construction, ultimately leading to the over-representation of either the 5 prime or 3 prime ends of transcripts.  "+
"Please see the CollectRnaSeqMetrics "+
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics'>definitions</a> "+
"for details on how these biases are calculated. </p>" +

"<p>The sequence input must be a valid SAM/BAM file containing RNAseq data aligned by an RNAseq-aware genome aligner such a "+
"<a href='http://github.com/alexdobin/STAR'>STAR</a> or <a href='http://ccb.jhu.edu/software/tophat/index.shtml'>TopHat</a>. "+
"The tool also requires a REF_FLAT file, a tab-delimited file containing information about the location of RNA transcripts, "+
"exon start and stop sites, etc. For an example refFlat file for GRCh38, see refFlat.txt.gz at "+
"<a href='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database'>http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database</a>.  "+
"The first five lines of the tab-limited text file appear as follows.</p>"+

"<pre>" +
"DDX11L1	NR_046018	chr1	+	11873	14409	14409	14409	3	11873,12612,13220,	12227,12721,14409," +
"WASH7P	NR_024540	chr1	-	14361	29370	29370	29370	11	14361,14969,15795,16606,16857,17232,17605,17914,18267,24737,29320,	14829,15038,15947,16765,17055,17368,17742,18061,18366,24891,29370," +
"DLGAP2-AS1	NR_103863	chr8_KI270926v1_alt	-	33083	35050	35050	35050	3	33083,33761,35028,	33281,33899,35050," +
"MIR570	NR_030296	chr3	+	195699400	195699497	195699497	195699497	1	195699400,	195699497," +
"MIR548A3	NR_030330	chr8	-	104484368	104484465	104484465	104484465	1	104484368,	104484465," +
"</pre>" +

"<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>"+
"<h4>Usage example:</h4>"+
"<pre>" +
"java -jar picard.jar CollectRnaSeqMetrics \\<br />" +
"      I=input.bam \\<br />" +
"      O=output.RNA_Metrics \\<br />" +
"      REF_FLAT=ref_flat.txt \\<br />" +
"      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \\<br />" +
"      RIBOSOMAL_INTERVALS=ribosomal.interval_list" +
"</pre>" +
"Please see the CollectRnaSeqMetrics " +
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics'>definitions</a> " +
"for a complete description of the metrics produced by this tool." +
"<hr />"
;

    private static final Log LOG = Log.getInstance(CollectRnaSeqMetrics.class);

    @Argument(doc="Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat")
    public File REF_FLAT;

    @Argument(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.  " +
            "Format described <a href=\"http://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/IntervalList.html\">here</a>:", optional = true)
    public File RIBOSOMAL_INTERVALS;

    @Argument(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public RnaSeqMetricsCollector.StrandSpecificity STRAND_SPECIFICITY;

    @Argument(doc="When calculating coverage based values (e.g. CV of coverage) only use transcripts of this length or greater.")
    public int MINIMUM_LENGTH = 500;

    @Argument(doc="The PDF file to write out a plot of normalized position vs. coverage.", shortName="CHART", optional = true)
    public File CHART_OUTPUT;

    @Argument(doc="If a read maps to a sequence specified with this option, all the bases in the read are counted as ignored bases.  " +
    "These reads are not counted as ")
    public Set<String> IGNORE_SEQUENCE = new HashSet<String>();

    @Argument(doc="This percentage of the length of a fragment must overlap one of the ribosomal intervals for a read or read pair to be considered rRNA.")
    public double RRNA_FRAGMENT_PERCENTAGE = 0.8;

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    private RnaSeqMetricsCollector collector;

    /**
     * A subtitle for the plot, usually corresponding to a library.
     */
    private String plotSubtitle = "";

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectRnaSeqMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected String[] customCommandLineValidation() {
        // No ribosomal intervals file and rRNA fragment percentage = 0
        if ( RIBOSOMAL_INTERVALS == null && RRNA_FRAGMENT_PERCENTAGE == 0 ) {
            throw new PicardException("Must use a RIBOSOMAL_INTERVALS file if RRNA_FRAGMENT_PERCENTAGE = 0.0");
        }
        return super.customCommandLineValidation();
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {

        if (CHART_OUTPUT != null) IOUtil.assertFileIsWritable(CHART_OUTPUT);

        final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadRefFlat(REF_FLAT, header.getSequenceDictionary());
        LOG.info("Loaded " + geneOverlapDetector.getAll().size() + " genes.");

        final Long ribosomalBasesInitialValue = RIBOSOMAL_INTERVALS != null ? 0L : null;
        final OverlapDetector<Interval> ribosomalSequenceOverlapDetector = RnaSeqMetricsCollector.makeOverlapDetector(samFile, header, RIBOSOMAL_INTERVALS, LOG);

        final HashSet<Integer> ignoredSequenceIndices = RnaSeqMetricsCollector.makeIgnoredSequenceIndicesSet(header, IGNORE_SEQUENCE);

        collector = new RnaSeqMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), ribosomalBasesInitialValue,
                geneOverlapDetector, ribosomalSequenceOverlapDetector, ignoredSequenceIndices, MINIMUM_LENGTH, STRAND_SPECIFICITY, RRNA_FRAGMENT_PERCENTAGE,
                true);

        // If we're working with a single library, assign that library's name as a suffix to the plot title
        final List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        if (readGroups.size() == 1) {
            this.plotSubtitle = readGroups.get(0).getLibrary();
            if (null == this.plotSubtitle) this.plotSubtitle = "";
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence refSeq) {
        collector.acceptRecord(rec, refSeq);
    }

    @Override
    protected void finish() {
        collector.finish();

        final MetricsFile<RnaSeqMetrics, Integer> file = getMetricsFile();
        collector.addAllLevelsToFile(file);
        file.write(OUTPUT);

        boolean atLeastOneHistogram = false;
        for (final Histogram<Integer> histo : file.getAllHistograms()) {
            atLeastOneHistogram = atLeastOneHistogram || !histo.isEmpty();
        }
        // Generate the coverage by position plot
        if (CHART_OUTPUT != null && atLeastOneHistogram) {
            final int rResult = RExecutor.executeFromClasspath("picard/analysis/rnaSeqCoverage.R",
                                                               OUTPUT.getAbsolutePath(),
                                                               CHART_OUTPUT.getAbsolutePath(),
                                                               INPUT.getName(),
                                                               this.plotSubtitle);

            if (rResult != 0) {
                throw new PicardException("Problem invoking R to generate plot.");
            }
        }
    }

}
