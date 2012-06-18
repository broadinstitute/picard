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
package net.sf.picard.analysis;

import net.sf.picard.PicardException;
import net.sf.picard.analysis.directed.RnaSeqMetricsCollector;
import net.sf.picard.annotation.Gene;
import net.sf.picard.annotation.GeneAnnotationReader;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.metrics.*;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.*;
import net.sf.samtools.*;

import java.io.File;
import java.util.*;

public class CollectRnaSeqMetrics extends SinglePassSamProgram {
    private static final Log LOG = Log.getInstance(CollectRnaSeqMetrics.class);

    @Usage
    public final String USAGE = getStandardUsagePreamble() +
            "Program to collect metrics about the alignment of RNA to various functional classes of loci in the genome:" +
            " coding, intronic, UTR, intergenic, ribosomal.\n" +
            "Also determines strand-specificity for strand-specific libraries.";


    @Option(doc="Gene annotations in refFlat form.  Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat")
    public File REF_FLAT;

    @Option(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.  " +
            "Format described here: http://picard.sourceforge.net/javadoc/net/sf/picard/util/IntervalList.html", optional = true)
    public File RIBOSOMAL_INTERVALS;

    @Option(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public RnaSeqMetricsCollector.StrandSpecificity STRAND_SPECIFICITY;

    @Option(doc="When calculating coverage based values (e.g. CV of coverage) only use transcripts of this length or greater.")
    public int MINIMUM_LENGTH = 500;

    @Option(doc="The PDF file to write out a plot of normalized position vs. coverage.", shortName="CHART", optional = true)
    public File CHART_OUTPUT;

    @Option(doc="If a read maps to a sequence specified with this option, all the bases in the read are counted as ignored bases.  " +
    "These reads are not counted as ")
    public Set<String> IGNORE_SEQUENCE = new HashSet<String>();

    @Option(doc="This percentage of the length of a fragment must overlap one of the ribosomal intervals for a read or read pair by this must in order to be considered rRNA.")
    public double RRNA_FRAGMENT_PERCENTAGE = 0.8;

    @Option(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.  ")
    private Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    private RnaSeqMetricsCollector collector;

    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectRnaSeqMetrics().instanceMainWithExit(argv);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {

        if (CHART_OUTPUT != null) IoUtil.assertFileIsWritable(CHART_OUTPUT);

        final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadRefFlat(REF_FLAT, header.getSequenceDictionary());
        LOG.info("Loaded " + geneOverlapDetector.getAll().size() + " genes.");

        final Long ribosomalBasesInitialValue = RIBOSOMAL_INTERVALS != null ? 0L : null;
        final OverlapDetector<Interval> ribosomalSequenceOverlapDetector = RnaSeqMetricsCollector.makeOverlapDetector(samFile, header, RIBOSOMAL_INTERVALS);

        final HashSet<Integer> ignoredSequenceIndices = RnaSeqMetricsCollector.makeIgnoredSequenceIndicesSet(header, IGNORE_SEQUENCE);

        collector = new RnaSeqMetricsCollector(METRIC_ACCUMULATION_LEVEL, header.getReadGroups(), ribosomalBasesInitialValue,
                geneOverlapDetector, ribosomalSequenceOverlapDetector, ignoredSequenceIndices, MINIMUM_LENGTH, STRAND_SPECIFICITY, RRNA_FRAGMENT_PERCENTAGE);
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
        for (Histogram<Integer> histo : file.getAllHistograms()) {
            atLeastOneHistogram = atLeastOneHistogram || !histo.isEmpty();
        }
        // Generate the coverage by position plot
        if (CHART_OUTPUT != null && atLeastOneHistogram) {
            final int rResult = RExecutor.executeFromClasspath("net/sf/picard/analysis/rnaSeqCoverage.R",
                                                               OUTPUT.getAbsolutePath(),
                                                               CHART_OUTPUT.getAbsolutePath(),
                                                               INPUT.getName());

            if (rResult != 0) {
                throw new PicardException("Problem invoking R to generate plot.");
            }
        }
    }


}
