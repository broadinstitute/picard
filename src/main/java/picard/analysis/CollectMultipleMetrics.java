/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.analysis.artifacts.CollectSequencingArtifactMetrics;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.programgroups.Metrics;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.*;

/**
 * Class that is designed to instantiate and execute multiple metrics programs that extend
 * SinglePassSamProgram while making only a single pass through the SAM file and supplying
 * each program with the records as it goes.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(

        summary = CollectMultipleMetrics.USAGE_SUMMARY + CollectMultipleMetrics.USAGE_DETAILS,
        oneLineSummary = CollectMultipleMetrics.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class CollectMultipleMetrics extends CommandLineProgram {

    /**
     * This interface allows developers to create Programs to run in addition to the ones defined in the Program enum.
     * Includes a method for determining whether or not a Program explicitly needs a reference sequence (i.e. cannot be null)
     */

    static final String USAGE_SUMMARY ="Collect multiple classes of metrics.  ";
    static final String USAGE_DETAILS ="This 'meta-metrics' tool runs one or more of the metrics collection modules at the same" +
            " time to cut down on the time spent reading in data from input files. Available modules include " +
            "CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution,  MeanQualityByCycle, " +
            "CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, and CollectQualityYieldMetrics. " +
            "The tool produces outputs of '.pdf' and '.txt' files for each module, except for the " +
            "CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file. Output files are named by specifying a base name " +
            "(without any file extensions).<br /><br />" +
            "" +
            "<p>Currently all programs are run with default options and fixed output extensions, " +
            "but this may become more flexible in future. Specifying a reference sequence file is required.</p>" +

            "<p>Note: Metrics labeled as percentages are actually expressed as fractions!</p>" +
            "" +
            "<h4>Usage example (all modules on by default):</h4>" +
            "<pre>" +
            "java -jar picard.jar CollectMultipleMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=multiple_metrics \\<br />" +
            "      R=reference_sequence.fasta <br />" +
            "</pre>" +
            "<h4>Usage example (two modules only):</h4>" +
            "java -jar picard.jar CollectMultipleMetrics \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=multiple_metrics \\<br />" +
            "      R=reference_sequence.fasta \\<br />" +
            "      PROGRAM=null \\<br />" +
            "      PROGRAM=QualityScoreDistribution \\<br />" +
            "      PROGRAM=MeanQualityByCycle "+
            "</pre>" +
            "<hr />";
    public static interface ProgramInterface {
        /** By default, this method calls the {@link #makeInstance(String, String, File, File, Set, File, File)} method without 'includeUnpaired' parameter. */
        default  SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference,
                                                   final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals, final boolean includeUnpaired) {
            return makeInstance(outbase, outext, input, reference, metricAccumulationLevel, dbSnp, intervals);
        }
        SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference,
                                          final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals);
        public boolean needsReferenceSequence();
        public boolean supportsMetricAccumulationLevel();
    }

    public static enum Program implements ProgramInterface {
        CollectAlignmentSummaryMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final CollectAlignmentSummaryMetrics program = new CollectAlignmentSummaryMetrics();
                program.OUTPUT = new File(outbase + ".alignment_summary_metrics" + outext);

                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.METRIC_ACCUMULATION_LEVEL = metricAccumulationLevel;
                program.INPUT = input;
                program.setReferenceSequence(reference);

                return program;
            }
        },
        CollectInsertSizeMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final CollectInsertSizeMetrics program = new CollectInsertSizeMetrics();
                program.OUTPUT = new File(outbase + ".insert_size_metrics" + outext);
                program.HISTOGRAM_FILE = new File(outbase + ".insert_size_histogram.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.METRIC_ACCUMULATION_LEVEL = metricAccumulationLevel;
                program.INPUT = input;
                program.setReferenceSequence(reference);

                return program;
            }
        },
        QualityScoreDistribution {
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return false;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final QualityScoreDistribution program = new QualityScoreDistribution();
                program.OUTPUT = new File(outbase + ".quality_distribution_metrics" + outext);
                program.CHART_OUTPUT = new File(outbase + ".quality_distribution.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.setReferenceSequence(reference);

                return program;
            }
        },
        MeanQualityByCycle {
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return false;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final MeanQualityByCycle program = new MeanQualityByCycle();
                program.OUTPUT = new File(outbase + ".quality_by_cycle_metrics" + outext);
                program.CHART_OUTPUT = new File(outbase + ".quality_by_cycle.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.setReferenceSequence(reference);

                return program;
            }
        },
        CollectBaseDistributionByCycle {
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return false;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final CollectBaseDistributionByCycle program = new CollectBaseDistributionByCycle();
                program.OUTPUT = new File(outbase + ".base_distribution_by_cycle_metrics" + outext);
                program.CHART_OUTPUT = new File(outbase + ".base_distribution_by_cycle.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.setReferenceSequence(reference);

                return program;
            }
        },
        CollectGcBiasMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return true;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final CollectGcBiasMetrics program = new CollectGcBiasMetrics();
                program.OUTPUT = new File(outbase + ".gc_bias.detail_metrics" + outext);
                program.SUMMARY_OUTPUT = new File(outbase + ".gc_bias.summary_metrics" + outext);
                program.CHART_OUTPUT = new File(outbase + ".gc_bias.pdf");

                program.INPUT = input;
                // previously MetricAccumulationLevel.ALL_READS, MetricAccumulationLevel.LIBRARY
                program.METRIC_ACCUMULATION_LEVEL = metricAccumulationLevel;
                program.SCAN_WINDOW_SIZE = 100;
                program.MINIMUM_GENOME_FRACTION = 1.0E-5;
                program.IS_BISULFITE_SEQUENCED = false;
                program.ASSUME_SORTED = false;
                program.ALSO_IGNORE_DUPLICATES = false;

                //GC_Bias actually uses the class-level REFERENCE_SEQUENCE variable.
                program.setReferenceSequence(reference);

                return program;
            }
        },
        RnaSeqMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return true;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final CollectRnaSeqMetrics program = new CollectRnaSeqMetrics();
                program.OUTPUT       = new File(outbase + ".rna_metrics" + outext);
                program.CHART_OUTPUT = new File(outbase + ".rna_coverage.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.METRIC_ACCUMULATION_LEVEL = metricAccumulationLevel;
                program.INPUT = input;
                program.setReferenceSequence(reference);
                
                return program;
            }
        },
        CollectSequencingArtifactMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return true;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() { return false; }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                return makeInstance(outbase, outext, input, reference, metricAccumulationLevel, dbSnp, intervals, false);
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals, final boolean includeUnpaired) {
                final CollectSequencingArtifactMetrics program = new CollectSequencingArtifactMetrics();
                program.OUTPUT = new File(outbase);
                program.FILE_EXTENSION = outext;
                program.DB_SNP = dbSnp;
                program.INTERVALS = intervals;
                program.INCLUDE_UNPAIRED = includeUnpaired;
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.setReferenceSequence(reference);
                return program;
            }
        },
        CollectQualityYieldMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return false;
            }
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return false;
            }
            @Override
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input, final File reference, final Set<MetricAccumulationLevel> metricAccumulationLevel, final File dbSnp, final File intervals) {
                final CollectQualityYieldMetrics program = new CollectQualityYieldMetrics();
                program.OUTPUT = new File(outbase + ".quality_yield_metrics" + outext);
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;
                program.setReferenceSequence(reference);
                return program;
            }
        }
    }

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM file.")
    public File INPUT;

    @Argument(doc = "If true (default), then the sort order in the header file will be ignored.",
            shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME)
    public boolean ASSUME_SORTED = true;

    @Argument(doc = "Stop after processing N reads, mainly for debugging.")
    public int STOP_AFTER = 0;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Base name of output files.")
    public String OUTPUT;

    // create the default accumulation level as a variable. We'll use this to init the command-line arg and for validation later.
    private final Set<MetricAccumulationLevel> accumLevelDefault = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(shortName="LEVEL", doc="The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = new HashSet<>(accumLevelDefault);

    @Argument(shortName = "EXT", doc="Append the given file extension to all metric file names (ex. OUTPUT.insert_size_metrics.EXT). None if null", optional=true)
    public String FILE_EXTENSION = null;

    @Argument(doc = "Set of metrics programs to apply during the pass through the SAM file.")
    public Set<Program> PROGRAM = new LinkedHashSet<>(Arrays.asList(Program.CollectAlignmentSummaryMetrics, Program.CollectBaseDistributionByCycle,
            Program.CollectInsertSizeMetrics, Program.MeanQualityByCycle, Program.QualityScoreDistribution));

    @Argument(doc = "An optional list of intervals to restrict analysis to. Only pertains to some of the PROGRAMs. Programs whose stand-alone CLP does not " +
            "have an INTERVALS argument will silently ignore this argument.", optional = true)
    public File INTERVALS;

    @Argument(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis " +
            "by some PROGRAMs; PROGRAMs whose CLP doesn't allow for this argument will quietly ignore it.", optional = true)
    public File DB_SNP;

    @Argument(shortName = "UNPAIRED", doc = "Include unpaired reads in CollectSequencingArtifactMetrics. If set to true then all paired reads will be included as well - " +
            "MINIMUM_INSERT_SIZE and MAXIMUM_INSERT_SIZE will be ignored in CollectSequencingArtifactMetrics.")
    public boolean INCLUDE_UNPAIRED = false;
    /**
     * Contents of PROGRAM set is transferred to this set during command-line validation, so that an outside
     * developer can invoke this class programmatically and provide alternative Programs to run by calling
     * setProgramsToRun().
     */
    private Set<ProgramInterface> programsToRun;

    private static final Log log = Log.getInstance(CollectMultipleMetrics.class);

    // Stock main method
    public static void main(final String[] args) {
        new CollectMultipleMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (PROGRAM.isEmpty()) {
            return new String[]{"No programs specified with PROGRAM"};
        }
        programsToRun = new LinkedHashSet<>(PROGRAM);

        return super.customCommandLineValidation();
    }

    /**
     * Use this method when invoking CollectMultipleMetrics programmatically to run programs other than the ones
     * available via enum.  This must be called before doWork().
     */
    public void setProgramsToRun(final Collection<ProgramInterface> programsToRun) {
        this.programsToRun.clear();
        this.programsToRun.addAll(programsToRun);
    }

    @Override
    public int doWork() {
        if (OUTPUT.endsWith(".")) {
            OUTPUT = OUTPUT.substring(0, OUTPUT.length() - 1);
        }

        final List<SinglePassSamProgram> programs = new ArrayList<>();
        for (final ProgramInterface program : programsToRun) {
            if (program.needsReferenceSequence() && REFERENCE_SEQUENCE == null) {
                throw new PicardException("The " + program.toString() + " program needs a Reference Sequence, please set REFERENCE_SEQUENCE in the command line");
            }
            if (!accumLevelDefault.equals(METRIC_ACCUMULATION_LEVEL) && !program.supportsMetricAccumulationLevel()) {
                log.warn("The " + program.toString() + " program does not support a metric accumulation level, but METRIC_ACCUMULATION_LEVEL" +
                        " was overridden in the command line. " + program.toString() + " will be run against the entire input.");
            }

            final String outext = (null != FILE_EXTENSION) ? FILE_EXTENSION : ""; // Add a file extension if desired
            final SinglePassSamProgram instance = program.makeInstance(OUTPUT, outext, INPUT, REFERENCE_SEQUENCE, METRIC_ACCUMULATION_LEVEL, DB_SNP, INTERVALS, INCLUDE_UNPAIRED);

            // Generally programs should not be accessing these directly but it might make things smoother
            // to just set them anyway
            instance.INPUT = INPUT;
            instance.setReferenceSequence(REFERENCE_SEQUENCE);

            instance.setDefaultHeaders(getDefaultHeaders());

            programs.add(instance);
        }
        SinglePassSamProgram.makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, programs);

        return 0;
    }
}
