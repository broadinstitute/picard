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
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineParser;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.PicardException;
import picard.analysis.artifacts.CollectSequencingArtifactMetrics;
import picard.analysis.directed.RnaSeqMetricsCollector;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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
        programGroup = DiagnosticsAndQCProgramGroup.class
)
@DocumentedFeature
public class CollectMultipleMetrics extends CommandLineProgram {

    /**
     * This interface allows developers to create Programs to run in addition to the ones defined in the Program enum.
     * Includes a method for determining whether or not a Program explicitly needs a reference sequence (i.e. cannot be null)
     */

    static final String USAGE_SUMMARY = "Collect multiple classes of metrics.";
    static final String USAGE_DETAILS = "This 'meta-metrics' tool runs one or more of the metrics collection modules at the same" +
            " time to cut down on the time spent reading in data from input files. Available modules include " +
            "CollectAlignmentSummaryMetrics, CollectInsertSizeMetrics, QualityScoreDistribution,  MeanQualityByCycle, " +
            "CollectBaseDistributionByCycle, CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics" +
            " and CollectQualityYieldMetrics. " +
            "The tool produces outputs of '.pdf' and '.txt' files for each module, except for the " +
            "CollectAlignmentSummaryMetrics module, which outputs only a '.txt' file." +
            " Output files are named by specifying a base name (without any file extensions).<br /><br />" +
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
            "      PROGRAM=MeanQualityByCycle " +
            "</pre>" +
            "<hr />";

    public interface ProgramInterface {

        /**
         * By default, this method calls the
         * {@link #makeInstance(String, String, File, File, Set, File, File, File, Set)} method
         * without 'includeUnpaired' parameter.
         */
        default SinglePassSamProgram makeInstance(final String outbase,
                                                  final String outext,
                                                  final File input,
                                                  final File reference,
                                                  final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                  final File dbSnp, final File intervals,
                                                  final File refflat, Set<String> ignoreSequence,
                                                  final boolean includeUnpaired) {

            return makeInstance(outbase, outext, input,
                    reference,
                    metricAccumulationLevel,
                    dbSnp,
                    intervals,
                    refflat,
                    ignoreSequence);
        }

        SinglePassSamProgram makeInstance(final String outbase,
                                          final String outext,
                                          final File input,
                                          final File reference,
                                          final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                          final File dbSnp,
                                          final File intervals,
                                          final File refflat,
                                          final Set<String> ignoreSequence);

        default boolean needsReferenceSequence() {
            return false;
        }

        default boolean needsRefflatFile() {
            return false;
        }

        default boolean supportsMetricAccumulationLevel() {
            return false;
        }
    }

    public enum Program implements ProgramInterface {
        CollectAlignmentSummaryMetrics {
            @Override
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }

            @Override
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
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
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }

            @Override
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
                final CollectInsertSizeMetrics program = new CollectInsertSizeMetrics();
                program.OUTPUT = new File(outbase + ".insert_size_metrics" + outext);
                program.Histogram_FILE = new File(outbase + ".insert_size_histogram.pdf");
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
            public SinglePassSamProgram makeInstance(final String outbase, final String outext, final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
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
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
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
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
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
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
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
            public boolean needsRefflatFile() {
                return true;
            }

            @Override
            public boolean supportsMetricAccumulationLevel() {
                return true;
            }

            @Override
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
                final CollectRnaSeqMetrics program = new CollectRnaSeqMetrics();
                program.OUTPUT = new File(outbase + ".rna_metrics" + outext);
                program.CHART_OUTPUT = new File(outbase + ".rna_coverage.pdf");
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.METRIC_ACCUMULATION_LEVEL = metricAccumulationLevel;
                program.INPUT = input;
                program.RIBOSOMAL_INTERVALS = intervals;
                program.IGNORE_SEQUENCE = ignoreSequence;
                program.REF_FLAT = refflat;
                program.STRAND_SPECIFICITY = RnaSeqMetricsCollector.StrandSpecificity.SECOND_READ_TRANSCRIPTION_STRAND;

                return program;
            }
        },

        CollectSequencingArtifactMetrics {
            @Override
            public boolean needsReferenceSequence() {
                return true;
            }

            @Override
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {

                return makeInstance(outbase, outext, input,
                        reference,
                        metricAccumulationLevel,
                        dbSnp,
                        intervals,
                        refflat,
                        ignoreSequence,
                        false);
            }

            @Override
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence,
                                                     final boolean includeUnpaired) {
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
            public SinglePassSamProgram makeInstance(final String outbase,
                                                     final String outext,
                                                     final File input,
                                                     final File reference,
                                                     final Set<MetricAccumulationLevel> metricAccumulationLevel,
                                                     final File dbSnp,
                                                     final File intervals,
                                                     final File refflat,
                                                     final Set<String> ignoreSequence) {
                final CollectQualityYieldMetrics program = new CollectQualityYieldMetrics();
                program.OUTPUT = new File(outbase + ".quality_yield_metrics" + outext);
                // Generally programs should not be accessing these directly but it might make things smoother
                // to just set them anyway. These are set here to make sure that in case of a the derived class
                // overrides
                program.INPUT = input;

                return program;
            }
        }
    }

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input SAM or BAM file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.ASSUME_SORTED_SHORT_NAME,
            doc = "If true (default), then the sort order in the header file will be ignored.")
    public boolean ASSUME_SORTED = true;

    @Argument(doc = "Stop after processing N reads, mainly for debugging.")
    public int STOP_AFTER = 0;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Base name of output files.")
    public String OUTPUT;

    //create the default accumulation level as a variable.
    // We'll use this to init the command-line arg and for validation later.
    private final Set<MetricAccumulationLevel> accumLevelDefault = CollectionUtil.makeSet(MetricAccumulationLevel.ALL_READS);

    @Argument(shortName = "LEVEL",
            doc = "The level(s) at which to accumulate metrics.")
    public Set<MetricAccumulationLevel> METRIC_ACCUMULATION_LEVEL = new HashSet<>(accumLevelDefault);

    @Argument(shortName = "EXT",
            doc = "Append the given file extension to all metric file names (ex. OUTPUT.insert_size_metrics.EXT). None if null",
            optional = true)
    public String FILE_EXTENSION = null;

    @Argument(doc = "Set of metrics programs to apply during the pass through the SAM file.")
    public Set<Program> PROGRAM = new LinkedHashSet<>(Arrays.asList(
            Program.CollectAlignmentSummaryMetrics,
            Program.CollectBaseDistributionByCycle,
            Program.CollectInsertSizeMetrics,
            Program.MeanQualityByCycle,
            Program.QualityScoreDistribution
    ));

    @Argument(doc = "An optional list of intervals to restrict analysis to. Only pertains to some of the PROGRAMs. " +
            "Programs whose stand-alone CLP does not have an INTERVALS argument will silently ignore this argument.",
            optional = true)
    public File INTERVALS;

    @Argument(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis " +
            "by some PROGRAMs; PROGRAMs whose CLP doesn't allow for this argument will quietly ignore it.",
            optional = true)
    public File DB_SNP;

    @Argument(doc = "Gene annotations in refFlat form.  " +
            "Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat",
            optional = true)
    public File REF_FLAT;

    @Argument(doc = "If a read maps to a sequence specified with this option, " +
            "all the bases in the read are counted as ignored bases.",
            optional = true)
    public Set<String> IGNORE_SEQUENCE = new HashSet<>();

    @Argument(shortName = "UNPAIRED",
            doc = "Include unpaired reads in CollectSequencingArtifactMetrics. " +
                    "If set to true then all paired reads will be included as well - " +
                    "MINIMUM_INSERT_SIZE and MAXIMUM_INSERT_SIZE will be ignored in CollectSequencingArtifactMetrics.")
    public boolean INCLUDE_UNPAIRED = false;

    @Argument(doc="extra arguments to the various tools can be specified using the following format:" +
            "<PROGRAM>::<ARGUMENT_AND_VALUE> where <PROGRAM> is one of the programs specified in PROGRAM, " +
            "and <ARGUMENT_AND_VALUE> are the argument and value that you'd like to specify as you would on the command line. " +
            "For example, to change the HISTOGRAM_WIDTH in CollectInsertSizeMetrics to 200, use:\n " +
            "\"EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200\"\n " +
            "or, in the new parser:" +
            "--EXTRA_ARGUMENT \"CollectInsertSizeMetrics::--HISTOGRAM_WIDTH 200\"\n " +
            "(Quotes are required to avoid the shell from separating this into two arguments.) " +
            "Note that the following arguments cannot be modified on a per-program level: INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, and STOP_AFTER. " +
            "Providing them in an EXTRA_ARGUMENT will _not_ result in an error, but they will be silently ignored. " , optional = true)
    public List<String> EXTRA_ARGUMENT = null;

    /**
     * Contents of PROGRAM set is transferred to this set during command-line validation, so that an outside
     * developer can invoke this class programmatically and provide alternative Programs to run by calling
     * setProgramsToRun().
     */
    private Set<ProgramInterface> programsToRun;

    private static final Log log = Log.getInstance(CollectMultipleMetrics.class);

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

        final Map<ProgramInterface, List<String>> additionalArguments = processAdditionalArguments(EXTRA_ARGUMENT);

        final List<SinglePassSamProgram> programs = new ArrayList<>();
        for (final ProgramInterface program : programsToRun) {
            if (program.needsReferenceSequence() && REFERENCE_SEQUENCE == null) {
                throw new PicardException("The " + program.toString() + " program needs a REF Sequence, " +
                        "please set REFERENCE_SEQUENCE in the command line");
            }
            if (program.needsRefflatFile() && REF_FLAT == null) {
                throw new PicardException("The " + program.toString() + " program needs a gene annotations file, " +
                        "please set REF_FLAT in the command line");
            }
            if (!accumLevelDefault.equals(METRIC_ACCUMULATION_LEVEL) && !program.supportsMetricAccumulationLevel()) {
                log.warn("The " + program.toString() + " program does not support a metric accumulation level, " +
                        "but METRIC_ACCUMULATION_LEVEL was overridden in the command line. " +
                        program.toString() + " will be run against the entire input.");
            }
            final String outext = (null != FILE_EXTENSION) ? FILE_EXTENSION : ""; // Add a file extension if desired
            final SinglePassSamProgram instance = program.makeInstance(OUTPUT, outext, INPUT,
                    REFERENCE_SEQUENCE,
                    METRIC_ACCUMULATION_LEVEL,
                    DB_SNP,
                    INTERVALS,
                    REF_FLAT,
                    IGNORE_SEQUENCE,
                    INCLUDE_UNPAIRED);

            if (additionalArguments.containsKey(program)) {
                final CommandLineParser commandLineParser = getCommandLineParser(instance);
                final boolean success = commandLineParser.parseArguments(System.err, additionalArguments.get(program).toArray(new String[0]));
                if (!success) {
                    throw new CommandLineException("Failed to parse arguments ["+ String.join(",", additionalArguments.get(program))+ "] for " + program);
                }
                additionalArguments.remove(program);
            }

            instance.setDefaultHeaders(getDefaultHeaders());
            programs.add(instance);
        }
        if (!additionalArguments.isEmpty()) {
            throw new CommandLineException("EXTRA_ARGUMENT values were provided, but corresponding PROGRAM wasn't requested:" +
                    additionalArguments.entrySet().stream().map(e -> e.getKey().toString() + "::" + e.getValue().toString()).collect(Collectors.joining()));
        }
        SinglePassSamProgram.makeItSo(INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED, STOP_AFTER, programs);

        return 0;
    }

    private static Map<ProgramInterface, List<String>> processAdditionalArguments(final List<String> arguments) {
        final Map<ProgramInterface, List<String>> map = new HashMap<>();
        final Pattern pattern = Pattern.compile("(?<program>.*)::(?<argumentAndValue>.+?)( +(?<optionalValue>.+))?");
        for (String str : arguments) {
            final Matcher matcher = pattern.matcher(str);
            if (!matcher.matches()) {
                throw new CommandLineException("couldn't understand EXTRA_ARGUMENT " + str +
                        " it doesn't conform to the form '<PROGRAM>::<ARGUMENT_AND_VALUE>'.");
            }
            final String programName = matcher.group("program");
            final ProgramInterface program;
            try {
                program = Program.valueOf(programName);
            } catch (IllegalArgumentException e) {
                throw new CommandLineException("Couldn't find program with value " + programName, e);
            }
            final String argumentAndValue = matcher.group("argumentAndValue");

            map.computeIfAbsent(program, (k) -> new ArrayList<>()).add(argumentAndValue);

            final String optionalValue = matcher.group("optionalValue");
            if (!StringUtils.isEmpty(optionalValue)) {
                map.get(program).add(optionalValue);
            }
        }
        return map;
    }
}
