package picard.analysis.artifacts;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.analysis.CollectOxoGMetrics.*;
import picard.analysis.artifacts.SequencingArtifactMetrics.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

@CommandLineProgramProperties(
        summary = ConvertSequencingArtifactToOxoG.USAGE_SUMMARY + ConvertSequencingArtifactToOxoG.USAGE_DETAILS,
        oneLineSummary = ConvertSequencingArtifactToOxoG.USAGE_SUMMARY,
        programGroup = Metrics.class
)
public class ConvertSequencingArtifactToOxoG extends CommandLineProgram {
static final String USAGE_SUMMARY = "Extract OxoG metrics from generalized artifacts metrics.  ";
static final String USAGE_DETAILS = "<p>This tool extracts 8-oxoguanine (OxoG) artifact metrics from the output of " +
"CollectSequencingArtifactsMetrics (a tool that provides detailed information on a variety of artifacts found in sequencing " +
"libraries) and converts them to the CollectOxoGMetrics tool's output format. This conveniently eliminates the need to run " +
"CollectOxoGMetrics if we already ran CollectSequencingArtifactsMetrics in our pipeline. See the documentation for " +
"<a href='http://broadinstitute.github.io/picard/command-line-overview.html#CollectSequencingArtifactsMetrics'>CollectSequencingArtifactsMetrics</a> "+
"and <a href='http://broadinstitute.github.io/picard/command-line-overview.html#CollectOxoGMetrics'>CollectOxoGMetrics</a> "+
"for additional information on these tools.</p>." +

"<p>Note that only the base of the CollectSequencingArtifactsMetrics output file name is required for the (INPUT_BASE) "+
"parameter. For example, if the file name is artifact_metrics.txt.bait_bias_detail_metrics or "+
"artifact_metrics.txt.pre_adapter_detail_metrics, only the file name base 'artifact_metrics' is " +
"required on the command line for this parameter.  An output file called 'artifact_metrics.oxog_metrics' will be generated "+
"automatically.  Finally, to run this tool successfully, the REFERENCE_SEQUENCE must be provided.</p>"+
"<h4>Usage example:</h4>" +
"<pre>" +
"java -jar picard.jar ConvertSequencingArtifactToOxoG \\<br />" +
"     I=artifact_metrics \\<br />" +
"     R=reference.fasta" +
"</pre>" +
"Please see the metrics definitions page at " +
"<a href='http://broadinstitute.github.io/picard/picard-metric-definitions.html#CollectOxoGMetrics.CpcgMetrics'>ConvertSequencingArtifactToOxoG</a> "+
"for detailed descriptions of the output metrics produced by this tool."+
"<hr />"
;
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Basename of the input artifact metrics file (output by CollectSequencingArtifactMetrics)")
    public File INPUT_BASE;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Basename for output OxoG metrics. Defaults to same basename as input metrics",
            optional = true)
    public File OUTPUT_BASE;

    // Stock main method
    public static void main(final String[] args) {
        new ConvertSequencingArtifactToOxoG().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        if (OUTPUT_BASE == null) { OUTPUT_BASE = INPUT_BASE; }

        final File PRE_ADAPTER_IN = new File(INPUT_BASE + SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT);
        final File BAIT_BIAS_IN = new File(INPUT_BASE + SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT);
        final File OXOG_OUT = new File(OUTPUT_BASE + ".oxog_metrics");

        IOUtil.assertFileIsReadable(PRE_ADAPTER_IN);
        IOUtil.assertFileIsReadable(BAIT_BIAS_IN);
        IOUtil.assertFileIsWritable(OXOG_OUT);

        final List<PreAdapterDetailMetrics> preAdapterDetailMetricsList = MetricsFile.readBeans(PRE_ADAPTER_IN);
        final List<BaitBiasDetailMetrics> baitBiasDetailMetricsList = MetricsFile.readBeans(BAIT_BIAS_IN);

        // TODO should we validate that the two inputs match up as expected?

        /**
         * Determine output fields. Just copy these from the input for now.
         */
        final String oxogSampleAlias = preAdapterDetailMetricsList.get(0).SAMPLE_ALIAS;
        final Set<String> oxogLibraries = new HashSet<String>();
        final Set<String> oxogContexts = new HashSet<String>();
        for (final PreAdapterDetailMetrics preAdapter : preAdapterDetailMetricsList) {
            oxogLibraries.add(preAdapter.LIBRARY);
            // Remember that OxoG only reports on the 'C' contexts
            if (preAdapter.REF_BASE == 'C') {
                oxogContexts.add(preAdapter.CONTEXT);
            }
        }

        /**
         * Store the inputs in maps of {Library -> {Context, Metric}} for easy access.
         * Remember, we only care about two transitions - C>A and G>T! Thus, for each context we
         * will only store one metric.
         */
        final Map<String, Map<String, PreAdapterDetailMetrics>> preAdapterDetailMetricsMap = new HashMap<String, Map<String, PreAdapterDetailMetrics>>();
        final Map<String, Map<String, BaitBiasDetailMetrics>> baitBiasDetailMetricsMap = new HashMap<String, Map<String, BaitBiasDetailMetrics>>();
        for (final String library : oxogLibraries) {
            final Map<String, PreAdapterDetailMetrics> contextsToPreAdapter = new HashMap<String, PreAdapterDetailMetrics>();
            final Map<String, BaitBiasDetailMetrics> contextsToBaitBias = new HashMap<String, BaitBiasDetailMetrics>();
            preAdapterDetailMetricsMap.put(library, contextsToPreAdapter);
            baitBiasDetailMetricsMap.put(library, contextsToBaitBias);
        }
        for (final PreAdapterDetailMetrics preAdapter : preAdapterDetailMetricsList) {
            final Transition transition = Transition.transitionOf(preAdapter.REF_BASE, preAdapter.ALT_BASE);
            if (isOxoG(transition)) {
                preAdapterDetailMetricsMap.get(preAdapter.LIBRARY).put(preAdapter.CONTEXT, preAdapter);
            }
        }
        for (final BaitBiasDetailMetrics baitBias : baitBiasDetailMetricsList) {
            final Transition transition = Transition.transitionOf(baitBias.REF_BASE, baitBias.ALT_BASE);
            if (isOxoG(transition)) {
                baitBiasDetailMetricsMap.get(baitBias.LIBRARY).put(baitBias.CONTEXT, baitBias);
            }
        }

        /**
         * Create the OxoG metrics
         */
        final List<CpcgMetrics> oxogMetrics = new ArrayList<CpcgMetrics>();
        for (final String library : oxogLibraries) {
            for (final String context : oxogContexts) {
                final CpcgMetrics m = new CpcgMetrics();
                m.SAMPLE_ALIAS = oxogSampleAlias;
                m.LIBRARY = library;
                m.CONTEXT = context;
                m.TOTAL_SITES = 0; // not calculated in the input metrics

                /**
                 * Get the relevant input metrics. This is done in a somewhat confusing way:
                 *
                 * 1. For pre-adapter metrics: note that OxoG only reports 'C' contexts, even though the actual pre-adapter OxoG artifact
                 *    occurs when the reference-strand base is 'G'. This is because OxoG reverse-complements all the contexts for some reason.
                 *    Thus when we add an entry for 'ACA' in the output, we actually need to get that data from 'TGT' in the input.
                 *
                 * 2. For bait-bias metrics: for each context, we report two opposing error rates, C_REF and G_REF, because for this metric
                 *    the bias could really go in either direction (whereas for pre-adapter artifacts we only expect one direction: G>T, but
                 *    never C>A, on the original reference strand). C_REF corresponds to the actual context printed in that row, and G_REF
                 *    corresponds to its reverse complement. So for 'ACA' in the output, we need to take data from both 'ACA' and 'TGT' in the
                 *    input.
                 */

                final PreAdapterDetailMetrics preAdapter = preAdapterDetailMetricsMap.get(library).get(SequenceUtil.reverseComplement(context));
                final BaitBiasDetailMetrics baitBiasFwd = baitBiasDetailMetricsMap.get(library).get(context);
                final BaitBiasDetailMetrics baitBiasRev = baitBiasDetailMetricsMap.get(library).get(SequenceUtil.reverseComplement(context));

                // extract fields from PreAdapterDetailMetrics
                m.TOTAL_BASES = preAdapter.PRO_REF_BASES + preAdapter.PRO_ALT_BASES + preAdapter.CON_REF_BASES + preAdapter.CON_ALT_BASES;
                m.REF_TOTAL_BASES = preAdapter.PRO_REF_BASES + preAdapter.CON_REF_BASES;
                m.REF_NONOXO_BASES = preAdapter.CON_REF_BASES;
                m.REF_OXO_BASES = preAdapter.PRO_REF_BASES;
                m.ALT_NONOXO_BASES = preAdapter.CON_ALT_BASES;
                m.ALT_OXO_BASES = preAdapter.PRO_ALT_BASES;
                m.OXIDATION_ERROR_RATE = preAdapter.ERROR_RATE;
                m.OXIDATION_Q = preAdapter.QSCORE;

                // extract fields from BaitBiasDetailMetrics
                m.C_REF_REF_BASES = baitBiasFwd.FWD_CXT_REF_BASES;
                m.G_REF_REF_BASES = baitBiasFwd.REV_CXT_REF_BASES;
                m.C_REF_ALT_BASES = baitBiasFwd.FWD_CXT_ALT_BASES;
                m.G_REF_ALT_BASES = baitBiasFwd.REV_CXT_ALT_BASES;
                m.C_REF_OXO_ERROR_RATE = baitBiasFwd.ERROR_RATE;
                m.C_REF_OXO_Q = baitBiasFwd.QSCORE;
                m.G_REF_OXO_ERROR_RATE = baitBiasRev.ERROR_RATE;
                m.G_REF_OXO_Q = baitBiasRev.QSCORE;

                // add it
                oxogMetrics.add(m);
            }
        }

        final MetricsFile<CpcgMetrics, Integer> outputFile = getMetricsFile();
        for (final CpcgMetrics m : oxogMetrics) {
            outputFile.addMetric(m);
        }

        outputFile.write(OXOG_OUT);
        return 0;
    }

    private boolean isOxoG(final Transition t) {
        return t.equals(Transition.CtoA) || t.equals(Transition.GtoT);
    }
}
