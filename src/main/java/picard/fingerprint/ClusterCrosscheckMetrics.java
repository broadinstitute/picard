
/*
 * The MIT License
 *
 * Copyright (c) 2017 The Broad Institute
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

package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.DiagnosticsAndQCProgramGroup;
import picard.util.GraphUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * <h3>Summary</h3>
 * Clusters the results from a {@link CrosscheckFingerprints} run according to the LOD score. The resulting metric file
 * can be used to assist diagnosing results from {@link CrosscheckFingerprints}. It clusters the connectivity graph between the
 * different groups. Two groups are connected if they have a LOD score greater than the {@link #LOD_THRESHOLD}.
 * <br/>
 *
 * <h3>Details</h3>
 * The results of running {@link CrosscheckFingerprints} can be difficult to analyze, especially when many groups are
 * related (meaning LOD greater than {@link #LOD_THRESHOLD}) in non-transitive manner (A is related to B, B is related to C,
 * but A doesn't seem to be related to C.) {@link ClusterCrosscheckMetrics} clusters the metrics from {@link CrosscheckFingerprints}
 * so that all the groups in a cluster are related to each other either directly, or indirectly (thus A, B and C would
 * end up in one cluster.) Two samples can only be in two different clusters if all the samples from these two clusters
 * do not get high LOD scores when compared to each other.
 * <br/>
 * <h3>Example</h3>
 * <pre>
 *     java -jar picard.jar ClusterCrosscheckMetrics \
 *              INPUT=sample.crosscheck_metrics \
 *              LOD_THRESHOLD=3 \
 *              OUTPUT=sample.clustered.crosscheck_metrics
 * </pre>
 *
 * The resulting file, consists of the {@link ClusteredCrosscheckMetric} class and contains the original crosscheck metric
 * values, for groups that end-up in the same clusters (regardless of LOD score of each comparison). In addition it notes
 * the {@link ClusteredCrosscheckMetric#CLUSTER} identifier and the size of the cluster (in {@link ClusteredCrosscheckMetric#CLUSTER_SIZE}.)
 * Groups that do not have high LOD scores with any other group (including itself!) will not be included in the metric file.
 * Note that cross-group comparisons are not included in the metric file.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        summary =
                "<h3>Summary</h3>\n" +
                "Clusters the results from a {@link CrosscheckFingerprints} run according to the LOD score. The resulting metric file " +
                "can be used to assist diagnosing results from {@link CrosscheckFingerprints}. It clusters the connectivity graph between the " +
                "different groups. Two groups are connected if they have a LOD score greater than the {@link #LOD_THRESHOLD}.\n " +
                "\n" +
                "<h3>Details</h3>" +
                "The results of running {@link CrosscheckFingerprints} can be difficult to analyze, especially when many groups are " +
                "related (meaning LOD greater than {@link #LOD_THRESHOLD}) in non-transitive manner (A is related to B, B is related to C, " +
                "but A doesn't seem to be related to C.) {@link ClusterCrosscheckMetrics} clusters the metrics from {@link CrosscheckFingerprints} " +
                "so that all the groups in a cluster are related to each other either directly, or indirectly (thus A, B and C would " +
                "end up in one cluster.) Two samples can only be in two different clusters if all the samples from these two clusters " +
                "do not get high LOD scores when compared to each other.\n" +
                "\n" +
                "<h3>Example</h3>\n" +
                "<pre>\n" +
                "    java -jar picard.jar ClusterCrosscheckMetrics \\\n" +
                "             INPUT=sample.crosscheck_metrics \\\n" +
                "             LOD_THRESHOLD=3 \\\n" +
                "             OUTPUT=sample.clustered.crosscheck_metrics\n" +
                "</pre>\n" +
                "\n" +
                "The resulting file, consists of the {@link ClusteredCrosscheckMetric} class and contains the original crosscheck metric " +
                "values, for groups that end-up in the same clusters (regardless of LOD score of each comparison). In addition it notes " +
                "the {@link ClusteredCrosscheckMetric#CLUSTER} identifier and the size of the cluster (in {@link ClusteredCrosscheckMetric#CLUSTER_SIZE}.) " +
                "Groups that do not have high LOD scores with any other group (including " +
                "itself!) will not be included in the metric file. Note that cross-group comparisons are not included in the metric file. ",
        oneLineSummary = "Clusters the results of a CrosscheckFingerprints run by LOD score",
        programGroup = DiagnosticsAndQCProgramGroup.class

)
@DocumentedFeature
public class ClusterCrosscheckMetrics extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "The cross-check metrics file to be clustered.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true,
            doc = "Output file to write metrics to. Will write to stdout if null.")
    public File OUTPUT;

    @Argument(shortName = "LOD",
            doc = "LOD score to be used as the threshold for clustering.")
    public double LOD_THRESHOLD = 0;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        if(OUTPUT != null) IOUtil.assertFileIsWritable(OUTPUT);

        final MetricsFile<CrosscheckMetric, ?> metricsFile = getMetricsFile();

        try {
            metricsFile.read(new FileReader(INPUT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return 1;
        }

        clusterMetrics(metricsFile.getMetrics()).write(OUTPUT);

        return 0;
    }

    private MetricsFile<ClusteredCrosscheckMetric, ?> clusterMetrics(final List<CrosscheckMetric> metrics) {
        final GraphUtils.Graph<String> graph = new GraphUtils.Graph<>();
        metrics.stream()
                .filter(metric -> metric.LOD_SCORE > LOD_THRESHOLD)
                .forEach(metric -> {
                    final String lhsBy = metric.LEFT_GROUP_VALUE;
                    final String rhsBy = metric.RIGHT_GROUP_VALUE;

                    graph.addEdge(lhsBy, rhsBy);
                });

        final Map<String, Integer> clusters = graph.cluster();

        // invert map...get map from group integer to list of group_value
        final Map<Integer, Set<String>> collection = clusters.entrySet().stream()
                .collect(Collectors.groupingBy(Map.Entry::getValue))
                .entrySet()
                .stream()
                .collect(Collectors
                        .toMap(Map.Entry::getKey, entry -> entry.getValue()
                                .stream()
                                .map(Map.Entry::getKey)
                                .collect(Collectors.toSet())));

        final MetricsFile<ClusteredCrosscheckMetric, ?> clusteredMetrics = getMetricsFile();
        // for each cluster, find the metrics that compare groups that are both from the cluster
        // and add them to the metrics file
        for (final Map.Entry<Integer, Set<String>> cluster : collection.entrySet()) {

            clusteredMetrics.addAllMetrics(
                    metrics.stream()
                            .filter(metric ->
                                    cluster.getValue().contains(metric.LEFT_GROUP_VALUE) &&
                                            cluster.getValue().contains(metric.RIGHT_GROUP_VALUE))
                            .map(metric -> {
                                final ClusteredCrosscheckMetric clusteredCrosscheckMetric = new ClusteredCrosscheckMetric(metric);
                                clusteredCrosscheckMetric.CLUSTER = cluster.getKey();
                                clusteredCrosscheckMetric.CLUSTER_SIZE = cluster.getValue().size();

                                return clusteredCrosscheckMetric;
                            })
                            .collect(Collectors.toSet()));
        }
        return clusteredMetrics;
    }
}
