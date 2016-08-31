
/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Fingerprinting;
import picard.util.GraphUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * Program to check that all (read-)groups within the set of input files appear to come from the same
 * individual. Can be used to cross-check libraries, samples, or files.
 *
 * @author Yossi Farjoun
 */
@CommandLineProgramProperties(
        usage = "Clusters the results from a CrosscheckReadGroupFingerprints into groups that are connected according" +
                "to a large enough LOD score.",
        usageShort = "Clusters the results of a CrosscheckReadGroupFingerprints run by LOD score.",
        programGroup = Fingerprinting.class
)
public class ClusterCrosscheckMetrics extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "The cross-check metrics file to be clustered")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, optional = true,
            doc = "Optional output file to write metrics to. Default is to write to stdout.")
    public File OUTPUT;

    @Option(shortName = "LOD",
            doc = "LOD score to be used as the threshold for clustering.")
    public double LOD_THRESHOLD = 0;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        MetricsFile<CrosscheckMetric, ?> metricsFile = getMetricsFile();

        try {
            metricsFile.read(new FileReader(INPUT));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return 1;
        }

        clusterMetrics(metricsFile.getMetrics()).write(OUTPUT);

        return 0;
    }

    private MetricsFile<ClusteredCrosscheckMetric, ?> clusterMetrics(List<CrosscheckMetric> metrics) {
        final GraphUtils.Graph<String> graph = new GraphUtils.Graph<>();
        metrics.stream()
                .filter(metric -> metric.LOD_SCORE > LOD_THRESHOLD)
                .forEach(metric -> {
                    final String lhsBy = metric.LEFT_GROUP_VALUE;
                    final String rhsBy = metric.RIGHT_GROUP_VALUE;

                    graph.addEdge(lhsBy, rhsBy);
                });

        final Map<String, Integer> clusters = graph.cluster();

        //invert map...get map from group integer to list of group_value
        final Map<Integer, Set<String>> collection = clusters.entrySet().stream()
                .collect(Collectors.groupingBy(Map.Entry::getValue))
                .entrySet()
                .stream()
                .collect(Collectors
                        .toMap(Map.Entry::getKey, entry -> entry.getValue()
                                .stream()
                                .map(Map.Entry::getKey)
                                .collect(Collectors.toSet())));

        MetricsFile<ClusteredCrosscheckMetric, ?> clusteredMetrics = getMetricsFile();
        // for each cluster, find the metrics that compare groups that are both from the cluster
        // and add them to the metrics file
        for (Map.Entry<Integer, Set<String>> cluster : collection.entrySet()) {

            clusteredMetrics.addAllMetrics(
                    metrics.stream()
                            .filter(metric ->
                                    cluster.getValue().contains(metric.LEFT_GROUP_VALUE) &&
                                            cluster.getValue().contains(metric.RIGHT_GROUP_VALUE))
                            .map(metric -> {
                                final ClusteredCrosscheckMetric clusteredCrosscheckMetric = new ClusteredCrosscheckMetric(metric);
                                clusteredCrosscheckMetric.CLUSTER = cluster.getKey();
                                return clusteredCrosscheckMetric;
                            })
                            .collect(Collectors.toSet()));
        }

        return clusteredMetrics;
    }

}
