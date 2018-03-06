package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by farjoun on 12/6/17.
 */
public class ClusterCrosscheckMetricsTest extends CommandLineProgramTest {

    @Override
    public String getCommandLineProgramName() {
        return ClusterCrosscheckMetrics.class.getSimpleName();
    }

    private final File TEST_DIR = new File("testdata/picard/fingerprint/");
    private final File HAPLOTYPE_MAP = new File(TEST_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");

    private final File NA12891_r1_sam = new File(TEST_DIR, "NA12891.over.fingerprints.r1.sam");
    private final File NA12891_r2_sam = new File(TEST_DIR, "NA12891.over.fingerprints.r2.sam");
    private final File NA12892_r1_sam = new File(TEST_DIR, "NA12892.over.fingerprints.r1.sam");
    private final File NA12892_r2_sam = new File(TEST_DIR, "NA12892.over.fingerprints.r2.sam");


    @Test
    public void testSimpleCluster() throws IOException {

        File metrics = File.createTempFile("Fingerprinting", "NA1291_and_NA12892.RG.crosscheck_metrics");
        metrics.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + NA12891_r1_sam.getAbsolutePath(),
                "INPUT=" + NA12892_r1_sam.getAbsolutePath(),
                "INPUT=" + NA12891_r2_sam.getAbsolutePath(),
                "INPUT=" + NA12892_r2_sam.getAbsolutePath(),
                "OUTPUT=" + metrics.getAbsolutePath(),
                "HAPLOTYPE_MAP=" + HAPLOTYPE_MAP,
                "LOD_THRESHOLD=" + -2.0,
                "CROSSCHECK_BY=FILE",
                "EXPECT_ALL_GROUPS_TO_MATCH=false"
        };

        final CrosscheckFingerprints crossChecker = new CrosscheckFingerprints();
        Assert.assertEquals(crossChecker.instanceMain(args), 0);

        File clusteredMetrics = File.createTempFile("Fingerprinting", "NA1291_and_NA12892.RG.clustered_crosscheck_metrics");
        clusteredMetrics.deleteOnExit();

        final String[] clusterArgs = new String[]{
                "INPUT=" + metrics.getAbsolutePath(),
                "OUTPUT=" + clusteredMetrics.getAbsolutePath(),
                "LOD_THRESHOLD=1.5"
        };

        final ClusterCrosscheckMetrics crosscheckerClusterer = new ClusterCrosscheckMetrics();
        Assert.assertEquals(crosscheckerClusterer.instanceMain(clusterArgs), 0);

        final MetricsFile<ClusteredCrosscheckMetric, Comparable<?>> metricsOutput = new MetricsFile<>();
        metricsOutput.read(new FileReader(clusteredMetrics));

        Assert.assertEquals(metricsOutput.getMetrics().size(),8);
        metricsOutput.getMetrics().forEach(m->Assert.assertEquals(m.CLUSTER_SIZE, (Integer) 2));

        final Map<String, Integer> clusterMap = new HashMap<>();

        metricsOutput.getMetrics().forEach(m -> {
            if (!clusterMap.containsKey(m.LEFT_SAMPLE)) {
                clusterMap.put(m.LEFT_SAMPLE, m.CLUSTER);
            } else {
                Assert.assertEquals(clusterMap.get(m.LEFT_SAMPLE), m.CLUSTER);
            }
            if (!clusterMap.containsKey(m.RIGHT_SAMPLE)) {
                clusterMap.put(m.RIGHT_SAMPLE, m.CLUSTER);
            } else {
                Assert.assertEquals(clusterMap.get(m.RIGHT_SAMPLE), m.CLUSTER);
            }
        });
    }
}