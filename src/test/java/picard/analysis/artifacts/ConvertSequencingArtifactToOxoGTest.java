package picard.analysis.artifacts;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.analysis.CollectOxoGMetrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Created by farjoun on 10/29/17.
 */
public class ConvertSequencingArtifactToOxoGTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");
    private static final File SAM_FILE = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
    private static final File REFERENCE_SEQUENCE = new File(TEST_DATA_DIR, "merger.fasta");

    @Test
    public void testEquivalence() throws IOException, NoSuchFieldException {
        final File input = SAM_FILE;

        final File outputFileOxoG = File.createTempFile("test", ".oxo_g_metrics", TEST_DATA_DIR);
        outputFileOxoG.deleteOnExit();

        final File outputFileArtifacts = File.createTempFile("test", "", TEST_DATA_DIR);
        outputFileArtifacts.deleteOnExit();
        new File(outputFileArtifacts + SequencingArtifactMetrics.PRE_ADAPTER_SUMMARY_EXT).deleteOnExit();
        new File(outputFileArtifacts + SequencingArtifactMetrics.PRE_ADAPTER_DETAILS_EXT).deleteOnExit();
        new File(outputFileArtifacts + SequencingArtifactMetrics.BAIT_BIAS_SUMMARY_EXT).deleteOnExit();
        new File(outputFileArtifacts + SequencingArtifactMetrics.BAIT_BIAS_DETAILS_EXT).deleteOnExit();
        new File(outputFileArtifacts + SequencingArtifactMetrics.ERROR_SUMMARY_EXT).deleteOnExit();

        final File convertedArtifacts = File.createTempFile("test", "", TEST_DATA_DIR);
        convertedArtifacts.deleteOnExit();
        new File (convertedArtifacts+".oxog_metrics").deleteOnExit();


        final List<String> args = new ArrayList<>();
        args.add("OUTPUT=" + outputFileOxoG.getAbsolutePath());
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath());

        CollectOxoGMetrics collectOxoGMetrics = new CollectOxoGMetrics();
        Assert.assertEquals(collectOxoGMetrics.instanceMain(args.toArray(new String[args.size()])), 0,
                "CollectOxoGMetrics can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        args.clear();
        args.add("OUTPUT=" + outputFileArtifacts.getAbsolutePath());
        args.add("INPUT=" + input.getAbsolutePath());
        args.add("REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath());

        CollectSequencingArtifactMetrics collectArtifactMetrics = new CollectSequencingArtifactMetrics();
        Assert.assertEquals(collectArtifactMetrics.instanceMain(args.toArray(new String[args.size()])), 0,
                "CollectSequencingArtifactMetrics can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        args.clear();
        args.add("OUTPUT_BASE=" + convertedArtifacts.getAbsolutePath());
        args.add("INPUT_BASE=" + outputFileArtifacts.getAbsolutePath());

        ConvertSequencingArtifactToOxoG convertSequencingArtifactToOxoG = new ConvertSequencingArtifactToOxoG();
        Assert.assertEquals(convertSequencingArtifactToOxoG.instanceMain(args.toArray(new String[args.size()])), 0,
                "ConvertSequencingArtifactToOxoG can't process base input" + outputFileArtifacts.getAbsolutePath() + " correctly");

        String[] columnToCompare=new String[]{
                "SAMPLE_ALIAS", "LIBRARY", "CONTEXT", "TOTAL_BASES", "REF_NONOXO_BASES", "REF_OXO_BASES", "REF_TOTAL_BASES", "ALT_NONOXO_BASES", "ALT_OXO_BASES", "OXIDATION_ERROR_RATE", "OXIDATION_Q", "C_REF_REF_BASES", "G_REF_REF_BASES", "C_REF_ALT_BASES", "G_REF_ALT_BASES", "C_REF_OXO_ERROR_RATE", "C_REF_OXO_Q", "G_REF_OXO_ERROR_RATE", "G_REF_OXO_Q"};

        Assert.assertTrue(areMetricsEqual(outputFileOxoG, new File(convertedArtifacts + ".oxog_metrics"), columnToCompare), "Metrics differ");
    }

    /**
     * Compare the metrics in two files, ignoring headers and histograms.
     */
    public static boolean areMetricsEqual(final File file1, final File file2, final String[] columnsToCompare) throws NoSuchFieldException {
        try {
            final MetricsFile<MetricBase, ?> mf1 = new MetricsFile<>();
            final MetricsFile<MetricBase, ?> mf2 = new MetricsFile<>();
            mf1.read(new FileReader(file1));
            mf2.read(new FileReader(file2));

            return areMetricsEqual(mf1,mf2,columnsToCompare);
        } catch (FileNotFoundException e) {
            throw new SAMException(e.getMessage(), e);
        }
    }

    private static <T extends MetricBase> boolean areMetricsEqual(MetricsFile<T,?> metricFile1, MetricsFile<T,?> metricFile2, final String[] columnsToCompare) throws NoSuchFieldException {

        if (metricFile1.getMetrics().size()!=metricFile2.getMetrics().size())
            return false;
        if (metricFile1.getMetrics().isEmpty()) return true;

        T firstMetric1 = metricFile1.getMetrics().get(0);

        for( final String column : columnsToCompare) {
            if (!metricFile1.getMetricsColumnLabels().contains(column))
                return false;
            if (!metricFile2.getMetricsColumnLabels().contains(column))
                return false;

            final Field f = firstMetric1.getClass().getField(column);

            List<String> metric1Values = metricFile1.getMetrics().stream().map(m -> {
                try {
                    return f.get(m).toString();
                } catch (IllegalAccessException e) {
                    e.printStackTrace();
                }
                return null;
            }).collect(Collectors.toList());
            List<String> metric2Values = metricFile2.getMetrics().stream().map(m -> {
                try {
                    return f.get(m).toString();
                } catch (IllegalAccessException e) {
                    e.printStackTrace();
                }
                return null;
            }).collect(Collectors.toList());

            if (!metric1Values.equals(metric2Values))
                return false;
        }

        return true;
    }
}