package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

public class CollectOxoGMetricsTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");
    private static final File SAM_FILE = new File(CollectAlignmentSummaryMetricsTest.TEST_DATA_DIR, "summary_alignment_stats_test.sam");
    private static final File REFERENCE_SEQUENCE = new File(TEST_DATA_DIR, "merger.fasta");

    @Test
    public void testCollectOxoGMetrics() throws IOException {
        final File outputFile = File.createTempFile("test", ".oxo_g_metrics", TEST_DATA_DIR);
        outputFile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + SAM_FILE.getAbsolutePath(),
                "OUTPUT=" + outputFile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath()
        };
        CollectOxoGMetrics collectOxoGMetrics = new CollectOxoGMetrics();
        Assert.assertEquals(collectOxoGMetrics.instanceMain(args), 0,
                "Can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        final MetricsFile<CollectOxoGMetrics.CpcgMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outputFile));

        {
            final int metricsNumber = 4; // metrics number for testing (randomly chosen) that corresponds "TCT" context string.
            final CollectOxoGMetrics.CpcgMetrics metrics = output.getMetrics().get(metricsNumber);

            Assert.assertEquals(metrics.SAMPLE_ALIAS, "Hi,Momma!");
            Assert.assertEquals(metrics.LIBRARY, "whatever");
            Assert.assertEquals(metrics.CONTEXT, "TCT");
            Assert.assertEquals(metrics.TOTAL_SITES, 3);
            Assert.assertEquals(metrics.TOTAL_BASES, 3);
            Assert.assertEquals(metrics.REF_NONOXO_BASES, 0);
            Assert.assertEquals(metrics.REF_OXO_BASES, 3);
            Assert.assertEquals(metrics.REF_TOTAL_BASES, 3);
            Assert.assertEquals(metrics.ALT_NONOXO_BASES, 0);
            Assert.assertEquals(metrics.ALT_OXO_BASES, 0);
            Assert.assertEquals(metrics.OXIDATION_ERROR_RATE, 0.333333);
            Assert.assertEquals(metrics.OXIDATION_Q, 4.771213);
            Assert.assertEquals(metrics.C_REF_REF_BASES, 0);
            Assert.assertEquals(metrics.G_REF_REF_BASES, 3);
            Assert.assertEquals(metrics.C_REF_ALT_BASES, 0);
            Assert.assertEquals(metrics.G_REF_ALT_BASES, 0);
            Assert.assertEquals(metrics.C_REF_OXO_ERROR_RATE, Double.NaN);
            Assert.assertEquals(metrics.C_REF_OXO_Q, Double.NaN);
            Assert.assertEquals(metrics.G_REF_OXO_ERROR_RATE, Double.NaN);
            Assert.assertEquals(metrics.G_REF_OXO_Q, Double.NaN);

        }

        {
            // metrics number for testing that corresponds "GCA" context string
            // (which is at the end of chr8 and is covered by a read).
            final int metricsNumber = 3;
            final CollectOxoGMetrics.CpcgMetrics metrics = output.getMetrics().get(metricsNumber);

            Assert.assertEquals(metrics.SAMPLE_ALIAS, "Hi,Momma!");
            Assert.assertEquals(metrics.LIBRARY, "whatever");
            Assert.assertEquals(metrics.CONTEXT, "GCA");
            Assert.assertEquals(metrics.TOTAL_SITES, 1);
            Assert.assertEquals(metrics.TOTAL_BASES, 1);
            Assert.assertEquals(metrics.REF_NONOXO_BASES, 1);
            Assert.assertEquals(metrics.REF_OXO_BASES, 0);
            Assert.assertEquals(metrics.REF_TOTAL_BASES, 1);
            Assert.assertEquals(metrics.ALT_NONOXO_BASES, 0);
            Assert.assertEquals(metrics.ALT_OXO_BASES, 0);
            Assert.assertEquals(metrics.OXIDATION_ERROR_RATE, 1D);
            Assert.assertEquals(metrics.OXIDATION_Q, -0D);
            Assert.assertEquals(metrics.C_REF_REF_BASES, 1);
            Assert.assertEquals(metrics.G_REF_REF_BASES, 0);
            Assert.assertEquals(metrics.C_REF_ALT_BASES, 0);
            Assert.assertEquals(metrics.G_REF_ALT_BASES, 0);
            Assert.assertEquals(metrics.C_REF_OXO_ERROR_RATE, Double.NaN);
            Assert.assertEquals(metrics.C_REF_OXO_Q, Double.NaN);
            Assert.assertEquals(metrics.G_REF_OXO_ERROR_RATE, Double.NaN);
            Assert.assertEquals(metrics.G_REF_OXO_Q, Double.NaN);
        }
    }


    @Test
    public void testCollectOxoGMetricsShortContext() throws IOException {
        final File outputFile = File.createTempFile("test", ".oxo_g_metrics", TEST_DATA_DIR);
        outputFile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + SAM_FILE.getAbsolutePath(),
                "OUTPUT=" + outputFile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath(),
                "CONTEXT_SIZE=0"
        };
        CollectOxoGMetrics collectOxoGMetrics = new CollectOxoGMetrics();
        Assert.assertEquals(collectOxoGMetrics.instanceMain(args), 0,
                "Can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        final MetricsFile<CollectOxoGMetrics.CpcgMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outputFile));

        final int metricsNumber = 0; // metrics number that corresponds "C" context string.
        final CollectOxoGMetrics.CpcgMetrics metrics = output.getMetrics().get(metricsNumber);

        Assert.assertEquals(metrics.SAMPLE_ALIAS, "Hi,Momma!");
        Assert.assertEquals(metrics.LIBRARY, "whatever");
        Assert.assertEquals(metrics.CONTEXT, "C");
        Assert.assertEquals(metrics.TOTAL_SITES, 22);
        Assert.assertEquals(metrics.TOTAL_BASES, 22);
        Assert.assertEquals(metrics.REF_NONOXO_BASES, 9);
        Assert.assertEquals(metrics.REF_OXO_BASES, 13);
        Assert.assertEquals(metrics.REF_TOTAL_BASES, 22);
        Assert.assertEquals(metrics.ALT_NONOXO_BASES, 0);
        Assert.assertEquals(metrics.ALT_OXO_BASES, 0);
        Assert.assertEquals(metrics.OXIDATION_ERROR_RATE, .045455);
        Assert.assertEquals(metrics.OXIDATION_Q, 13.424227);
        Assert.assertEquals(metrics.C_REF_REF_BASES, 9);
        Assert.assertEquals(metrics.G_REF_REF_BASES, 13);
        Assert.assertEquals(metrics.C_REF_ALT_BASES, 0);
        Assert.assertEquals(metrics.G_REF_ALT_BASES, 0);
        Assert.assertEquals(metrics.C_REF_OXO_ERROR_RATE, 0D);
        Assert.assertEquals(metrics.C_REF_OXO_Q, 100D);
        Assert.assertEquals(metrics.G_REF_OXO_ERROR_RATE, 0D);
        Assert.assertEquals(metrics.G_REF_OXO_Q, 100D);
    }

    @Test
    public void testCollectOxoGMetricsLongContext() throws IOException {
        final File outputFile = File.createTempFile("test", ".oxo_g_metrics", TEST_DATA_DIR);
        outputFile.deleteOnExit();
        final String[] args = new String[]{
                "INPUT=" + SAM_FILE.getAbsolutePath(),
                "OUTPUT=" + outputFile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + REFERENCE_SEQUENCE.getAbsolutePath(),
                "CONTEXT_SIZE=2"
        };
        CollectOxoGMetrics collectOxoGMetrics = new CollectOxoGMetrics();
        Assert.assertEquals(collectOxoGMetrics.instanceMain(args), 0,
                "Can't process " + SAM_FILE.getAbsolutePath() + " correctly");

        final MetricsFile<CollectOxoGMetrics.CpcgMetrics, Comparable<?>> output = new MetricsFile<>();
        output.read(new FileReader(outputFile));

        final int metricsNumber = 141; // metrics number that corresponds "GGCTG" context string.
        final CollectOxoGMetrics.CpcgMetrics metrics = output.getMetrics().get(metricsNumber);

        Assert.assertEquals(metrics.SAMPLE_ALIAS, "Hi,Momma!");
        Assert.assertEquals(metrics.LIBRARY, "whatever");
        Assert.assertEquals(metrics.CONTEXT, "GGCTG");
        Assert.assertEquals(metrics.TOTAL_SITES, 1);
        Assert.assertEquals(metrics.TOTAL_BASES, 1);
        Assert.assertEquals(metrics.REF_NONOXO_BASES, 0);
        Assert.assertEquals(metrics.REF_OXO_BASES, 1);
        Assert.assertEquals(metrics.REF_TOTAL_BASES, 1);
        Assert.assertEquals(metrics.ALT_NONOXO_BASES, 0);
        Assert.assertEquals(metrics.ALT_OXO_BASES, 0);
        Assert.assertEquals(metrics.OXIDATION_ERROR_RATE, 1D);
        Assert.assertEquals(metrics.OXIDATION_Q, -0D);
        Assert.assertEquals(metrics.C_REF_REF_BASES, 0);
        Assert.assertEquals(metrics.G_REF_REF_BASES, 1);
        Assert.assertEquals(metrics.C_REF_ALT_BASES, 0);
        Assert.assertEquals(metrics.G_REF_ALT_BASES, 0);
        Assert.assertEquals(metrics.C_REF_OXO_ERROR_RATE, Double.NaN);
        Assert.assertEquals(metrics.C_REF_OXO_Q, Double.NaN);
        Assert.assertEquals(metrics.G_REF_OXO_ERROR_RATE, Double.NaN);
        Assert.assertEquals(metrics.G_REF_OXO_Q, Double.NaN);
    }


    @DataProvider(name = "RightOptions")
    public static Object[][] rightOptions() {
        final HashSet<String> rightContext1 = new HashSet<>();
        rightContext1.add("ACC");

        final HashSet<String> rightContext2 = new HashSet<>();
        rightContext2.add("AACAA");
        rightContext2.add("ATCAT");
        return new Object[][] {
                {5, 10, 1, rightContext1}, //contextSize = 1
                {10, 10, 2, rightContext2} //contextSize = 2
        };
    }

    @Test(dataProvider = "RightOptions")
    public void testPositiveCustomCommandLineValidation(final int minimumInsertSize,
                                                        final int maximumInsertSize,
                                                        final int contextSize,
                                                        final HashSet<String> context) throws Exception {
        final CollectOxoGMetrics collectOxoGMetrics = getCollectOxoGMetrics(minimumInsertSize, maximumInsertSize, contextSize, context);
        Assert.assertNull(collectOxoGMetrics.customCommandLineValidation());
        Assert.assertEquals(collectOxoGMetrics.MINIMUM_INSERT_SIZE, minimumInsertSize);
        Assert.assertEquals(collectOxoGMetrics.MAXIMUM_INSERT_SIZE, maximumInsertSize);
        Assert.assertEquals(collectOxoGMetrics.CONTEXT_SIZE, contextSize);
        Assert.assertEquals(collectOxoGMetrics.CONTEXTS, context);
    }

    @DataProvider(name = "WrongOptions")
    public static Object[][] wrongOptions() {
        //Middle base of context sequence must be C
        final HashSet<String> wrongContext1 = new HashSet<>();
        wrongContext1.add("AAC");

        //Middle base of context sequence must be C
        final HashSet<String> wrongContext2 = new HashSet<>();
        wrongContext2.add("AAGAA");
        wrongContext2.add("ATCAT");

        final HashSet<String> rightContext1 = new HashSet<>();
        rightContext1.add("ACC");

        return new Object[][] {
                {5, 10, 1, wrongContext1, "Middle base of context sequence AAC must be C"},
                {10, 10, 1, wrongContext2, "Context AAGAA is not 3 long as implied by CONTEXT_SIZE=1"},
                {10, 5, 1, rightContext1, "MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE"}, //min insert size > max insert size
                {-5, 10, 1, rightContext1, "MINIMUM_INSERT_SIZE cannot be negative"} //negative insert size
        };
    }

    @Test(dataProvider = "WrongOptions")
    public void testNegativeCustomCommandLineValidation(final int minimumInsertSize,
                                                        final int maximumInsertSize,
                                                        final int contextSize,
                                                        final HashSet<String> context,
                                                        final String expectedMessage) throws Exception {
        final CollectOxoGMetrics collectOxoGMetrics = getCollectOxoGMetrics(minimumInsertSize, maximumInsertSize, contextSize, context);
        Assert.assertNotNull(collectOxoGMetrics.customCommandLineValidation());
        Assert.assertEquals(collectOxoGMetrics.customCommandLineValidation()[0], expectedMessage);
    }

    private static CollectOxoGMetrics getCollectOxoGMetrics(final int minimumInsertSize,
                                                     final int maximumInsertSize,
                                                     final int contextSize,
                                                     final HashSet<String> context) {
        final CollectOxoGMetrics collectOxoGMetrics = new CollectOxoGMetrics();
        collectOxoGMetrics.MINIMUM_INSERT_SIZE = minimumInsertSize;
        collectOxoGMetrics.MAXIMUM_INSERT_SIZE = maximumInsertSize;
        collectOxoGMetrics.CONTEXT_SIZE = contextSize;
        collectOxoGMetrics.CONTEXTS = context;
        return collectOxoGMetrics;
    }
}
