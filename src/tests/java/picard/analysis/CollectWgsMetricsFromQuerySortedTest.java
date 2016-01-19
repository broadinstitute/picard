package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.List;

/**
 * Tests for methods in CollectWgsMetricsFromQuerySorted
 *
 *
 * @author Eric Banks
 */

public class CollectWgsMetricsFromQuerySortedTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    public String getCommandLineProgramName() {
        return CollectWgsMetricsFromQuerySorted.class.getSimpleName();
    }


    @Test
    public void testMetricsFromClippedOverhangs() throws IOException {
        final File input = new File(TEST_DATA_DIR, "namesorted.test.sam");
        final File outfile   = File.createTempFile("metrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));
        validateMetrics(output.getMetrics(), 3095693981L);
    }

    @Test
    public void testPassingInGenomeTerritory() throws IOException {
        final File input = new File(TEST_DATA_DIR, "namesorted.test.sam");
        final File outfile   = File.createTempFile("metrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "GENOME_TERRITORY=1000"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics, Comparable<?>> output = new MetricsFile<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));
        validateMetrics(output.getMetrics(), 1000L);
    }

    private void validateMetrics(final List<CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics> metrics, final long genomeSize) {
        for (final CollectWgsMetricsFromQuerySorted.QuerySortedSeqMetrics row : metrics) {
            final boolean isRaw = row.TYPE == CollectWgsMetricsFromQuerySorted.FILTERING_STRINGENCY.RAW;

            Assert.assertEquals(row.GENOME_TERRITORY, genomeSize);
            Assert.assertEquals(row.PF_BASES, 606);
            Assert.assertEquals(row.PF_PASSING_BASES, isRaw ? 238 : 200);
            Assert.assertEquals(row.PCT_EXC_OVERLAP, isRaw ? 0.085809 : 0.013201);  // raw: 52/606, usable: 8/606
            Assert.assertEquals(row.PCT_EXC_BASEQ, isRaw ? 0.188119 : 0.156766);    // raw: 114/606, usable 95/606
            Assert.assertEquals(row.PCT_EXC_MAPQ, isRaw ? 0.0 : 0.166667);          // raw: 0/606, usable:101/606
            Assert.assertEquals(row.PCT_EXC_DUPE, 0.333333);                        // both: 202/606
            Assert.assertEquals(row.PF_READ_PAIRS, 3);
            Assert.assertEquals(row.PF_DUPE_PAIRS, 1);
            Assert.assertEquals(row.PF_READS_ALIGNED, 6);
            Assert.assertEquals(row.PF_ORIENTED_PAIRS, 2);
            Assert.assertEquals(row.MEAN_INSERT_SIZE, 118.0);

            final BigDecimal meanCov = new BigDecimal((double)row.PF_PASSING_BASES / genomeSize).setScale(6, RoundingMode.HALF_UP);
            Assert.assertEquals(Double.compare(row.MEAN_COVERAGE, meanCov.doubleValue()), 0);
        }
    }
}