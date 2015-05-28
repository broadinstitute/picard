package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created by kbergin on 3/26/15 to test GcBias MultiLevel Collector.
 * Note: 'accLevel' needs to be changed to test various accumulation levels with this test class
 */
public class CollectGcBiasMetricsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/");

    public String getCommandLineProgramName() {
        return CollectGcBiasMetrics.class.getSimpleName();
    }

    @Test
    public void test() throws IOException {
        final File input = new File(TEST_DATA_DIR, "gc_bias_metrics_test.bam");
        final File Soutfile   = File.createTempFile("test", ".gc_bias_summary_metrics");
        final File Doutfile = File.createTempFile("test", ".gc_bias_detail_metrics");
        final File pdf   = File.createTempFile("test", ".pdf");
        final String referenceFile = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
        //final String referenceFile = "Homo_sapiens_assembly19.fasta";
        final String accLevel = "ALL_READS"; //Options: ALL_READS, LIBRARY, SAMPLE, READ_GROUP or a combination
        final int windowSize = 100;
        final double minGenFraction = 1.0E-5;
        final boolean biSulfiteSeq = false;
        final boolean assumeSorted = false;
        Soutfile.deleteOnExit();
        Doutfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + Doutfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + referenceFile,
                "SUMMARY_OUTPUT=" + Soutfile.getAbsolutePath(),
                "CHART_OUTPUT=" + pdf.getAbsolutePath(),
                "WINDOW_SIZE=" + windowSize,
                "MINIMUM_GENOME_FRACTION=" + minGenFraction,
                "IS_BISULFITE_SEQUENCED=" + biSulfiteSeq,
                "METRIC_ACCUMULATION_LEVEL=" + accLevel,
                "ASSUME_SORTED=" + assumeSorted
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<GcBiasSummaryMetrics, Comparable<?>> output = new MetricsFile<GcBiasSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(Soutfile));

        for (final GcBiasSummaryMetrics metrics : output.getMetrics()) {
            if (metrics.ACCUMULATION_LEVEL.equals("All Reads")) { //ALL_READS level
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 26);
                Assert.assertEquals(metrics.ALIGNED_READS, 39);
                Assert.assertEquals(metrics.AT_DROPOUT, 76.137572);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.LIBRARY != null && metrics.LIBRARY.equals("Solexa-41753")) { //One Library over one read group
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 9);
                Assert.assertEquals(metrics.ALIGNED_READS, 11);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.625488);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.LIBRARY!=null && metrics.LIBRARY.equals("Solexa-41748")) { //RG 6 and 7
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 14);
                Assert.assertEquals(metrics.ALIGNED_READS, 23);
                Assert.assertEquals(metrics.AT_DROPOUT, 76.137572);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.LIBRARY!=null && metrics.LIBRARY.equals("Solexa-41734")) { // RG 3 and 5
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 3);
                Assert.assertEquals(metrics.ALIGNED_READS, 5);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.625488);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.READ_GROUP!=null && metrics.READ_GROUP.equals("62A79AAXX100907.7")) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 5);
                Assert.assertEquals(metrics.ALIGNED_READS, 9);
                Assert.assertEquals(metrics.AT_DROPOUT, 76.137572);
                Assert.assertEquals(metrics.GC_DROPOUT, 19.03004);
            } else if (metrics.READ_GROUP!=null && metrics.READ_GROUP.equals("62A79AAXX100907.6")) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 9);
                Assert.assertEquals(metrics.ALIGNED_READS, 14);
                Assert.assertEquals(metrics.AT_DROPOUT, 76.137572);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.READ_GROUP!=null && metrics.READ_GROUP.equals("62A79AAXX100907.5")) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 2);
                Assert.assertEquals(metrics.ALIGNED_READS, 3);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.625488);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.READ_GROUP!=null && metrics.READ_GROUP.equals("62A79AAXX100907.3")) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 1);
                Assert.assertEquals(metrics.ALIGNED_READS, 2);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.625488);
                Assert.assertEquals(metrics.GC_DROPOUT, 19.03004);
            } else if (metrics.READ_GROUP!=null && metrics.READ_GROUP.equals("62A79AAXX100907.8")) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 9);
                Assert.assertEquals(metrics.ALIGNED_READS, 11);
                Assert.assertEquals(metrics.AT_DROPOUT, 78.625488);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else if (metrics.SAMPLE.equals("NA12878")) {
                Assert.assertEquals(metrics.TOTAL_CLUSTERS, 26);
                Assert.assertEquals(metrics.ALIGNED_READS, 39);
                Assert.assertEquals(metrics.AT_DROPOUT, 76.137572);
                Assert.assertEquals(metrics.GC_DROPOUT, 16.853326);
            } else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }
}

