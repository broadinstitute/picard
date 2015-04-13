package picard.analysis;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Tests the two default "programs" that have tests in CollectMultipleMetrics
 *
 *
 * @author Yossi farjoun
 */

public class CollectMultipleMetricsTest extends CommandLineProgramTest {

    private static final File TEST_DATA_DIR = new File("testdata/picard/sam");

    public String getCommandLineProgramName() {
        return CollectMultipleMetrics.class.getSimpleName();
    }


    @Test
    public void testAlignmentSummaryViaMultipleMetrics() throws IOException {
        final File input = new File(TEST_DATA_DIR, "summary_alignment_stats_test.sam");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File outfile   = File.createTempFile("alignmentMetrics", ".txt");
        outfile.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "PROGRAM=null",
                "PROGRAM="+CollectMultipleMetrics.Program.CollectAlignmentSummaryMetrics.name(),
                "PROGRAM="+CollectMultipleMetrics.Program.CollectInsertSizeMetrics.name()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<AlignmentSummaryMetrics, Comparable<?>> output = new MetricsFile<AlignmentSummaryMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final AlignmentSummaryMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.MEAN_READ_LENGTH, 101.0);
            switch (metrics.CATEGORY) {
                case FIRST_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 7);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 3);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 59);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 19.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 303);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*58D/303D*/0.191419);
                    Assert.assertEquals(metrics.BAD_CYCLES, 19);
                    break;
                case SECOND_OF_PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 9);
                    Assert.assertEquals(metrics.PF_READS, 9);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 1);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 7);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 239);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 707);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*19D/707D*/0.026874);
                    Assert.assertEquals(metrics.BAD_CYCLES, 3);
                    break;
                case PAIR:
                    Assert.assertEquals(metrics.TOTAL_READS, 18);
                    Assert.assertEquals(metrics.PF_READS, 16);
                    Assert.assertEquals(metrics.PF_NOISE_READS, 2);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_READS, 10);
                    Assert.assertEquals(metrics.PF_HQ_ALIGNED_Q20_BASES, 298);
                    Assert.assertEquals(metrics.PF_HQ_MEDIAN_MISMATCHES, 3.0);
                    Assert.assertEquals(metrics.PF_ALIGNED_BASES, 1010);
                    Assert.assertEquals(metrics.PF_MISMATCH_RATE, /*77D/1010D*/0.076238);
                    Assert.assertEquals(metrics.BAD_CYCLES, 22);
                    break;
                case UNPAIRED:
                default:
                    Assert.fail("Data does not contain this category: " + metrics.CATEGORY);
            }
        }
    }

    @Test
    public void testInsertSize() throws IOException {
        final File input = new File(TEST_DATA_DIR, "insert_size_metrics_test.sam");
        final File outfile   = File.createTempFile("test", ".insert_size_metrics");
        final File reference = new File(TEST_DATA_DIR, "summary_alignment_stats_test.fasta");
        final File pdf   = File.createTempFile("test", ".pdf");
        outfile.deleteOnExit();
        pdf.deleteOnExit();
        final String[] args = new String[] {
                "INPUT="  + input.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + reference.getAbsolutePath(),
                "PROGRAM=null",
                "PROGRAM="+CollectMultipleMetrics.Program.CollectAlignmentSummaryMetrics.name(),
                "PROGRAM="+CollectMultipleMetrics.Program.CollectInsertSizeMetrics.name()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        output.read(new FileReader(outfile));

        for (final InsertSizeMetrics metrics : output.getMetrics()) {
            Assert.assertEquals(metrics.PAIR_ORIENTATION.name(), "FR");
            if (metrics.LIBRARY==null) {  // SAMPLE or ALL_READS level
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 13);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            }
            else if (metrics.LIBRARY.equals("Solexa-41753")) { // one LIBRARY and one READ_GROUP
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 44);
                Assert.assertEquals(metrics.READ_PAIRS, 2);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);

            }
            else if (metrics.LIBRARY.equals("Solexa-41748") && metrics.READ_GROUP == null) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 40);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 9);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);

            }
            else if (metrics.LIBRARY.equals("Solexa-41734") && metrics.READ_GROUP == null) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 26);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 9);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.7")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 4);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 5);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 11);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 11);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.6")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 38);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 45);
                Assert.assertEquals(metrics.READ_PAIRS, 5);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 3);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 7);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 9);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 9);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.5")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 41);
                Assert.assertEquals(metrics.READ_PAIRS, 1);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);
            }
            else if (metrics.READ_GROUP.equals("62A79AAXX100907.3")) {
                Assert.assertEquals((int)metrics.MEDIAN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MIN_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.MAX_INSERT_SIZE, 36);
                Assert.assertEquals(metrics.READ_PAIRS, 1);
                Assert.assertEquals(metrics.WIDTH_OF_10_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_20_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_30_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_40_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_50_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_60_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_70_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_80_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_90_PERCENT, 1);
                Assert.assertEquals(metrics.WIDTH_OF_99_PERCENT, 1);
            }
            else {
                Assert.fail("Unexpected metric: " + metrics);
            }
        }
    }
}