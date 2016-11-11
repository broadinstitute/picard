package picard.fingerprint;

import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.analysis.CollectInsertSizeMetrics;
import picard.analysis.InsertSizeMetrics;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;


/**
 * Created by hogstrom on 10/31/16.
 */

public class CrosscheckReadGroupFingerprintsTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/fingerprint/");

    public String getCommandLineProgramName() {
        return CrosscheckReadGroupFingerprints.class.getSimpleName();
    }

    /**
     * verify that output file is genereated
     */
    @Test
    public void testCrossCheckOutput() throws IOException {
        final File input = new File("testdata/picard/sam/insert_size_metrics_test.sam");
        final File input1 = new File(TEST_DATA_DIR, "1026_d14.1_kneaddata_paired_hg19_bowtie2_contam.bam");
        final File input2 = new File(TEST_DATA_DIR, "1006_d2.1_kneaddata_paired_hg19_bowtie2_contam.bam");
        final File input3 = new File(TEST_DATA_DIR, "1006_d14.1_kneaddata_paired_hg19_bowtie2_contam.bam");

        final File input4 = new File(TEST_DATA_DIR, "1026_d14.1_kneaddata_paired_hg19_bowtie2_contam_PUmod.bam");
        final File input5 = new File(TEST_DATA_DIR, "1006_d14.1_kneaddata_paired_hg19_bowtie2_contam_PUmod.bam");
        final File input6 = new File(TEST_DATA_DIR, "1002_d14.1_kneaddata_paired_hg19_bowtie2_contam_PUmod.bam");

        final File SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING = new File(TEST_DATA_DIR, "Homo_sapiens_assembly19.haplotype_database.subset.txt");
        final File outfile = File.createTempFile("test", ".CrossCheckOut");
        //outfile.deleteOnExit();
        System.out.print(outfile.getAbsolutePath());
        final String[] args = new String[]{
                "I=" + input5.getAbsolutePath(),
                "I=" + input6.getAbsolutePath(),
                "H=" + SUBSETTED_HAPLOTYPE_DATABASE_FOR_TESTING.getAbsolutePath(),
                "OUTPUT=" + outfile.getAbsolutePath()
                //"EXPECT_ALL_READ_GROUPS_TO_MATCH=true"
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

        //final MetricsFile<InsertSizeMetrics, Comparable<?>> output = new MetricsFile<InsertSizeMetrics, Comparable<?>>();
        //output.read(new FileReader(outfile));
        //Assert.assertEquals(output.getAllHistograms().size(), 5);
    }

}