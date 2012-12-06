package net.sf.picard.util;

import net.sf.picard.fastq.FastqReader;
import net.sf.samtools.SAMFileReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class QualityEncodingDetectorTest {

    private static class Testcase {
        private final File f;
        private final FastqQualityFormat q;

        Testcase(final File file, final FastqQualityFormat qualityFormat) {
            this.f = file;
            this.q = qualityFormat;
        }
    }

    final static List<Testcase> FASTQ_TESTCASES = Arrays.asList(
            // Need to use full-range quality here, as Solexa and Illumina are near indistinguishable
            new Testcase(new File("./testdata/net/sf/picard/sam/fastq2bam/fastq-solexa/solexa_full_range_as_solexa.fastq"), FastqQualityFormat.Solexa),
            new Testcase(new File("./testdata/net/sf/picard/sam/fastq2bam/fastq-illumina/s_1_sequence.txt"), FastqQualityFormat.Illumina),
            new Testcase(new File("./testdata/net/sf/picard/sam/fastq2bam/fastq-sanger/5k-30BB2AAXX.3.aligned.sam.fastq"), FastqQualityFormat.Standard)
    );
    final static List<Testcase> BAM_TESTCASES = Arrays.asList(
            new Testcase(new File("./testdata/net/sf/picard/sam/unmapped.sam"), FastqQualityFormat.Standard),
            new Testcase(new File("./testdata/net/sf/samtools/BAMFileIndexTest/index_test.bam"), FastqQualityFormat.Standard),
            new Testcase(new File("./testdata/net/sf/picard/util/QualityEncodingDetectorTest/solexa-as-standard.bam"), FastqQualityFormat.Solexa),
            new Testcase(new File("./testdata/net/sf/picard/util/QualityEncodingDetectorTest/illumina-as-standard.bam"), FastqQualityFormat.Illumina)
            
    );

    Object[][] renderObjectArrayArray(final List<Testcase> testcaseList) {
        final Object[][] data = new Object[testcaseList.size()][];
        for (int i = 0; i < data.length; i++) {
            final Testcase testcase = testcaseList.get(i);
            data[i] = new Object[]{testcase.f, testcase.q};
        }
        return data;
    }

    @DataProvider(name = "BAM_TESTCASES")
    Object[][] bamTestcases() {
        return renderObjectArrayArray(BAM_TESTCASES);
    }

    @DataProvider(name = "FASTQ_TESTCASES")
    Object[][] fastqTestcases() {
        return renderObjectArrayArray(FASTQ_TESTCASES);
    }

    @Test(dataProvider = "FASTQ_TESTCASES", groups = {"unix"})
    public void testFastqQualityInference(final File input, final FastqQualityFormat expectedQualityFormat) {
        final FastqReader reader = new FastqReader(input);
        Assert.assertEquals(QualityEncodingDetector.detect(reader), expectedQualityFormat);
        reader.close();
    }

    @Test(dataProvider = "BAM_TESTCASES", groups = {"unix"})
    public void testBamQualityInference(final File input, final FastqQualityFormat expectedQualityFormat) {
        final SAMFileReader reader = new SAMFileReader(input);
        Assert.assertEquals(QualityEncodingDetector.detect(reader), expectedQualityFormat);
        reader.close();
    }
}
