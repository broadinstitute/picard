package picard.vcf;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;

import java.io.File;
import java.util.Queue;

public class SplitVcfsTest extends CommandLineProgramTest {

	private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("SplitVcfsTest", null);
	private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");

    public String getCommandLineProgramName() {
        return SplitVcfs.class.getSimpleName();
    }

	@AfterClass
	public void teardown() {
		IOUtil.deleteDirectoryTree(OUTPUT_DATA_PATH);
	}

	@Test
	public void testSplit() {
		final File indelOutputFile = new File(OUTPUT_DATA_PATH, "split-vcfs-test-indels-delete-me.vcf");
		final File snpOutputFile = new File(OUTPUT_DATA_PATH, "split-vcfs-test-snps-delete-me.vcf");
		final File input = new File(TEST_DATA_PATH, "CEUTrio-merged-indels-snps.vcf");

		indelOutputFile.deleteOnExit();
		snpOutputFile.deleteOnExit();

        final String[] args = new String[]{
                "INPUT=" + input.getAbsolutePath(),
                "SNP_OUTPUT=" + snpOutputFile.getAbsolutePath(),
                "INDEL_OUTPUT=" + indelOutputFile.getAbsolutePath()
        };
        Assert.assertEquals(runPicardCommandLine(args), 0);

		final Queue<String> indelContigPositions = AbstractVcfMergingClpTester.loadContigPositions(indelOutputFile);
		final Queue<String> snpContigPositions = AbstractVcfMergingClpTester.loadContigPositions(snpOutputFile);

		final VCFFileReader reader = new VCFFileReader(input);
        for (final VariantContext inputContext : reader) {
            if (inputContext.isIndel())
                Assert.assertEquals(AbstractVcfMergingClpTester.getContigPosition(inputContext), indelContigPositions.poll());
            if (inputContext.isSNP())
                Assert.assertEquals(AbstractVcfMergingClpTester.getContigPosition(inputContext), snpContigPositions.poll());
        }

		// We should have polled everything off the indel (snp) queues
		Assert.assertEquals(indelContigPositions.size(), 0);
		Assert.assertEquals(snpContigPositions.size(), 0);
	}
}
