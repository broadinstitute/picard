package picard.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.Test;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Queue;

public class SplitVcfsTest {

	private static final File OUTPUT_DATA_PATH = IOUtil.createTempDir("SplitVcfsTest", null);
	private static final File TEST_DATA_PATH = new File("testdata/picard/vcf/");

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

		final SplitVcfs splitVcfs = new SplitVcfs();
		splitVcfs.SNP_OUTPUT = snpOutputFile;
		splitVcfs.INDEL_OUTPUT = indelOutputFile;
		splitVcfs.INPUT = input;

		final int returnCode = splitVcfs.instanceMain(new String[0]);
		Assert.assertEquals(returnCode, 0);

		final Queue<String> indelContigPositions = MergeVcfsTest.loadContigPositions(indelOutputFile);
		final Queue<String> snpContigPositions = MergeVcfsTest.loadContigPositions(snpOutputFile);

		final VCFFileReader reader = new VCFFileReader(input);
		final CloseableIterator<VariantContext> iterator = reader.iterator();
		while (iterator.hasNext()) {
			final VariantContext inputContext = iterator.next();
			if (inputContext.isIndel()) Assert.assertEquals(MergeVcfsTest.getContigPosition(inputContext), indelContigPositions.poll());
			if (inputContext.isSNP()) Assert.assertEquals(MergeVcfsTest.getContigPosition(inputContext), snpContigPositions.poll());
		}

		// We should have polled everything off the indel (snp) queues
		Assert.assertEquals(indelContigPositions.size(), 0);
		Assert.assertEquals(snpContigPositions.size(), 0);
	}

	@Test (enabled = false)
	public void sampleVCF() {

		final int SAMPLE_FREQUENCY = 20000;

		int totalIn = 0;
		int totalOut = 0;
		final Map<Type, Integer> inputCounts = new HashMap<Type, Integer>();
		final Map<Type, Integer> outputCounts = new HashMap<Type, Integer>();
		final File INPUT = new File("/Volumes/Disko Segundo/splitvcfs/CEUTrio.HiSeq.WGS.b37.snps_and_indels.recalibrated.filtered.phased.CURRENT.vcf.gz");
		final VCFFileReader reader = new VCFFileReader(INPUT);

		final VariantContextWriter OUTPUT = new VariantContextWriterBuilder()
                .setOutputFile("/Volumes/shm/CEUTrio-REDUCED.vcf")
                .clearOptions()
                .build();
		OUTPUT.writeHeader(reader.getFileHeader());

		final CloseableIterator<VariantContext> iterator = reader.iterator();
		while (iterator.hasNext()) {
			final VariantContext variantContext = iterator.next();
			totalIn++;

			final Integer inputCount = inputCounts.get(variantContext.getType());
			if (inputCount == null) inputCounts.put(variantContext.getType(), 1);
			else inputCounts.put(variantContext.getType(), inputCount + 1);

			if ((totalIn % SAMPLE_FREQUENCY) == 0) {
				OUTPUT.add(variantContext);
				totalOut++;
				final Integer outputCount = outputCounts.get(variantContext.getType());
				if (outputCount == null) outputCounts.put(variantContext.getType(), 1);
				else outputCounts.put(variantContext.getType(), outputCount + 1);
			}
		}

		reader.close();
		OUTPUT.close();

		System.out.println("INPUT: " + totalIn + "; OUTPUT: " + totalOut);
		System.out.println("INPUT: " + inputCounts.toString());
		System.out.println("OUTPUT: " + outputCounts.toString());
	}
}
