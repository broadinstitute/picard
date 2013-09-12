package net.sf.picard.vcf;

import net.sf.samtools.util.CloseableIterator;
import org.broad.tribble.TribbleException;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextComparator;
import org.broadinstitute.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public class MergeVcfsTest {

	private static final String TEST_DATA_PATH = "testdata/net/sf/picard/vcf/";

	@Test(enabled = false)
	public void testLong() {
		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = new File("/Users/jrose/development/long-merge-test.vcf.gz");
		mergeVcfs.CREATE_INDEX = false;
		mergeVcfs.INPUT = Arrays.asList(
				new File("/Volumes/Disko Segundo/mergevcfs/t2d_genes_contam_test4_per_sample_plus_five.snps.recalibrated.vcf"),
				new File("/Volumes/Disko Segundo/mergevcfs/t2d_genes_contam_test4_per_sample_plus_five.indels.filtered.vcf"));

		final int returnCode = mergeVcfs.instanceMain(new String[0]);
		Assert.assertEquals(returnCode, 0);
	}

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnDissimilarContigLists() {
		final File dissimilarContigs = new File(TEST_DATA_PATH, "CEUTrio-indels-dissimilar-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = new File("/dev/null/blah");
		mergeVcfs.CREATE_INDEX = false;
		mergeVcfs.INPUT = Arrays.asList(dissimilarContigs, snpInputFile);

		mergeVcfs.instanceMain(new String[0]);
	}

	@Test (expectedExceptions = TribbleException.class)
	public void testFailsOnNoContigList() {
		final File contiglessIndelFile = new File(TEST_DATA_PATH + "CEUTrio-indels-no-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = new File("/dev/null/blah");
		mergeVcfs.INPUT = Arrays.asList(contiglessIndelFile, snpInputFile);

		mergeVcfs.instanceMain(new String[0]);
	}

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnDissimilarSampleLists() {
		final File badSampleIndelFile = new File(TEST_DATA_PATH + "CEUTrio-indels-bad-samples.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = new File("/dev/null/blah");
		mergeVcfs.INPUT = Arrays.asList(badSampleIndelFile, snpInputFile);

		mergeVcfs.instanceMain(new String[0]);
	}

	@Test
	public void testMergeIndelsSnps() {
		final File indelInputFile = new File(TEST_DATA_PATH + "CEUTrio-indels.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH + "CEUTrio-snps.vcf");
		final File output = new File(TEST_DATA_PATH + "merge-indels-snps-test-output-delete-me.vcf");

		final Queue<String> indelContigPositions = loadContigPositions(indelInputFile);
		final Queue<String> snpContigPositions = loadContigPositions(snpInputFile);

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = output;
		mergeVcfs.INPUT = Arrays.asList(indelInputFile, snpInputFile);

		final int returnCode = mergeVcfs.instanceMain(new String[0]);
		Assert.assertEquals(returnCode, 0);

		// Make sure that the order of the output file is identical to the order
		// of the input files by iterating through the output, making sure that,
		// if the context is an indel (snp), the next genomic position in the indel
		// (snp) queue is the same. Also make sure that the context is in the order
		// specified by the input files.
		final VCFFileReader outputReader = new VCFFileReader(output);
		final VariantContextComparator outputComparator = outputReader.getFileHeader().getVCFRecordComparator();
		VariantContext last = null;
		final CloseableIterator<VariantContext> iterator = outputReader.iterator();
		while (iterator.hasNext()) {
			final VariantContext outputContext = iterator.next();
			if (outputContext.isIndel()) Assert.assertEquals(getContigPosition(outputContext), indelContigPositions.poll());
			if (outputContext.isSNP()) Assert.assertEquals(getContigPosition(outputContext), snpContigPositions.poll());
			if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) < 0);
			last = outputContext;
		}

		// We should have polled everything off the indel (snp) queues
		Assert.assertEquals(indelContigPositions.size(), 0);
		Assert.assertEquals(snpContigPositions.size(), 0);

		output.deleteOnExit();
	}

	@Test
	public void testMergeRandomScatter() {
		final File zero = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-0.vcf");
		final File one = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-1.vcf");
		final File two = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-2.vcf");
		final File three = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-3.vcf");
		final File four = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-4.vcf");
		final File five = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-5.vcf");

		final List<Queue<String>> positionQueues = new ArrayList<Queue<String>>(6);
		positionQueues.add(0, loadContigPositions(zero));
		positionQueues.add(1, loadContigPositions(one));
		positionQueues.add(2, loadContigPositions(two));
		positionQueues.add(3, loadContigPositions(three));
		positionQueues.add(4, loadContigPositions(four));
		positionQueues.add(5, loadContigPositions(five));

		final File output = new File(TEST_DATA_PATH + "random-scatter-test-output-delete-me.vcf");
		output.deleteOnExit();

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = output;
		mergeVcfs.INPUT = Arrays.asList(zero, one, two, three, four, five);

		final int returnCode = mergeVcfs.instanceMain(new String[0]);
		Assert.assertEquals(returnCode, 0);

		final VCFFileReader outputReader = new VCFFileReader(output);
		final VariantContextComparator outputComparator = outputReader.getFileHeader().getVCFRecordComparator();
		VariantContext last = null;
		final CloseableIterator<VariantContext> iterator = outputReader.iterator();
		while (iterator.hasNext()) {
			final VariantContext outputContext = iterator.next();
			final String position = getContigPosition(outputContext);
			for (final Queue<String> positionQueue : positionQueues) {
				if (position.equals(positionQueue.peek())) {
					positionQueue.poll();
					break;
				}
			}

			if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) < 0);
			last = outputContext;
		}

		for (final Queue<String> positionQueue : positionQueues) {
			Assert.assertEquals(positionQueue.size(), 0);
		}
	}

	static Queue<String> loadContigPositions(final File inputFile) {
		final VCFFileReader reader = new VCFFileReader(inputFile);
		final Queue<String> contigPositions = new LinkedList<String>();
		final CloseableIterator<VariantContext> iterator = reader.iterator();
		while (iterator.hasNext()) contigPositions.add(getContigPosition(iterator.next()));
		iterator.close();
		reader.close();
		return contigPositions;
	}

	static String getContigPosition(final VariantContext context) {
		return context.getChr() + "-" + Integer.toString(context.getStart());
	}
}
