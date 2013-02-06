package net.sf.picard.vcf;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.broadinstitute.variant.vcf.VCFContigHeaderLine;
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

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnNoContigList() {
		final File contiglessIndelFile = new File(TEST_DATA_PATH + "CEUTrio-indels-no-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = new File("/dev/null/blah");
		mergeVcfs.CREATE_INDEX = false;
		mergeVcfs.INPUT = Arrays.asList(contiglessIndelFile, snpInputFile);

		mergeVcfs.instanceMain(new String[0]);
	}

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnDissimilarSampleLists() {
		final File badSampleIndelFile = new File(TEST_DATA_PATH + "CEUTrio-indels-bad-samples.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = new File("/dev/null/blah");
		mergeVcfs.CREATE_INDEX = false;
		mergeVcfs.INPUT = Arrays.asList(badSampleIndelFile, snpInputFile);

		mergeVcfs.instanceMain(new String[0]);
	}

	@Test
	public void testMergeIndelsSnps() {
		final File indelInputFile = new File(TEST_DATA_PATH + "CEUTrio-indels.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH + "CEUTrio-snps.vcf");
		final File output = new File(TEST_DATA_PATH + "merge-indels-snps-test-output-delete-me.vcf");

		final VariantContextIterator indelIterator = new VariantContextIterator(indelInputFile);
		final VariantContextIterator snpIterator = new VariantContextIterator(snpInputFile);

		final Queue<String> indelContigPositions = getContigPositions(indelIterator);
		final Queue<String> snpContigPositions = getContigPositions(snpIterator);

		indelIterator.close();
		snpIterator.close();

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = output;
		mergeVcfs.CREATE_INDEX = false;
		mergeVcfs.INPUT = Arrays.asList(indelInputFile, snpInputFile);

		final int returnCode = mergeVcfs.instanceMain(new String[0]);
		Assert.assertEquals(returnCode, 0);

		// Make sure that the order of the output file is identical to the order
		// of the input files by iterating through the output, making sure that,
		// if the context is an indel (snp), the next genomic position in the indel
		// (snp) queue is the same. Also make sure that the context is in the order
		// specified by the input files.
		final VariantContextIterator outputIterator = new VariantContextIterator(output);
		final VariantContextComparator outputComparator = new VariantContextComparator(getContigs(outputIterator));
		VariantContext last = null;
		while (outputIterator.hasNext()) {
			final VariantContext outputContext = outputIterator.next();
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
		positionQueues.add(0, getContigPositions(new VariantContextIterator(zero)));
		positionQueues.add(1, getContigPositions(new VariantContextIterator(one)));
		positionQueues.add(2, getContigPositions(new VariantContextIterator(two)));
		positionQueues.add(3, getContigPositions(new VariantContextIterator(three)));
		positionQueues.add(4, getContigPositions(new VariantContextIterator(four)));
		positionQueues.add(5, getContigPositions(new VariantContextIterator(five)));

		final File output = new File(TEST_DATA_PATH + "random-scatter-test-output-delete-me.vcf");

		final MergeVcfs mergeVcfs = new MergeVcfs();
		mergeVcfs.OUTPUT = output;
		mergeVcfs.CREATE_INDEX = false;
		mergeVcfs.INPUT = Arrays.asList(zero, one, two, three, four, five);

		final int returnCode = mergeVcfs.instanceMain(new String[0]);
		Assert.assertEquals(returnCode, 0);

		final VariantContextIterator outputIterator = new VariantContextIterator(output);
		final VariantContextComparator outputComparator = new VariantContextComparator(getContigs(outputIterator));
		VariantContext last = null;
		while (outputIterator.hasNext()) {
			final VariantContext outputContext = outputIterator.next();
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

		output.deleteOnExit();
	}

	@Test (enabled = false)
	public void dumpHeaders() {
		final File[] files = new File[] {
				new File("/Volumes/Disko Segundo/mergevcfs/t2d_genes_contam_test4_per_sample_plus_five.snps.recalibrated.vcf"),
				new File("/Volumes/Disko Segundo/mergevcfs/t2d_genes_contam_test4_per_sample_plus_five.indels.filtered.vcf"),
				new File("/Volumes/Disko Segundo/mergevcfs/t2d_genes_contam_test4_per_sample_plus_five.unannotated.vcf"),
				new File("/Users/jrose/development/long-merge-test.vcf")
		};
		for (final File file : files) {
			final VariantContextIterator iterator = new VariantContextIterator(file);
			final File output = new File("/Volumes/Disko Segundo/mergevcfs/", file.getName() + ".header");
			final VariantContextWriter writer = VariantContextWriterFactory.create(output, null, VariantContextWriterFactory.NO_OPTIONS);
			writer.writeHeader(iterator.getHeader());
			writer.close();
		}
	}

	static Queue<String> getContigPositions(final VariantContextIterator iterator) {
		final Queue<String> contigPositions = new LinkedList<String>();
		while (iterator.hasNext()) contigPositions.add(getContigPosition(iterator.next()));
		return contigPositions;
	}

	static String getContigPosition(final VariantContext context) {
		return context.getChr() + "-" + Integer.toString(context.getStart());
	}

	static List<String> getContigs(final VariantContextIterator iterator) {
		final List<String> contigList = new ArrayList<String>();
		for (final VCFContigHeaderLine contigHeaderLine : iterator.getHeader().getContigLines()) {
			contigList.add(contigHeaderLine.getID());
		}

		return contigList;
	}
}
