package picard.vcf;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public abstract class AbstractVcfMergingClpTester {

	protected static final String TEST_DATA_PATH = "testdata/picard/vcf/";

    protected abstract CommandLineProgram getProgram();

    protected void runClp(final List<File> inputs, final File output, final int expectedReturnCode) {
        runClp(inputs, output, Collections.<String>emptyList(), expectedReturnCode);
    }

    protected void runClp(final List<File> inputs, final File output, final List<String> otherArguments, final int expectedReturnCode) {
        final ArrayList<String> arguments = new ArrayList<String>();
        for (final File input : inputs) {
            arguments.add("INPUT=" + input);
        }
        arguments.add("OUTPUT=" + output);
        for (final String argument : otherArguments) {
            arguments.add(argument);
        }
        Assert.assertEquals(getProgram().instanceMain(arguments.toArray(new String[arguments.size()])), expectedReturnCode);
    }

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnDissimilarContigLists() {
		final File dissimilarContigs = new File(TEST_DATA_PATH, "CEUTrio-indels-dissimilar-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = new File("/dev/null/blah");
        final List<String> indexing = Arrays.asList("CREATE_INDEX=false");

        runClp(Arrays.asList(dissimilarContigs, snpInputFile), output, indexing, 0);
	}

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnNoContigList() {
		final File contiglessIndelFile = new File(TEST_DATA_PATH, "CEUTrio-indels-no-contigs.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = new File("/dev/null");

        runClp(Arrays.asList(contiglessIndelFile, snpInputFile), output, 1);
	}

	@Test (expectedExceptions = IllegalArgumentException.class)
	public void testFailsOnDissimilarSampleLists() {
		final File badSampleIndelFile = new File(TEST_DATA_PATH, "CEUTrio-indels-bad-samples.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
        final File output = new File("/dev/null");

        runClp(Arrays.asList(badSampleIndelFile, snpInputFile), output, 1);
	}

	@Test
	public void testMergeIndelsSnps() throws IOException {
		final File indelInputFile = new File(TEST_DATA_PATH, "CEUTrio-indels.vcf");
		final File snpInputFile = new File(TEST_DATA_PATH, "CEUTrio-snps.vcf");
		final File output = File.createTempFile("merge-indels-snps-test-output.", ".vcf");
        final List<String> indexing = Arrays.asList("CREATE_INDEX=false");
        output.deleteOnExit();

		final Queue<String> indelContigPositions = loadContigPositions(indelInputFile);
		final Queue<String> snpContigPositions = loadContigPositions(snpInputFile);

        runClp(Arrays.asList(indelInputFile, snpInputFile), output, indexing, 0);
        validateSnpAndIndelResults(output, indelContigPositions, snpContigPositions);
    }

    /**
     * Make sure that the order of the output file is identical to the order
     * of the input files by iterating through the output, making sure that,
     * if the context is an indel (snp), the next genomic position in the indel
     * (snp) queue is the same. Also make sure that the context is in the order
     * specified by the input files.
     */
    private void validateSnpAndIndelResults(final File output, final Queue<String> indelContigPositions, final Queue<String> snpContigPositions) {
        final VCFFileReader outputReader = new VCFFileReader(output, false);
        final VariantContextComparator outputComparator = outputReader.getFileHeader().getVCFRecordComparator();
        VariantContext last = null;
        final CloseableIterator<VariantContext> iterator = outputReader.iterator();
        while (iterator.hasNext()) {
            final VariantContext outputContext = iterator.next();
            if (outputContext.isIndel()) Assert.assertEquals(getContigPosition(outputContext), indelContigPositions.poll());
            if (outputContext.isSNP()) Assert.assertEquals(getContigPosition(outputContext), snpContigPositions.poll());
            if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) <= 0);
            last = outputContext;
        }
        iterator.close();

        // We should have polled everything off the indel (snp) queues
        Assert.assertEquals(indelContigPositions.size(), 0);
        Assert.assertEquals(snpContigPositions.size(), 0);
    }

    @Test
	public void testMergeRandomScatter() throws IOException {
		final File zero = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-0.vcf");
		final File one = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-1.vcf");
		final File two = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-2.vcf");
		final File three = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-3.vcf");
		final File four = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-4.vcf");
		final File five = new File(TEST_DATA_PATH, "CEUTrio-random-scatter-5.vcf");
        final List<File> inputs = Arrays.asList(zero, one, two, three, four, five);

		final List<Queue<String>> positionQueues = new ArrayList<Queue<String>>(6);
		positionQueues.add(0, loadContigPositions(zero));
		positionQueues.add(1, loadContigPositions(one));
		positionQueues.add(2, loadContigPositions(two));
		positionQueues.add(3, loadContigPositions(three));
		positionQueues.add(4, loadContigPositions(four));
		positionQueues.add(5, loadContigPositions(five));

        final List<String> indexing = Arrays.asList("CREATE_INDEX=false");
        final File output = File.createTempFile("random-scatter-test-output.", ".vcf");
		output.deleteOnExit();

        runClp(inputs, output, indexing, 0);
        validateResultsForMultipleInputs(output, positionQueues);
	}

    private void validateResultsForMultipleInputs(final File output, final List<Queue<String>> positionQueues) {
        final VCFFileReader outputReader = new VCFFileReader(output, false);
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
            if (last != null) Assert.assertTrue(outputComparator.compare(last, outputContext) <= 0);
            last = outputContext;
        }
        iterator.close();

        for (final Queue<String> positionQueue : positionQueues) {
            Assert.assertEquals(positionQueue.size(), 0);
        }
    }

    static Queue<String> loadContigPositions(final File inputFile) {
		final VCFFileReader reader = new VCFFileReader(inputFile, false);
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
