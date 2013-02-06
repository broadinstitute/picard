package net.sf.picard.vcf;

import net.sf.samtools.util.CloseableIterator;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class VariantContextIteratorTest {

	private static final String TEST_DATA_PATH = "testdata/net/sf/picard/vcf/";

	@Test
	public void testNext() {
		final File input = new File(TEST_DATA_PATH + "CEUTrio-merged-indels-snps.vcf");
		final CloseableIterator<VariantContext> variantContextIterator = new VariantContextIterator(input);

		Assert.assertNotNull(variantContextIterator.next());
		Assert.assertTrue(variantContextIterator.hasNext());
		Assert.assertNotNull(variantContextIterator.next());
	}
}
