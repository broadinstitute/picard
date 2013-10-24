package org.broadinstitute.variant.vcf;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;

public class VCFEncoderTest {

	@DataProvider(name = "VCFWriterDoubleFormatTestData")
	public Object[][] makeVCFWriterDoubleFormatTestData() {
		List<Object[]> tests = new ArrayList<Object[]>();
		tests.add(new Object[]{1.0, "1.00"});
		tests.add(new Object[]{10.1, "10.10"});
		tests.add(new Object[]{10.01, "10.01"});
		tests.add(new Object[]{10.012, "10.01"});
		tests.add(new Object[]{10.015, "10.02"});
		tests.add(new Object[]{0.0, "0.00"});
		tests.add(new Object[]{0.5, "0.500"});
		tests.add(new Object[]{0.55, "0.550"});
		tests.add(new Object[]{0.555, "0.555"});
		tests.add(new Object[]{0.5555, "0.556"});
		tests.add(new Object[]{0.1, "0.100"});
		tests.add(new Object[]{0.050, "0.050"});
		tests.add(new Object[]{0.010, "0.010"});
		tests.add(new Object[]{0.012, "0.012"});
		tests.add(new Object[]{0.0012, "1.200e-03"});
		tests.add(new Object[]{1.2e-4, "1.200e-04"});
		tests.add(new Object[]{1.21e-4, "1.210e-04"});
		tests.add(new Object[]{1.212e-5, "1.212e-05"});
		tests.add(new Object[]{1.2123e-6, "1.212e-06"});
		tests.add(new Object[]{Double.POSITIVE_INFINITY, "Infinity"});
		tests.add(new Object[]{Double.NEGATIVE_INFINITY, "-Infinity"});
		tests.add(new Object[]{Double.NaN, "NaN"});
		return tests.toArray(new Object[][]{});
	}

	@Test(dataProvider = "VCFWriterDoubleFormatTestData")
	public void testVCFWriterDoubleFormatTestData(final double d, final String expected) {
		Assert.assertEquals(VCFEncoder.formatVCFDouble(d), expected, "Failed to pretty print double in VCFWriter");
	}

}
