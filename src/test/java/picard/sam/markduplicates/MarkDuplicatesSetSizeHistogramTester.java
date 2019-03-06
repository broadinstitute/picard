/*
 * The MIT License
 *
 * Copyright (c) 2018 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.sam.markduplicates;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import picard.cmdline.CommandLineProgram;
import picard.sam.DuplicationMetrics;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test MarkDuplicatesWithMateCigar with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as see AbstractMarkDuplicatesCommandLineProgramTest and MarkDuplicatesWithMateCigarTest.
 */
public class MarkDuplicatesSetSizeHistogramTester extends AbstractMarkDuplicatesCommandLineProgramTester {

    final public Map<List<String>, Double> expectedSetSizeMap = new HashMap<>(); // key=(Histogram Label, histogram bin), value=histogram entry

    @Override
    protected CommandLineProgram getProgram() {
        return new MarkDuplicates();
    }

    @Override
    public void test() throws IOException {
        final MetricsFile<DuplicationMetrics, Double> metricsOutput = testMetrics();

        // Check contents of set size bin against expected values
        if (!expectedSetSizeMap.isEmpty()) {
            boolean checked = false;
            for (final Histogram<Double> histo : metricsOutput.getAllHistograms()) {
                final String label = histo.getValueLabel();
                for (final Double bin : histo.keySet()) {
                    final String binStr = String.valueOf(bin);
                    final List<String> labelBinStr = Arrays.asList(label, binStr);
                    if (expectedSetSizeMap.containsKey(labelBinStr)) {
                        checked = true;
                        Histogram.Bin<Double> binValue = histo.get(bin);
                        final double actual = binValue.getValue();
                        final double expected = expectedSetSizeMap.get(labelBinStr);
                        Assert.assertEquals(actual, expected);
                    }
                }
            }
            if (!checked) {
                Assert.fail("Could not not find matching entry for expectedSetSizeMap in metrics.");
            }
        }
    }

    @AfterClass
    public void afterTest() {
        IOUtil.recursiveDelete(getOutputDir().toPath());
    }

}
