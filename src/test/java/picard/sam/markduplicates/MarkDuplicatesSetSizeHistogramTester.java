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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.TestUtil;
import org.apache.commons.jexl2.UnifiedJEXL;
import org.testng.Assert;
import org.testng.annotations.AfterClass;
import picard.cmdline.CommandLineProgram;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

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
        updateExpectedDuplicationMetrics();

        // Read the output and check the duplicate flag
        int outputRecords = 0;

        try (SamReader reader = SamReaderFactory.makeDefault().open(getOutput())) {
            for (final SAMRecord record : reader) {
                outputRecords++;
                final String key = samRecordToDuplicatesFlagsKey(record);
                if (!this.duplicateFlags.containsKey(key)) {
                    System.err.println("DOES NOT CONTAIN KEY: " + key);
                }
                Assert.assertTrue(this.duplicateFlags.containsKey(key));
                final boolean value = this.duplicateFlags.get(key);
                this.duplicateFlags.remove(key);
                if (value != record.getDuplicateReadFlag()) {

                    System.err.println("Mismatching read:");
                    System.err.println(record.getSAMString());
                }
                Assert.assertEquals(record.getDuplicateReadFlag(), value);
            }
        }

        // Ensure the program output the same number of records as were read in
        Assert.assertEquals(outputRecords, this.getNumberOfRecords());

        // Check the values written to metrics.txt against our input expectations
        final MetricsFile<DuplicationMetrics, Double> metricsOutput = new MetricsFile<DuplicationMetrics, Double>();
        try {
            metricsOutput.read(new FileReader(metricsFile));
        } catch (final FileNotFoundException ex) {
            // Rethrows as unchecked exception to avoid having to change existing tests
            throw new RuntimeException(ex);
        }
        Assert.assertEquals(metricsOutput.getMetrics().size(), 1);
        final DuplicationMetrics observedMetrics = metricsOutput.getMetrics().get(0);
        Assert.assertEquals(observedMetrics.UNPAIRED_READS_EXAMINED, expectedMetrics.UNPAIRED_READS_EXAMINED, "UNPAIRED_READS_EXAMINED does not match expected");
        Assert.assertEquals(observedMetrics.READ_PAIRS_EXAMINED, expectedMetrics.READ_PAIRS_EXAMINED, "READ_PAIRS_EXAMINED does not match expected");
        Assert.assertEquals(observedMetrics.UNMAPPED_READS, expectedMetrics.UNMAPPED_READS, "UNMAPPED_READS does not match expected");
        Assert.assertEquals(observedMetrics.UNPAIRED_READ_DUPLICATES, expectedMetrics.UNPAIRED_READ_DUPLICATES, "UNPAIRED_READ_DUPLICATES does not match expected");
        Assert.assertEquals(observedMetrics.READ_PAIR_DUPLICATES, expectedMetrics.READ_PAIR_DUPLICATES, "READ_PAIR_DUPLICATES does not match expected");
        Assert.assertEquals(observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES, "READ_PAIR_OPTICAL_DUPLICATES does not match expected");
        Assert.assertEquals(observedMetrics.PERCENT_DUPLICATION, expectedMetrics.PERCENT_DUPLICATION, "PERCENT_DUPLICATION does not match expected");
        Assert.assertEquals(observedMetrics.ESTIMATED_LIBRARY_SIZE, expectedMetrics.ESTIMATED_LIBRARY_SIZE, "ESTIMATED_LIBRARY_SIZE does not match expected");
        Assert.assertEquals(observedMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS, expectedMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS, "SECONDARY_OR_SUPPLEMENTARY_RDS does not match expected");

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
                throw new RuntimeException("Could not not find matching entry for expectedSetSizeMap in metrics.");
            }
        }
    }

    @AfterClass
    public void afterTest() {
        TestUtil.recursiveDelete(getOutputDir());
    }

}
