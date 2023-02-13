/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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

import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import picard.cmdline.CommandLineProgram;
import picard.sam.DuplicationMetrics;
import picard.sam.DuplicationMetricsFactory;
import picard.sam.testers.SamFileTester;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * This class is an extension of SamFileTester used to test AbstractMarkDuplicatesCommandLineProgram's with SAM files generated on the fly.
 * This performs the underlying tests defined by classes such as AbstractMarkDuplicatesCommandLineProgramTest.
 */
abstract public class AbstractMarkDuplicatesCommandLineProgramTester extends SamFileTester {

    final File metricsFile;
    final DuplicationMetrics expectedMetrics;

    boolean testOpticalDuplicateDTTag = false;
    final public Map<List<String>, Double> expectedSetSizeMap = new HashMap<>();

    public AbstractMarkDuplicatesCommandLineProgramTester(final ScoringStrategy duplicateScoringStrategy, SAMFileHeader.SortOrder sortOrder) {
        this(duplicateScoringStrategy, sortOrder, true);
    }

    public AbstractMarkDuplicatesCommandLineProgramTester(final ScoringStrategy duplicateScoringStrategy, SAMFileHeader.SortOrder sortOrder, boolean recordNeedSorting) {
        super(50, true, SAMRecordSetBuilder.DEFAULT_CHROMOSOME_LENGTH, duplicateScoringStrategy, sortOrder, recordNeedSorting);

        expectedMetrics = DuplicationMetricsFactory.createMetrics();
        expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES = 0;

        metricsFile = new File(getOutputDir(), "metrics.txt");
        addArg("METRICS_FILE=" + metricsFile);
        addArg("DUPLICATE_SCORING_STRATEGY=" + duplicateScoringStrategy.name());
        recordOpticalDuplicatesMarked();
    }

    public AbstractMarkDuplicatesCommandLineProgramTester(final ScoringStrategy duplicateScoringStrategy) {
        this(duplicateScoringStrategy, SAMFileHeader.SortOrder.coordinate);
    }

    public AbstractMarkDuplicatesCommandLineProgramTester() {
        this(SAMRecordSetBuilder.DEFAULT_DUPLICATE_SCORING_STRATEGY);
    }

    @Override
    public String getCommandLineProgramName() { return getProgram().getClass().getSimpleName(); }

    /**
     * Tells MarkDuplicates to record which reads are optical duplicates
     *
     * NOTE: this should be overridden as a blank method for inheriting classes where the tested tool doesn't support the 'TAGGING_POLICY' argument
     */
    public void recordOpticalDuplicatesMarked() {
        testOpticalDuplicateDTTag = true;
        addArg("TAGGING_POLICY=" + MarkDuplicates.DuplicateTaggingPolicy.OpticalOnly);
    }

    /**
     * Fill in expected duplication metrics directly from the input records given to this tester
     */
    public void updateExpectedDuplicationMetrics() {

        final FormatUtil formatter = new FormatUtil();

        try (final CloseableIterator<SAMRecord> inputRecordIterator = this.getRecordIterator()) {
            while (inputRecordIterator.hasNext()) {
                final SAMRecord record = inputRecordIterator.next();
                if (record.isSecondaryOrSupplementary()) {
                    ++expectedMetrics.SECONDARY_OR_SUPPLEMENTARY_RDS;
                } else {
                    final String key = samRecordToDuplicatesFlagsKey(record);
                    if (!this.duplicateFlags.containsKey(key)) {
                        System.err.println("DOES NOT CONTAIN KEY: " + key);
                    }
                    final boolean isDuplicate = this.duplicateFlags.get(key);

                    // First bring the simple metricsFile up to date
                    if (record.getReadUnmappedFlag()) {
                        ++expectedMetrics.UNMAPPED_READS;
                    } else if (!record.getReadPairedFlag() || record.getMateUnmappedFlag()) {
                        ++expectedMetrics.UNPAIRED_READS_EXAMINED;
                        if (isDuplicate) {
                            ++expectedMetrics.UNPAIRED_READ_DUPLICATES;
                        }
                    } else {
                        ++expectedMetrics.READ_PAIRS_EXAMINED; // will need to be divided by 2 at the end
                        if (isDuplicate) {
                            ++expectedMetrics.READ_PAIR_DUPLICATES; // will need to be divided by 2 at the end
                        }
                    }
                }
            }
        }
        expectedMetrics.READ_PAIR_DUPLICATES = expectedMetrics.READ_PAIR_DUPLICATES / 2;
        expectedMetrics.READ_PAIRS_EXAMINED = expectedMetrics.READ_PAIRS_EXAMINED / 2;
        expectedMetrics.calculateDerivedFields();

        // Have to run this Double value through the same format/parsing operations as during a file write/read
        expectedMetrics.PERCENT_DUPLICATION = formatter.parseDouble(formatter.format(expectedMetrics.PERCENT_DUPLICATION));
    }

    public void setExpectedOpticalDuplicate(final int expectedOpticalDuplicatePairs) {
        expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES = expectedOpticalDuplicatePairs;
    }


    /**
    * This method is called before iterating through output records in test method.
    * Should be used to update any expectations to be used in implementation
    * specific tests.
     **/
    void updateExpectationsHook() {
    }

    /**
    * This method is called for each record in the output.  Should be used
    * to perform any additional implementation specific tests not included
     * in default test.
     **/
    void testRecordHook(final SAMRecord record){
    }

    /**
     * This method is called after iterating through output records.
     * Should be used to test any implementation specific methods not
     * included in default test
     */
    void testOutputsHook() {
    }

    @Override
    public void test() throws IOException {
        try {
            updateExpectedDuplicationMetrics();
            updateExpectationsHook();

            // Read the output and check the duplicate flag
            int outputRecords = 0;
            final Set<String> sequencingDTErrorsSeen = new HashSet<>();
            try(final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(fastaFiles.get(samRecordSetBuilder.getHeader())).open(getOutput())) {
                for (final SAMRecord record : reader) {
                    outputRecords++;
                    testRecordHook(record);
                    final String key = samRecordToDuplicatesFlagsKey(record);
                    Assert.assertTrue(this.duplicateFlags.containsKey(key), "DOES NOT CONTAIN KEY: " + key);
                    final boolean value = this.duplicateFlags.get(key);
                    this.duplicateFlags.remove(key);
                    Assert.assertEquals(record.getDuplicateReadFlag(), value, "Mismatching read: " + record.getSAMString());
                    if (testOpticalDuplicateDTTag && MarkDuplicates.DUPLICATE_TYPE_SEQUENCING.equals(record.getAttribute("DT"))) {
                        sequencingDTErrorsSeen.add(record.getReadName());
                    }
                    Assert.assertEquals(record.getDuplicateReadFlag(), value);
                }
            }

            // Ensure the program output the same number of records as were read in
            Assert.assertEquals(outputRecords, this.getNumberOfRecords(), ("saw " + outputRecords + " output records, vs. " + this.getNumberOfRecords() + " input records"));

            // Check the values written to metrics.txt against our input expectations
            final MetricsFile<DuplicationMetrics, Double> metricsOutput = new MetricsFile<>();
            try{
                metricsOutput.read(new FileReader(metricsFile));
            }
            catch (final FileNotFoundException ex) {
                Assert.fail("Metrics file not found: " + ex.getMessage());
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
            if (testOpticalDuplicateDTTag) {
                Assert.assertEquals(sequencingDTErrorsSeen.size(), expectedMetrics.READ_PAIR_OPTICAL_DUPLICATES, "READ_PAIR_OPTICAL_DUPLICATES does not match duplicate groups observed in the file");
                Assert.assertEquals(sequencingDTErrorsSeen.size(), observedMetrics.READ_PAIR_OPTICAL_DUPLICATES, "READ_PAIR_OPTICAL_DUPLICATES does not match duplicate groups observed in the file");
            }

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

            testOutputsHook();
        } finally {
            IOUtil.recursiveDelete(getOutputDir().toPath());
        }
    }

    abstract protected CommandLineProgram getProgram();
}
