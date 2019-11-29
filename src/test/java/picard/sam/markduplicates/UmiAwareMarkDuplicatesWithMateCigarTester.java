/*
 * The MIT License
 *
 * Copyright (c) 2016 The Broad Institute
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
import org.apache.commons.lang3.StringUtils;
import org.testng.Assert;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test
 * AbstractMarkDuplicatesCommandLineProgram's with SAM files generated on the fly.  This performs the underlying tests
 * defined by classes such as AbstractMarkDuplicatesCommandLineProgramTest.
 * @author fleharty
 */

public class UmiAwareMarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicatesCommandLineProgramTester {
    private int readNameCounter = 0;
    private List<String> expectedAssignedUmis;
    private UmiMetrics expectedMetrics;
    private File umiMetricsFile = new File(getOutputDir(), "umi_metrics.txt");

    // This tag is only used for testing, it indicates what we expect to see in the inferred UMI tag.
    private final String expectedUmiTag = "RE";

    // This default constructor is intended to be used by tests inherited from
    // AbstractMarkDuplicatesCommandLineProgramTester.  Since those tests use
    // reads that don't have UMIs we enable the ALLOW_MISSING_UMIS option.
    UmiAwareMarkDuplicatesWithMateCigarTester() {
        addArg("UMI_METRICS_FILE=" + umiMetricsFile);
        addArg("ALLOW_MISSING_UMIS=" + true);
    }

    UmiAwareMarkDuplicatesWithMateCigarTester(final boolean allowMissingUmis) {
        addArg("UMI_METRICS_FILE=" + umiMetricsFile);

        if (allowMissingUmis) {
            addArg("ALLOW_MISSING_UMIS=" + true);
        }
    }

    @Override
    public void recordOpticalDuplicatesMarked() {}

    public void addMatePairWithUmi(final String library, final String umi, final String assignedUMI, final boolean isDuplicate1, final boolean isDuplicate2) {

        final String readName = "READ" + readNameCounter++;
        final String cigar1 = null;
        final String cigar2 = null;
        final boolean strand1 = false;
        final boolean strand2 = true;

        final int referenceSequenceIndex1 = 0;
        final int referenceSequenceIndex2 = 0;
        final int alignmentStart1 = 20;
        final int alignmentStart2 = 20;

        final boolean record1Unmapped = false;
        final boolean record2Unmapped = false;

        final boolean firstOnly = false;
        final boolean record1NonPrimary = false;
        final boolean record2NonPrimary = false;

        final int defaultQuality = 10;

        addMatePairWithUmi(library, readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2, record1Unmapped,
                record2Unmapped, isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary,
                defaultQuality, umi, assignedUMI);

    }

    public void addMatePairWithUmi(final String library, final String umi, final String assignedUMI, final boolean isDuplicate1, final boolean isDuplicate2,
                                   final int alignmentStart1, final int alignmentStart2, final boolean strand1, final boolean strand2) {

        final String readName = "READ" + readNameCounter++;
        final String cigar1 = null;
        final String cigar2 = null;

        final int referenceSequenceIndex1 = 0;
        final int referenceSequenceIndex2 = 0;

        final boolean record1Unmapped = false;
        final boolean record2Unmapped = false;

        final boolean firstOnly = false;
        final boolean record1NonPrimary = false;
        final boolean record2NonPrimary = false;

        final int defaultQuality = 10;

        addMatePairWithUmi(library, readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2, record1Unmapped,
                record2Unmapped, isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary,
                defaultQuality, umi, assignedUMI);
    }

    public void addMatePairWithUmi(final String library,
                                   final String readName,
                                   final int referenceSequenceIndex1,
                                   final int referenceSequenceIndex2,
                                   final int alignmentStart1,
                                   final int alignmentStart2,
                                   final boolean record1Unmapped,
                                   final boolean record2Unmapped,
                                   final boolean isDuplicate1,
                                   final boolean isDuplicate2,
                                   final String cigar1,
                                   final String cigar2,
                                   final boolean strand1,
                                   final boolean strand2,
                                   final boolean firstOnly,
                                   final boolean record1NonPrimary,
                                   final boolean record2NonPrimary,
                                   final int defaultQuality,
                                   final String umi,
                                   final String assignedUMI) {
        final List<SAMRecord> samRecordList = samRecordSetBuilder.addPair(readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2,
                record1Unmapped, record2Unmapped, cigar1, cigar2, strand1, strand2, record1NonPrimary, record2NonPrimary, defaultQuality);

        final SAMRecord record1 = samRecordList.get(0);
        final SAMRecord record2 = samRecordList.get(1);

        record1.getReadGroup().setLibrary(library);
        record2.getReadGroup().setLibrary(library);
        if (this.noMateCigars) {
            record1.setAttribute("MC", null);
            record2.setAttribute("MC", null);
        }

        if (firstOnly) {
            samRecordSetBuilder.getRecords().remove(record2);
        }

        final String key1 = samRecordToDuplicatesFlagsKey(record1);
        Assert.assertFalse(this.duplicateFlags.containsKey(key1));
        this.duplicateFlags.put(key1, isDuplicate1);

        final String key2 = samRecordToDuplicatesFlagsKey(record2);
        Assert.assertFalse(this.duplicateFlags.containsKey(key2));
        this.duplicateFlags.put(key2, isDuplicate2);

        if (umi != null) {
            // TODO: Replace "RX" with SAMTag.RX once this tag is available in HTSJDK
            record1.setAttribute("RX", umi);
            record2.setAttribute("RX", umi);
        }
        if (assignedUMI != null) {
            // Set the expected UMI, this is a special tag used only for testing.
            record1.setAttribute(expectedUmiTag, assignedUMI);
            record2.setAttribute(expectedUmiTag, assignedUMI);
        }
    }

    UmiAwareMarkDuplicatesWithMateCigarTester setExpectedAssignedUmis(final List<String> expectedAssignedUmis) {
        this.expectedAssignedUmis = expectedAssignedUmis;
        return this;
    }

    UmiAwareMarkDuplicatesWithMateCigarTester setExpectedMetrics(final UmiMetrics expectedMetrics) {
        this.expectedMetrics = expectedMetrics;
        return this;
    }

    @Override
    public void test() {
        final SamReader reader = SamReaderFactory.makeDefault().open(getOutput());
        for (final SAMRecord record : reader) {
            // If there are expected assigned UMIs, check to make sure they match
            if (expectedAssignedUmis != null) {
                Assert.assertEquals(getAssignedUmi(record.getStringAttribute("MI")), record.getAttribute(expectedUmiTag));
            }
        }

        if (expectedMetrics != null) {
            // Check the values written to metrics.txt against our input expectations
            final MetricsFile<UmiMetrics, Comparable<?>> metricsOutput = new MetricsFile<UmiMetrics, Comparable<?>>();
            try {
                metricsOutput.read(new FileReader(umiMetricsFile));
            }
            catch (final FileNotFoundException ex) {
                System.err.println("Metrics file not found: " + ex);
            }
            final double tolerance = 1e-6;
            Assert.assertEquals(metricsOutput.getMetrics().size(), 1);
            final UmiMetrics observedMetrics = metricsOutput.getMetrics().get(0);

            Assert.assertEquals(observedMetrics.LIBRARY, expectedMetrics.LIBRARY, "LIBRARY does not match expected");
            Assert.assertEquals(observedMetrics.MEAN_UMI_LENGTH, expectedMetrics.MEAN_UMI_LENGTH, "UMI_LENGTH does not match expected");
            Assert.assertEquals(observedMetrics.OBSERVED_UNIQUE_UMIS, expectedMetrics.OBSERVED_UNIQUE_UMIS, "OBSERVED_UNIQUE_UMIS does not match expected");
            Assert.assertEquals(observedMetrics.INFERRED_UNIQUE_UMIS, expectedMetrics.INFERRED_UNIQUE_UMIS, "INFERRED_UNIQUE_UMIS does not match expected");
            Assert.assertEquals(observedMetrics.OBSERVED_BASE_ERRORS, expectedMetrics.OBSERVED_BASE_ERRORS, "OBSERVED_BASE_ERRORS does not match expected");
            Assert.assertEquals(observedMetrics.DUPLICATE_SETS_IGNORING_UMI, expectedMetrics.DUPLICATE_SETS_IGNORING_UMI, "DUPLICATE_SETS_IGNORING_UMI does not match expected");
            Assert.assertEquals(observedMetrics.DUPLICATE_SETS_WITH_UMI, expectedMetrics.DUPLICATE_SETS_WITH_UMI, "DUPLICATE_SETS_WITH_UMI does not match expected");
            Assert.assertEquals(observedMetrics.INFERRED_UMI_ENTROPY, expectedMetrics.INFERRED_UMI_ENTROPY, tolerance, "INFERRED_UMI_ENTROPY does not match expected");
            Assert.assertEquals(observedMetrics.OBSERVED_UMI_ENTROPY, expectedMetrics.OBSERVED_UMI_ENTROPY, tolerance, "OBSERVED_UMI_ENTROPY does not match expected");
            Assert.assertEquals(observedMetrics.UMI_BASE_QUALITIES, expectedMetrics.UMI_BASE_QUALITIES, tolerance, "UMI_BASE_QUALITIES does not match expected");
            Assert.assertEquals(observedMetrics.PCT_UMI_WITH_N, expectedMetrics.PCT_UMI_WITH_N, tolerance,"PERCENT_UMI_WITH_N does not match expected" );
        }

        // Also do tests from AbstractMarkDuplicatesCommandLineProgramTester
        try {
            super.test();
        } catch (IOException ex) {
            Assert.fail("Could not open metrics file: ", ex);
        }
    }

    @Override
    protected CommandLineProgram getProgram() {
        UmiAwareMarkDuplicatesWithMateCigar uamdwmc = new UmiAwareMarkDuplicatesWithMateCigar();
        return uamdwmc;
    }


    /**
     *
     * @param molecularIndex
     * @return
     */
    private String getAssignedUmi(final String molecularIndex) {
        if (molecularIndex == null) {
            return null;
        }

        if (StringUtils.countMatches(molecularIndex, UmiUtil.UMI_NAME_SEPARATOR) == 2) {
            return StringUtils.substringBetween(molecularIndex, UmiUtil.UMI_NAME_SEPARATOR, UmiUtil.UMI_NAME_SEPARATOR);
        } else {
            return StringUtils.substringAfter(molecularIndex, UmiUtil.UMI_NAME_SEPARATOR);
        }
    }

}
