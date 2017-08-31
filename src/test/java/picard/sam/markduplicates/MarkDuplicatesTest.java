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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesTester (see getTester).
 */
public class MarkDuplicatesTest extends AbstractMarkDuplicatesCommandLineProgramTest {
    protected static String TEST_BASE_NAME = null;
    protected static File TEST_DATA_DIR = null;

    @BeforeClass
    public void setUp() {
        TEST_BASE_NAME = "MarkDuplicates";
        TEST_DATA_DIR = new File("testdata/picard/sam/MarkDuplicates");
    }

    protected AbstractMarkDuplicatesCommandLineProgramTester getTester() {
        return new MarkDuplicatesTester();
    }

    // NB: this test should return different results than MarkDuplicatesWithMateCigar
    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnly() {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(76);
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(0, 12, 46, false, false, "6S42M28S", "3S73M", true, 50); // only add the first one
        // NB: this next record should not be a duplicate in MarkDuplicates
        tester.addMappedPair(0, 12, 51, false, false, "6S42M28S", "8S68M", true, 50); // only add the first one
        tester.runTest();
    }

    /**
     * Test that PG header records are created & chained appropriately (or not created), and that the PG record chains
     * are as expected.  MarkDuplicates is used both to merge and to mark dupes in this case.
     * @param suppressPg If true, do not create PG header record.
     * @param expectedPnVnByReadName For each read, info about the expect chain of PG records.
     */
    @Test(dataProvider = "pgRecordChainingTest")
    public void pgRecordChainingTest(final boolean suppressPg,
                                     final Map<String, List<ExpectedPnAndVn>> expectedPnVnByReadName) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        try {
            // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
            // record creation according to suppressPg.
            final MarkDuplicates markDuplicates = new MarkDuplicates();
            final ArrayList<String> args = new ArrayList<String>();
            for (int i = 1; i <= 3; ++i) {
                args.add("INPUT=" + new File(TEST_DATA_DIR, "merge" + i + ".sam").getAbsolutePath());
            }
            final File outputSam = new File(outputDir, TEST_BASE_NAME + ".sam");
            args.add("OUTPUT=" + outputSam.getAbsolutePath());
            args.add("METRICS_FILE=" + new File(outputDir, TEST_BASE_NAME + ".duplicate_metrics").getAbsolutePath());
            args.add("ADD_PG_TAG_TO_READS=true");
            if (suppressPg) args.add("PROGRAM_RECORD_ID=null");

            // I generally prefer to call doWork rather than invoking the argument parser, but it is necessary
            // in this case to initialize the command line.
            // Note that for the unit test, version won't come through because it is obtained through jar
            // manifest, and unit test doesn't run code from a jar.
            Assert.assertEquals(markDuplicates.instanceMain(args.toArray(new String[args.size()])), 0);

            // Read the MarkDuplicates output file, and get the PG ID for each read.  In this particular test,
            // the PG ID should be the same for both ends of a pair.
            final SamReader reader = SamReaderFactory.makeDefault().open(outputSam);

            final Map<String, String> pgIdForReadName = new HashMap<String, String>();
            for (final SAMRecord rec : reader) {
                final String existingPgId = pgIdForReadName.get(rec.getReadName());
                final String thisPgId = rec.getStringAttribute(SAMTag.PG.name());
                if (existingPgId != null) {
                    Assert.assertEquals(thisPgId, existingPgId);
                } else {
                    pgIdForReadName.put(rec.getReadName(), thisPgId);
                }
            }
            final SAMFileHeader header = reader.getFileHeader();
            CloserUtil.close(reader);

            // Confirm that for each read name, the chain of PG records contains exactly the number that is expected,
            // and that values in the PG chain are as expected.
            for (final Map.Entry<String, List<ExpectedPnAndVn>> entry : expectedPnVnByReadName.entrySet()) {
                final String readName = entry.getKey();
                final List<ExpectedPnAndVn> expectedList = entry.getValue();
                String pgId = pgIdForReadName.get(readName);
                for (final ExpectedPnAndVn expected : expectedList) {
                    final SAMProgramRecord programRecord = header.getProgramRecord(pgId);
                    if (expected.expectedPn != null) Assert.assertEquals(programRecord.getProgramName(), expected.expectedPn);
                    if (expected.expectedVn != null) Assert.assertEquals(programRecord.getProgramVersion(), expected.expectedVn);
                    pgId = programRecord.getPreviousProgramGroupId();
                }
                Assert.assertNull(pgId);
            }

        } finally {
            TestUtil.recursiveDelete(outputDir);
        }
    }

    /**
     * Represents an expected PN value and VN value for a PG record.  If one of thexe is null, any value is allowed
     * in the PG record being tested.
     */
    private static class ExpectedPnAndVn {
        final String expectedPn;
        final String expectedVn;

        private ExpectedPnAndVn(final String expectedPn, final String expectedVn) {
            this.expectedPn = expectedPn;
            this.expectedVn = expectedVn;
        }
    }

    @DataProvider(name = "pgRecordChainingTest")
    public Object[][] pgRecordChainingTestDataProvider() {
        // Two test cases: One in which PG record generation is enabled, the other in which it is turned off.
        final Map<String, List<ExpectedPnAndVn>> withPgMap = new HashMap<String, List<ExpectedPnAndVn>>();
        withPgMap.put("1AAXX.1.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, null), new ExpectedPnAndVn(TEST_BASE_NAME, "1"), new ExpectedPnAndVn("bwa", "1")));
        withPgMap.put("1AAXX.2.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, null), new ExpectedPnAndVn("bwa", "2")));
        withPgMap.put("1AAXX.3.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, null)));

        final Map<String, List<ExpectedPnAndVn>> suppressPgMap = new HashMap<String, List<ExpectedPnAndVn>>();
        suppressPgMap .put("1AAXX.1.1", Arrays.asList(new ExpectedPnAndVn(TEST_BASE_NAME, "1"), new ExpectedPnAndVn("bwa", "1")));
        suppressPgMap .put("1AAXX.2.1", Arrays.asList(new ExpectedPnAndVn("bwa", "2")));
        suppressPgMap .put("1AAXX.3.1", new ArrayList<ExpectedPnAndVn>(0));
        return new Object[][] {
                { false, withPgMap},
                { true, suppressPgMap}
        };
    }

    @Test(dataProvider = "testOpticalDuplicateDetectionDataProvider")
    public void testOpticalDuplicateDetection(final File sam, final long expectedNumOpticalDuplicates) {
        final File outputDir = IOUtil.createTempDir(TEST_BASE_NAME + ".", ".tmp");
        outputDir.deleteOnExit();
        final File outputSam = new File(outputDir, TEST_BASE_NAME + ".sam");
        outputSam.deleteOnExit();
        final File metricsFile = new File(outputDir, TEST_BASE_NAME + ".duplicate_metrics");
        metricsFile.deleteOnExit();
        // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
        // record creation according to suppressPg.
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.setupOpticalDuplicateFinder();
        markDuplicates.INPUT = CollectionUtil.makeList(sam.getAbsolutePath());
        markDuplicates.OUTPUT = outputSam;
        markDuplicates.METRICS_FILE = metricsFile;
        markDuplicates.TMP_DIR = CollectionUtil.makeList(outputDir);
        // Needed to suppress calling CommandLineProgram.getVersion(), which doesn't work for code not in a jar
        markDuplicates.PROGRAM_RECORD_ID = null;
        Assert.assertEquals(markDuplicates.doWork(), 0);
        Assert.assertEquals(markDuplicates.numOpticalDuplicates(), expectedNumOpticalDuplicates);
        TestUtil.recursiveDelete(outputDir);

    }

    @DataProvider(name="testOpticalDuplicateDetectionDataProvider")
    public Object[][] testOpticalDuplicateDetectionDataProvider() {
        return new Object[][] {
                {new File(TEST_DATA_DIR, "optical_dupes.sam"), 1L},
                {new File(TEST_DATA_DIR, "optical_dupes_casava.sam"), 1L},
        };
    }

    @Test
    public void testWithBarcodeFragmentDuplicate() {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.addMappedFragment(2, 41212324, false, "50M", DEFAULT_BASE_QUALITY);
        tester.addMappedFragment(2, 41212324, true, "50M", DEFAULT_BASE_QUALITY);
        final String barcodeTag = "BC";
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(tester.getRecordIterator())) {
            record.setAttribute(barcodeTag, "Barcode1");
        }
        tester.addArg("BARCODE_TAG=" + barcodeTag);
        tester.runTest();
    }

    public void addMappedFragment(final int referenceSequenceIndex, final int alignmentStart, final boolean isDuplicate, final String cigar,
                                  final int defaultQualityScore) {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.addMatePair("RUNID:1:1:15993:13361", 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:2:2:15993:13362", 2, 41212324, 41212310, false, false, true, true, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        final String barcodeTag = "BC";
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(tester.getRecordIterator())) {
            record.setAttribute(barcodeTag, "Barcode1");
        }
        tester.addArg("BARCODE_TAG=" + barcodeTag);
        tester.runTest();
    }

    @Test
    public void testWithBarcodeDuplicate() {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);
        tester.addMatePair("RUNID:1:1:15993:13361", 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair("RUNID:2:2:15993:13362", 2, 41212324, 41212310, false, false, true, true, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        final String barcodeTag = "BC";
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(tester.getRecordIterator())) {
            record.setAttribute(barcodeTag, "Barcode1");
        }
        tester.addArg("BARCODE_TAG=" + barcodeTag);
        tester.runTest();
    }

    @Test
    public void testWithBarcodeComplex() {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);
        final String readNameOne = "RUNID:1:1:15993:13361";
        final String readNameTwo = "RUNID:2:2:15993:13362";
        final String readNameThree = "RUNID:3:3:15993:13362";

        // first two reads have the same barcode, third read has a different barcode
        tester.addMatePair(readNameOne, 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair(readNameTwo, 2, 41212324, 41212310, false, false, true, true, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY); // same barcode as the first
        tester.addMatePair(readNameThree, 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);

        final String barcodeTag = "BC";
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(tester.getRecordIterator())) {
            if (record.getReadName().equals(readNameOne) || record.getReadName().equals(readNameTwo)) {
                record.setAttribute(barcodeTag, "Barcode1");
            }
            else if (record.getReadName().equals(readNameThree)) {
                record.setAttribute(barcodeTag, "Barcode2");
            }
        }
        tester.addArg("BARCODE_TAG=" + barcodeTag);
        tester.runTest();
    }

    @Test
    public void testWithIndividualReadBarcodes() {
        final AbstractMarkDuplicatesCommandLineProgramTester tester = getTester();
        tester.getSamRecordSetBuilder().setReadLength(68);
        final String readNameOne = "RUNID:1:1:15993:13361";
        final String readNameTwo = "RUNID:2:2:15993:13362";
        final String readNameThree = "RUNID:3:3:15993:13362";

        // first two reads have the same barcode (all three), third read has a different barcode for the second end
        tester.addMatePair(readNameOne, 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.addMatePair(readNameTwo, 2, 41212324, 41212310, false, false, true, true, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY); // same barcode as the first
        tester.addMatePair(readNameThree, 2, 41212324, 41212310, false, false, false, false, "33S35M", "19S49M", true, true, false, false, false, DEFAULT_BASE_QUALITY);

        final String barcodeTag = "BC";
        final String readOneBarcodeTag = "BX"; // want the same tag as the second end, since this is allowed
        final String readTwoBarcodeTag = "BX";
        for (final SAMRecord record : new IterableAdapter<SAMRecord>(tester.getRecordIterator())) {
            record.setAttribute(barcodeTag, "Barcode1"); // same barcode
            if (record.getFirstOfPairFlag()) { // always the same value for the first end
                record.setAttribute(readOneBarcodeTag, "readOne1");
            }
            else { // second end
                if (record.getReadName().equals(readNameOne) || record.getReadName().equals(readNameTwo)) {
                    record.setAttribute(readTwoBarcodeTag, "readTwo1");
                } else if (record.getReadName().equals(readNameThree)) {
                    record.setAttribute(readTwoBarcodeTag, "readTwo2");
                }
            }
        }
        tester.addArg("BARCODE_TAG=" + barcodeTag);
        tester.addArg("READ_ONE_BARCODE_TAG=" + readOneBarcodeTag);
        tester.addArg("READ_TWO_BARCODE_TAG=" + readTwoBarcodeTag);

        tester.runTest();
    }
}
