/*
 * The MIT License
 *
 * Copyright (c) 2012 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.testers.MarkDuplicatesTester;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class MarkDuplicatesTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/MarkDuplicates");


    /**
     * Test that PG header records are created & chained appropriately (or not created), and that the PG record chains
     * are as expected.  MarkDuplicates is used both to merge and to mark dupes in this case.
     * @param suppressPg If true, do not create PG header record.
     * @param expectedPnVnByReadName For each read, info about the expect chain of PG records.
     */
    @Test(dataProvider = "pgRecordChainingTest")
    public void pgRecordChainingTest(final boolean suppressPg,
                                     final Map<String, List<ExpectedPnAndVn>> expectedPnVnByReadName) {
        final File outputDir = IOUtil.createTempDir("MarkDuplicatesTest.", ".tmp");
        outputDir.deleteOnExit();
        try {
            // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
            // record creation according to suppressPg.
            final MarkDuplicates markDuplicates = new MarkDuplicates();
            final ArrayList<String> args = new ArrayList<String>();
            for (int i = 1; i <= 3; ++i) {
                args.add("INPUT=" + new File(TEST_DATA_DIR, "merge" + i + ".sam").getAbsolutePath());
            }
            final File outputSam = new File(outputDir, "markDuplicatesTest.sam");
            args.add("OUTPUT=" + outputSam.getAbsolutePath());
            args.add("METRICS_FILE=" + new File(outputDir, "markDuplicatesTest.duplicate_metrics").getAbsolutePath());
            if (suppressPg) args.add("PROGRAM_RECORD_ID=null");

            // I generally prefer to call doWork rather than invoking the argument parser, but it is necessary
            // in this case to initialize the command line.
            // Note that for the unit test, version won't come through because it is obtained through jar
            // manifest, and unit test doesn't run code from a jar.
            Assert.assertEquals(markDuplicates.instanceMain(args.toArray(new String[args.size()])), 0);

            // Read the MarkDuplicates output file, and get the PG ID for each read.  In this particular test,
            // the PG ID should be the same for both ends of a pair.
            final SAMFileReader reader = new SAMFileReader(outputSam);

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
            reader.close();

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
        withPgMap.put("1AAXX.1.1", Arrays.asList(new ExpectedPnAndVn("MarkDuplicates", null), new ExpectedPnAndVn("MarkDuplicates", "1"), new ExpectedPnAndVn("bwa", "1")));
        withPgMap.put("1AAXX.2.1", Arrays.asList(new ExpectedPnAndVn("MarkDuplicates", null), new ExpectedPnAndVn("bwa", "2")));
        withPgMap.put("1AAXX.3.1", Arrays.asList(new ExpectedPnAndVn("MarkDuplicates", null)));


        final Map<String, List<ExpectedPnAndVn>> suppressPgMap = new HashMap<String, List<ExpectedPnAndVn>>();
        suppressPgMap .put("1AAXX.1.1", Arrays.asList(new ExpectedPnAndVn("MarkDuplicates", "1"), new ExpectedPnAndVn("bwa", "1")));
        suppressPgMap .put("1AAXX.2.1", Arrays.asList(new ExpectedPnAndVn("bwa", "2")));
        suppressPgMap .put("1AAXX.3.1", new ArrayList<ExpectedPnAndVn>(0));
        return new Object[][] {
                { false, withPgMap},
                { true, suppressPgMap}
        };
    }


    @Test
    public void testSingleUnmappedFragment() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addUnmappedFragment(-1, 50);
        tester.runTest();
    }

    @Test
    public void testSingleUnmappedPair() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addUnmappedPair(-1, 50);
        tester.runTest();
    }


    @Test
    public void testSingleMappedFragment() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedFragment(1, 1, false, 50);
        tester.runTest();
    }

    @Test
    public void testSingleMappedPair() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 1, 100, false, false, 50);
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndSingleMappedPair() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedFragment(1, 1, true, 30); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, 50);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairs() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 1, 100, false, false, 50);
        tester.addMappedPair(1, 1, 100, true, true, 30); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsOpticalDupeDetectionDisabled() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addArg("READ_NAME_REGEX=null");
        tester.addMappedPair(1, 1, 100, false, false, 50);
        tester.addMappedPair(1, 1, 100, true, true, 30); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairs() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 1, 100, false, false, 50);
        tester.addMappedPair(1, 1, 100, true, true, 30); // duplicate!!!
        tester.addMappedPair(1, 1, 100, true, true, 30); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testSingleMappedFragmentAndTwoMappedPairs() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedFragment(1, 1, true, 30); // duplicate!!!
        tester.addMappedPair(1, 1, 100, false, false, 50);
        tester.addMappedPair(1, 1, 100, true, true, 30); // duplicate!!!
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsMatesSoftClipped() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 10022, 10051, false, false, "76M", "8S68M", false, true, false, 50);
        tester.addMappedPair(1, 10022, 10063, false, false, "76M", "5S71M", false, true, false, 50);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClipping() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 2, 46, false, false, "6S42M28S", "3S73M", false, 50);
        tester.addMappedPair(1, 2, 51, true, true, "6S42M28S", "8S68M", false, 50);
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingFirstOfPairOnlyNoMateCigar() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.setNoMateCigars(true);
        // NB: no duplicates
        // 5'1: 2, 5'2:46+73M=118
        // 5'1: 2, 5'2:51+68M=118
        tester.addMappedPair(1, 12, 46, false, false, "6S42M28S", null, true, 50); // only add the first one
        tester.addMappedPair(1, 12, 51, false, false, "6S42M28S", null, true, 50); // only add the first one
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithSoftClippingBoth() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 10046, 10002, false, false, "3S73M", "6S42M28S", true, false, false, 50);
        tester.addMappedPair(1, 10051, 10002, true, true, "8S68M", "6S48M22S", true, false, false, 50);
        tester.runTest();
    }

    @Test
    public void testMatePairFirstUnmapped() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, 50);
        tester.runTest();
    }

    @Test
    public void testMatePairSecondUnmapped() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMatePair(1, 10056, 10056, true, false, false, false, null, "54M22S", false, false, false, 50);
        tester.runTest();
    }

    @Test
    public void testMappedFragmentAndMatePairOneUnmapped() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMatePair(1, 10049, 10049, false, true, false, false, "11M2I63M", null, false, false, false, 50);
        tester.addMappedFragment(1, 10049, true, 30); // duplicate
        tester.runTest();
    }

    @Test
    public void testMappedPairAndMatePairOneUnmapped() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMatePair(1, 10040, 10040, false, true, true, false, "76M", null, false, false, false, 30); // first a duplicate,
        // second end unmapped
        tester.addMappedPair(1, 10189, 10040, false, false, "41S35M", "65M11S", true, false, false, 50); // mapped OK
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientations() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 10182, 10038, false, false, "32S44M", "66M10S", true, false, false, 50); // -/+
        tester.addMappedPair(1, 10038, 10182, true, true, "70M6S", "32S44M", false, true, false, 50); // +/-, both are duplicates
        tester.runTest();
    }

    @Test
    public void testTwoMappedPairsWithOppositeOrientationsNumberTwo() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 10038, 10182, false, false, "70M6S", "32S44M", false, true, false, 50); // +/-, both are duplicates
        tester.addMappedPair(1, 10182, 10038, true, true, "32S44M", "66M10S", true, false, false, 50); // -/+
        tester.runTest();
    }

    @Test
    public void testThreeMappedPairsWithMatchingSecondMate() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        // Read0 and Read2 are duplicates
        // 10181+35=10216, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "41S35M", "47M29S", true, false, false, 50); // -/+
        // 10181+37=10218, 10058
        tester.addMappedPair(1, 10181, 10058, false, false, "37S39M", "44M32S", true, false, false, 50); // -/+
        // 10180+36=10216, 10058
        tester.addMappedPair(1, 10180, 10058, true, true, "36S40M", "50M26S", true, false, false, 50); // -/+, both are duplicates
        tester.runTest();
    }

    @Test
    public void testMappedPairWithSamePosition() {
        final MarkDuplicatesTester tester = new MarkDuplicatesTester();
        tester.addMappedPair(1, 4914, 4914, false, false, "37M39S", "73M3S", false, false, false, 50); // +/+
        tester.runTest();
    }

    @Test(dataProvider = "testOpticalDuplicateDetectionDataProvider")
    public void testOpticalDuplicateDetection(final File sam, final long expectedNumOpticalDuplicates) {
        final File outputDir = IOUtil.createTempDir("MarkDuplicatesTest.", ".tmp");
        outputDir.deleteOnExit();
        final File outputSam = new File(outputDir, "markDuplicatesTest.sam");
        outputSam.deleteOnExit();
        final File metricsFile = new File(outputDir, "markDuplicatesTest.duplicate_metrics");
        metricsFile.deleteOnExit();
        // Run MarkDuplicates, merging the 3 input files, and either enabling or suppressing PG header
        // record creation according to suppressPg.
        final MarkDuplicates markDuplicates = new MarkDuplicates();
        markDuplicates.INPUT = CollectionUtil.makeList(sam);
        markDuplicates.OUTPUT = outputSam;
        markDuplicates.METRICS_FILE = metricsFile;
        markDuplicates.TMP_DIR = CollectionUtil.makeList(outputDir);
        // Needed to suppress calling CommandLineProgram.getVersion(), which doesn't work for code not in a jar
        markDuplicates.PROGRAM_RECORD_ID = null;
        Assert.assertEquals(markDuplicates.doWork(), 0);
        Assert.assertEquals(markDuplicates.numOpticalDuplicates(), expectedNumOpticalDuplicates);
    }

    @DataProvider(name="testOpticalDuplicateDetectionDataProvider")
    public Object[][] testOpticalDuplicateDetectionDataProvider() {
        return new Object[][] {
            {new File(TEST_DATA_DIR, "optical_dupes.sam"), 1L},
            {new File(TEST_DATA_DIR, "optical_dupes_casava.sam"), 1L},
        };
    }
}
