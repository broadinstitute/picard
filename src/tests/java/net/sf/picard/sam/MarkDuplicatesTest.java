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
package net.sf.picard.sam;

import net.sf.picard.io.IoUtil;
import net.sf.samtools.*;
import net.sf.samtools.util.TestUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class MarkDuplicatesTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam/MarkDuplicates");


    /**
     * Test that PG header records are created & chained appropriately (or not created), and that the PG record chains
     * are as expected.  MarkDuplicates is used both to merge and to mark dupes in this case.
     * @param suppressPg If true, do not create PG header record.
     * @param expectedPnVnByReadName For each read, info about the expect chain of PG records.
     */
    @Test(dataProvider = "pgRecordChainingTest")
    public void pgRecordChainingTest(final boolean suppressPg,
                                     final Map<String, List<ExpectedPnAndVn>> expectedPnVnByReadName) {
        final File outputDir = IoUtil.createTempDir("MarkDuplicatesTest.", ".tmp");
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
}
