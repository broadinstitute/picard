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
import org.testng.Assert;
import picard.cmdline.CommandLineProgram;

import java.util.List;

/**
 * This class is an extension of AbstractMarkDuplicatesCommandLineProgramTester used to test
 * AbstractMarkDuplicatesCommandLineProgram's with SAM files generated on the fly.  This performs the underlying tests
 * defined by classes such as AbstractMarkDuplicatesCommandLineProgramTest.
 * @author fleharty
 */

public class UmiAwareMarkDuplicatesWithMateCigarTester extends AbstractMarkDuplicatesCommandLineProgramTester {
    private int readNameCounter = 0;
    private List<String> expectedInferredUmis;

    public void addMatePairWithUmi(final String umi, final String inferredUmi, final boolean isDuplicate1, final boolean isDuplicate2) {

        String readName = "READ" + readNameCounter++;
        String cigar1 = null;
        String cigar2 = null;
        boolean strand1 = false;
        boolean strand2 = true;

        int referenceSequenceIndex1 = 0;
        int referenceSequenceIndex2 = 0;
        int alignmentStart1 = 20;
        int alignmentStart2 = 20;

        boolean record1Unmapped = false;
        boolean record2Unmapped = false;

        boolean firstOnly = false;
        boolean record1NonPrimary = false;
        boolean record2NonPrimary = false;

        int defaultQuality = 10;

        addMatePairWithUmi(readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2, record1Unmapped,
                record2Unmapped, isDuplicate1, isDuplicate2, cigar1, cigar2, strand1, strand2, firstOnly, record1NonPrimary, record2NonPrimary,
                defaultQuality, umi, inferredUmi);

    }

    public void addMatePairWithUmi(final String readName,
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
                            final String inferredUmi) {
        final List<SAMRecord> samRecordList = samRecordSetBuilder.addPair(readName, referenceSequenceIndex1, referenceSequenceIndex2, alignmentStart1, alignmentStart2,
                record1Unmapped, record2Unmapped, cigar1, cigar2, strand1, strand2, record1NonPrimary, record2NonPrimary, defaultQuality);

        final SAMRecord record1 = samRecordList.get(0);
        final SAMRecord record2 = samRecordList.get(1);

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

        if(umi != null) {
            record1.setAttribute("RX", umi);
            record2.setAttribute("RX", umi);
        }
        if(inferredUmi != null) {
            record1.setAttribute("RE", inferredUmi);
            record2.setAttribute("RE", inferredUmi);
        }
    }

    public void setExpectedInferredUmis(final List<String> expectedInferredUmis) {
        this.expectedInferredUmis = expectedInferredUmis;
    }

    @Override
    public void test() {
        final SamReader reader = SamReaderFactory.makeDefault().open(getOutput());
        for (final SAMRecord record : reader) {
            // If there are expected inferred UMIs, check to make sure they match
            if (expectedInferredUmis != null) {
                Assert.assertEquals(record.getAttribute("RI"), record.getAttribute("RE"));
            }
        }
        // Also do tests from AbstractMarkDuplicatesCommandLineProgramTester
        super.test();
    }

    @Override
    protected CommandLineProgram getProgram() { return new UmiAwareMarkDuplicatesWithMateCigar(); }
}
