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

import htsjdk.samtools.SAMFileHeader;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesWithMateCigarTester (see getTester).
 * @author hogstrom@broadinstitute.org
 */
public class MarkDuplicatesTagRepresentativeReadIndexTest extends AbstractMarkDuplicatesCommandLineProgramTest {
    protected MarkDuplicatesTagRepresentativeReadIndexTester getTester() {
        return new MarkDuplicatesTagRepresentativeReadIndexTester();
    }

    protected MarkDuplicatesTagRepresentativeReadIndexTester getTester(final SAMFileHeader.SortOrder sortOrder) {
        return new MarkDuplicatesTagRepresentativeReadIndexTester(sortOrder);
    }

    // This tests the creation of a single duplicate pair. The expected size of the duplicate set is 2
    // and the expected representative read is the first read defined in the text. The test fails if
    // the 'DI' and 'DS' tags do not match the expectation.
    @DataProvider(name = "sortOrderDataProvider")
    public Object[][] sortOrderDataProvider() {
        return new Object[][]{
                {SAMFileHeader.SortOrder.coordinate},
                {SAMFileHeader.SortOrder.queryname}
        };
    }

    @Test(dataProvider = "sortOrderDataProvider")
    public void testRepresentativeReadTag(final SAMFileHeader.SortOrder sortOrder) {
        final MarkDuplicatesTagRepresentativeReadIndexTester tester = getTester(sortOrder);
        tester.getSamRecordSetBuilder().setReadLength(45);
        tester.testRepresentativeReads = true;
        tester.setExpectedOpticalDuplicate(1);
        final String representativeReadName = "RUNID:1:1:16020:13352";
        tester.addMatePair(representativeReadName, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(representativeReadName, representativeReadName);
        tester.expectedSetSizeMap.put(representativeReadName,2);
        final String duplicateReadName = "RUNID:1:1:15993:13361";
        tester.addMatePair(duplicateReadName, 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(duplicateReadName, representativeReadName);
        tester.expectedSetSizeMap.put(duplicateReadName,2);
        tester.runTest();
    }

    // This tests the creation of a two duplicate sets of size 2 and 3 respectively. The test fails if
    // the 'DI' and 'DS' tags do not match the expectation.

    @Test(dataProvider = "sortOrderDataProvider")
    public void testMultiRepresentativeReadTags(final SAMFileHeader.SortOrder sortOrder) {
        final MarkDuplicatesTagRepresentativeReadIndexTester tester = getTester(sortOrder);
        tester.getSamRecordSetBuilder().setReadLength(45);
        tester.testRepresentativeReads = true;
        tester.setExpectedOpticalDuplicate(3);
        // Duplicate set: size 2 - all optical
        final String representativeReadName1 = "RUNID:1:1:16020:13352";
        tester.addMatePair(representativeReadName1, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(representativeReadName1, representativeReadName1);
        tester.expectedSetSizeMap.put(representativeReadName1,2);
        final String duplicateSet1DuplicateReadName = "RUNID:1:1:15993:13361";
        tester.addMatePair(duplicateSet1DuplicateReadName, 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(duplicateSet1DuplicateReadName, representativeReadName1);
        tester.expectedSetSizeMap.put(duplicateSet1DuplicateReadName,2);

        // Duplicate set: size 3 - all optical
        final String representativeReadName2 = "RUNID:1:1:15993:13360";
        tester.addMatePair(representativeReadName2, 1, 485299, 485299, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(representativeReadName2, representativeReadName2);
        tester.expectedSetSizeMap.put(representativeReadName2,3);
        final String duplicateSet2DuplicateReadName1 = "RUNID:1:1:15993:13365";
        tester.addMatePair(duplicateSet2DuplicateReadName1, 1, 485299, 485299, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(duplicateSet2DuplicateReadName1, representativeReadName2);
        tester.expectedSetSizeMap.put(duplicateSet2DuplicateReadName1,3);
        final String duplicateSet2DuplicateReadName2 = "RUNID:1:1:15993:13370";
        tester.addMatePair(duplicateSet2DuplicateReadName2, 1, 485299, 485299, false, false, true, true, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(duplicateSet2DuplicateReadName2, representativeReadName2);
        tester.expectedSetSizeMap.put(duplicateSet2DuplicateReadName2,3);

        // Add non-duplicate read
        final String nonDuplicateReadName = "RUNID:1:1:15993:13375";
        tester.addMatePair(nonDuplicateReadName, 1, 485255, 485255, false, false, false, false, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeReadMap.put(nonDuplicateReadName, null);
        tester.expectedSetSizeMap.put(nonDuplicateReadName,null);

        tester.runTest();
    }

}
