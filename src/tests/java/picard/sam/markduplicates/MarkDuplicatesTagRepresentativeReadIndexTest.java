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

import org.testng.annotations.Test;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesWithMateCigarTester (see getTester).
 * @author nhomer@broadinstitute.org
 */
public class MarkDuplicatesTagRepresentativeReadIndexTest extends AbstractMarkDuplicatesCommandLineProgramTest {
    protected MarkDuplicatesTagRepresentativeeadIndexTester getTester() {
        return new MarkDuplicatesTagRepresentativeeadIndexTester();
    }

    @Test
    public void testRepresentativeReadTag() {
        final MarkDuplicatesTagRepresentativeeadIndexTester tester = getTester();
        tester.testRepresentativeReads = true;
        tester.setExpectedOpticalDuplicate(1);
        String representativeReadName = "RUNID:1:1:16020:13352";
        Integer representativeReadIndexInFileForward = 1;
        tester.addMatePair(representativeReadName, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(representativeReadIndexInFileForward, representativeReadIndexInFileForward);
        tester.expectedRepresentativeIndexMap.put(3, representativeReadIndexInFileForward);
        tester.expectedSetSizeMap.put(representativeReadName,2);
        tester.addMatePair("RUNID:1:1:15993:13361", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(0, representativeReadIndexInFileForward);
        tester.expectedRepresentativeIndexMap.put(2, representativeReadIndexInFileForward);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13361",2);
        tester.runTest();
    }

    @Test
    public void testMultiRepresentativeReadTags() {
        final MarkDuplicatesTagRepresentativeeadIndexTester tester = getTester();
        tester.testRepresentativeReads = true;
        tester.setExpectedOpticalDuplicate(3);
        // Duplicate set: size 2 - all optical
        String representativeReadName1 = "RUNID:1:1:16020:13352";
        Integer representativeReadIndexInFileForward = 1;
        tester.addMatePair(representativeReadName1, 1, 485253, 485253, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(representativeReadIndexInFileForward, representativeReadIndexInFileForward);
        tester.expectedRepresentativeIndexMap.put(3, representativeReadIndexInFileForward);
        tester.expectedSetSizeMap.put(representativeReadName1,2);
        tester.addMatePair("RUNID:1:1:15993:13361", 1, 485253, 485253, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(0, representativeReadIndexInFileForward);
        tester.expectedRepresentativeIndexMap.put(2, representativeReadIndexInFileForward);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13361",2);

        // Duplicate set: size 3 - all optical
        String representativeReadName2 = "RUNID:1:1:15993:13360";
        Integer representativeReadIndexInFileForward2 = 6;
        tester.addMatePair(representativeReadName2, 1, 485299, 485299, false, false, false, false, "45M", "45M", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(representativeReadIndexInFileForward2, representativeReadIndexInFileForward2);
        tester.expectedRepresentativeIndexMap.put(9, representativeReadIndexInFileForward2);
        tester.expectedSetSizeMap.put(representativeReadName2,3);
        tester.addMatePair("RUNID:1:1:15993:13365", 1, 485299, 485299, false, false, true, true, "44M1S", "44M1S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(7, representativeReadIndexInFileForward2);
        tester.expectedRepresentativeIndexMap.put(10, representativeReadIndexInFileForward2);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13365",3);
        tester.addMatePair("RUNID:1:1:15993:13370", 1, 485299, 485299, false, false, true, true, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(8, representativeReadIndexInFileForward2);
        tester.expectedRepresentativeIndexMap.put(11, representativeReadIndexInFileForward2);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13370",3);

        // Add non-duplicate read
        tester.addMatePair("RUNID:1:1:15993:13375", 1, 485255, 485255, false, false, false, false, "43M2S", "43M2S", false, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.expectedRepresentativeIndexMap.put(4, null);
        tester.expectedRepresentativeIndexMap.put(5, null);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13375",null);

        tester.runTest();
    }

}
