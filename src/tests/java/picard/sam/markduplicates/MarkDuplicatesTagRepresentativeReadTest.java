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

import picard.PicardException;
import org.testng.annotations.Test;

/**
 * This class defines the individual test cases to run. The actual running of the test is done
 * by MarkDuplicatesWithMateCigarTester (see getTester).
 * @author nhomer@broadinstitute.org
 */
public class MarkDuplicatesTagRepresentativeReadTest extends AbstractMarkDuplicatesCommandLineProgramTest {
    protected MarkDuplicatesTagRepresentativeReadTester getTester() {
        return new MarkDuplicatesTagRepresentativeReadTester();
    }

    @Test
    public void testRepresentativeReadTag() {
        final MarkDuplicatesTagRepresentativeReadTester tester = getTester();
        tester.testRepresentativeReads = true;
        tester.setExpectedOpticalDuplicate(2);
        int HIGH_BASE_QUALITY = DEFAULT_BASE_QUALITY+5;
        String representativeReadName = "RUNID:1:1:16020:13352";
        tester.addMatePair("RUNID:1:1:15993:13361", 2, 41212324, 41212310, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.readToRepReadMap.put("RUNID:1:1:15993:13361", representativeReadName);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13361",3);
        System.out.println(tester.readToRepReadMap.get("RUNID:1:1:15993:13361"));
        //the representative read in the set should have the highest base quality
        tester.addMatePair(representativeReadName, 2, 41212324, 41212310, false, false, false, false, "33S35M", "28S40M", true, true, false, false, false, HIGH_BASE_QUALITY);
        tester.readToRepReadMap.put(representativeReadName, representativeReadName);
        tester.expectedSetSizeMap.put(representativeReadName,3);
        tester.addMatePair("RUNID:1:1:15994:13363", 2, 41212324, 41212310, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.readToRepReadMap.put("RUNID:1:1:15994:13363", representativeReadName);
        tester.expectedSetSizeMap.put("RUNID:1:1:15994:13363",3);
        tester.runTest();
    }

    @Test
    public void testMultiRepresentativeReadTags() {
        final MarkDuplicatesTagRepresentativeReadTester tester = getTester();
        tester.testRepresentativeReads = true;
        tester.setExpectedOpticalDuplicate(3);
        int HIGH_BASE_QUALITY = DEFAULT_BASE_QUALITY+5;

        // Duplicate set: size 2
        String representativeReadName1 = "RUNID:1:1:16020:14000";
        tester.addMatePair(representativeReadName1, 2, 31212324, 31212310, false, false, false, false, "33S35M", "28S40M", true, true, false, false, false, HIGH_BASE_QUALITY);
        tester.readToRepReadMap.put(representativeReadName1, representativeReadName1);
        tester.expectedSetSizeMap.put(representativeReadName1,2);
        tester.addMatePair("RUNID:1:1:15993:14001", 2, 31212324, 31212310, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.readToRepReadMap.put("RUNID:1:1:15993:14001", representativeReadName1);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:14001",2);
        //System.out.println(tester.readToRepReadMap.get("RUNID:1:1:15993:13361"));

        // Duplicate set: size 3
        String representativeReadName2 = "RUNID:1:1:16020:13352";
        tester.addMatePair("RUNID:1:1:15993:13361", 2, 41212324, 41212310, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.readToRepReadMap.put("RUNID:1:1:15993:13361", representativeReadName2);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:13361",3);
        System.out.println(tester.readToRepReadMap.get("RUNID:1:1:15993:13361"));
        //the representative read in the set should have the highest base quality
        tester.addMatePair(representativeReadName2, 2, 41212324, 41212310, false, false, false, false, "33S35M", "28S40M", true, true, false, false, false, HIGH_BASE_QUALITY);
        tester.readToRepReadMap.put(representativeReadName2, representativeReadName2);
        tester.expectedSetSizeMap.put(representativeReadName2,3);
        tester.addMatePair("RUNID:1:1:15994:13363", 2, 41212324, 41212310, false, false, true, true, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.readToRepReadMap.put("RUNID:1:1:15994:13363", representativeReadName2);
        tester.expectedSetSizeMap.put("RUNID:1:1:15994:13363",3);

        // Add non-duplicate read
        tester.addMatePair("RUNID:1:1:15993:15000", 2, 21212324, 21212310, false, false, false, false, "33S35M", "28S40M", true, true, false, false, false, DEFAULT_BASE_QUALITY);
        tester.readToRepReadMap.put("RUNID:1:1:15993:15000", null);
        tester.expectedSetSizeMap.put("RUNID:1:1:15993:15000",null);

        tester.runTest();
    }

}
