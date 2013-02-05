/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.vcf;

import org.broadinstitute.variant.VariantBaseTest;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.testng.Assert;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: aaron
 * Date: Jun 30, 2010
 * Time: 3:32:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class VCFStandardHeaderLinesUnitTest extends VariantBaseTest {
    @DataProvider(name = "getStandardLines")
    public Object[][] makeGetStandardLines() {
        List<Object[]> tests = new ArrayList<Object[]>();

        // info
        tests.add(new Object[]{"AC", "info", true});
        tests.add(new Object[]{"AN", "info", true});
        tests.add(new Object[]{"AF", "info", true});
        tests.add(new Object[]{"DP", "info", true});
        tests.add(new Object[]{"DB", "info", true});
        tests.add(new Object[]{"END", "info", true});

        // format
        tests.add(new Object[]{"GT", "format", true});
        tests.add(new Object[]{"GQ", "format", true});
        tests.add(new Object[]{"DP", "format", true});
        tests.add(new Object[]{"AD", "format", true});
        tests.add(new Object[]{"PL", "format", true});

        tests.add(new Object[]{"NOT_STANDARD", "info", false});
        tests.add(new Object[]{"NOT_STANDARD", "format", false});

        return tests.toArray(new Object[][]{});
    }


    @Test(dataProvider = "getStandardLines")
    public void getStandardLines(final String key, final String type, final boolean expectedToBeStandard) {
        VCFCompoundHeaderLine line = null;
        if ( type.equals("info") )
            line = VCFStandardHeaderLines.getInfoLine(key, false);
        else if ( type.equals("format") )
            line = VCFStandardHeaderLines.getFormatLine(key, false);
        else
            throw new IllegalArgumentException("Unexpected type in getStandardLines " + type);

        if ( expectedToBeStandard ) {
            Assert.assertNotNull(line);
            Assert.assertEquals(line.getID(), key);
        } else
            Assert.assertNull(line);
    }

    private class RepairHeaderTest {
        final VCFCompoundHeaderLine original, expectedResult;

        private RepairHeaderTest(final VCFCompoundHeaderLine original) {
            this(original, original);
        }

        private RepairHeaderTest(final VCFCompoundHeaderLine original, final VCFCompoundHeaderLine expectedResult) {
            this.original = original;
            this.expectedResult = expectedResult;
        }

        public String toString() {
            return "RepairHeaderTest: Original: " + original.toStringEncoding() + "  Expected: " + expectedResult.toStringEncoding();
        }
    }

    @DataProvider(name = "RepairHeaderTest")
    public Object[][] makeRepairHeaderTest() {
        final VCFInfoHeaderLine standardAC = VCFStandardHeaderLines.getInfoLine("AC");
        final VCFInfoHeaderLine goodAC = new VCFInfoHeaderLine("AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "x");

        final VCFFormatHeaderLine standardGT = VCFStandardHeaderLines.getFormatLine("GT");
        final VCFFormatHeaderLine goodGT = new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "x");

        List<Object[]> tests = new ArrayList<Object[]>();

        tests.add(new Object[]{new RepairHeaderTest( standardGT, standardGT)});
        tests.add(new Object[]{new RepairHeaderTest( goodGT, goodGT )});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("GT", 2, VCFHeaderLineType.String, "x"), standardGT)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.Integer, "x"), standardGT)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.Float, "x"), standardGT)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("GT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Float, "x"), standardGT)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("GT", VCFHeaderLineCount.G, VCFHeaderLineType.String, "x"), standardGT)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("GT", VCFHeaderLineCount.A, VCFHeaderLineType.String, "x"), standardGT)});

        tests.add(new Object[]{new RepairHeaderTest( standardAC, standardAC)});
        tests.add(new Object[]{new RepairHeaderTest( goodAC, goodAC )});
        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("AC", 1, VCFHeaderLineType.Integer, "x"), standardAC)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("AC", VCFHeaderLineCount.G, VCFHeaderLineType.Integer, "x"), standardAC)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("AC", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer, "x"), standardAC)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("AC", 1, VCFHeaderLineType.Float, "x"), standardAC)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("AC", 1, VCFHeaderLineType.String, "x"), standardAC)});
        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("AC", 0, VCFHeaderLineType.Flag, "x"), standardAC)});

        tests.add(new Object[]{new RepairHeaderTest( new VCFInfoHeaderLine("NON_STANDARD_INFO", 1, VCFHeaderLineType.String, "x"))});
        tests.add(new Object[]{new RepairHeaderTest( new VCFFormatHeaderLine("NON_STANDARD_FORMAT", 1, VCFHeaderLineType.String, "x"))});

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "RepairHeaderTest")
    public void testRepairHeaderTest(RepairHeaderTest cfg) {
        final VCFHeader toRepair = new VCFHeader(Collections.singleton((VCFHeaderLine)cfg.original));
        final VCFHeader repaired = VCFStandardHeaderLines.repairStandardHeaderLines(toRepair);

        VCFCompoundHeaderLine repairedLine = (VCFCompoundHeaderLine)repaired.getFormatHeaderLine(cfg.original.getID());
        if ( repairedLine == null ) repairedLine = (VCFCompoundHeaderLine)repaired.getInfoHeaderLine(cfg.original.getID());

        Assert.assertNotNull(repairedLine, "Repaired header didn't contain the expected line");
        Assert.assertEquals(repairedLine.getID(), cfg.expectedResult.getID());
        Assert.assertEquals(repairedLine.getType(), cfg.expectedResult.getType());
        Assert.assertEquals(repairedLine.getCountType(), cfg.expectedResult.getCountType());
        if ( repairedLine.getCountType() == VCFHeaderLineCount.INTEGER )
            Assert.assertEquals(repairedLine.getCount(), cfg.expectedResult.getCount());
    }
}
