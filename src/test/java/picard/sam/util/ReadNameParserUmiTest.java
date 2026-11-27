/*
 * The MIT License
 *
 * Copyright (c) 2026 The Broad Institute
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
package picard.sam.util;

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

/**
 * Checks that the regex-free UMI parser ({@link ReadNameParser#ILLUMINA_UMI_READ_NAME_REGEX}) returns exactly what the
 * equivalent regex does, on names that match, names that do not, and edge cases around the UMI field.
 */
public class ReadNameParserUmiTest {

    /** The regex that the regex-free parser is specified to be equivalent to. */
    private static final String EQUIVALENT_REGEX = "(?:.*:)?([0-9]+):([0-9]+):([0-9]+)(:[ACGTN+-]+)?$";

    @DataProvider(name = "names")
    public Object[][] names() {
        return new Object[][]{
                // no UMI, default Illumina shapes
                {"FC:1:1203:2886:82292"},
                {"RUNID:7:1203:2884:16834"},
                {"A:B:1:2:3"},
                {"1:2:3"},
                {"0:0:0"},
                // UMI present
                {"FC:1:1101:12345:67890:ACGT"},
                {"FC:1:1101:12345:67890:AC+GT"},
                {"FC:1:1101:12345:67890:ACG-T"},
                {"FC:1:1101:12345:67890:N"},
                {"1:2:3:ACGT"},
                {"A:1:2:3:TTGCA"},
                // malformed UMI or coordinate fields
                {"FC:1:1101:12345:67890:"},
                {"FC:1:1101:12345:67890:acgt"},
                {"FC:1:1101:12345:67890:ACGT:"},
                {"FC:1:1101:12345:67890:ACGT:TTGC"},
                {"FC:1:1101:12345:67890:ACXT"},
                {"FC:1:1101:12345:67890:ACGT:1"},
                {"FC:1:1101:x12345:67890:ACGT"},
                {"FC:1:1101:12345:67890x:ACGT"},
                {"FC:1:+1101:12345:67890"},
                {"FC:1:1101:12345:-67890"},
                {"FC:1:1101:12345:"},
                {"FC:1:1101:12345:67890:"},
                {"FC:1:1101:12345"},
                {"1101:12345:ACGT"},
                {"FC:1:1101:12345:ACGT"},
                {"ACGT"},
                {":"},
                {"::"},
                {":::"},
                {""},
                // out of int range
                {"FC:1:99999999999:1:2:ACGT"},
                {"FC:1:1:99999999999:2"},
                {"FC:1:1:2:99999999999"},
                // digits in the UMI-looking last field are not a UMI
                {"FC:1:1101:12345:67890:1"},
                {"FC:1:1101:12345:67890:0ACGT"},
        };
    }

    @Test(dataProvider = "names")
    public void testUmiParserMatchesRegex(final String readName) {
        final ReadNameParser regexParser = new ReadNameParser(EQUIVALENT_REGEX);
        final ReadNameParser umiParser = new ReadNameParser(ReadNameParser.ILLUMINA_UMI_READ_NAME_REGEX);

        final PhysicalLocationInt fromRegex = new PhysicalLocationInt();
        final PhysicalLocationInt fromUmi = new PhysicalLocationInt();

        final boolean regexMatched = regexParser.addLocationInformation(readName, fromRegex);
        final boolean umiMatched = umiParser.addLocationInformation(readName, fromUmi);

        Assert.assertEquals(umiMatched, regexMatched, "match result differs for " + readName);
        if (regexMatched) {
            Assert.assertEquals(fromUmi.getTile(), fromRegex.getTile(), "tile differs for " + readName);
            Assert.assertEquals(fromUmi.getX(), fromRegex.getX(), "x differs for " + readName);
            Assert.assertEquals(fromUmi.getY(), fromRegex.getY(), "y differs for " + readName);
        }
    }

    @Test
    public void testUmiParserOnlyEnabledByItsOwnConstant() {
        // The default regex path must not change: a UMI suffix still fails to parse there.
        final ReadNameParser defaultParser = new ReadNameParser(ReadNameParser.DEFAULT_READ_NAME_REGEX);
        Assert.assertFalse(defaultParser.addLocationInformation("FC:1:1101:12345:67890:ACGT",
                new PhysicalLocationInt()));
    }
}
