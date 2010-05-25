/*
 * The MIT License
 *
 * Copyright (c) 2010 The Broad Institute
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
package net.sf.picard.util;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.TextCigarCodec;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Basic positive tests for testing cigar string clipping
 *
 * @author Martha Borkan  mborkan@broadinstitute.org
 */
public class CigarUtilTest {

    TextCigarCodec codec = TextCigarCodec.getSingleton();

   @Test(dataProvider="clipData")
    public void basicTest(final String testName, final int start, final String inputCigar, final boolean negativeStrand,
                          final int clipPosition,
                          final String expectedCigar, final int expectedAdjustedStart) throws IOException {
      List<CigarElement> cigar =  codec.decode(inputCigar).getCigarElements();
      if (negativeStrand){
          List<CigarElement> copiedList = new ArrayList<CigarElement>(cigar);
          Collections.reverse(copiedList);
          cigar = copiedList;
      }
      List<CigarElement> result = CigarUtil.softClipEndOfRead( clipPosition, cigar);
      if (negativeStrand){
          Collections.reverse(result);

          int oldLength = (new Cigar(cigar)).getReferenceLength();
          int newLength = (new Cigar(result)).getReferenceLength();
          int sizeChange = oldLength - newLength;
          //Assert.assertEquals(sizeChange, numClippedBases + adjustment, testName + " sizeChange == numClippedBases");
          Assert.assertEquals(start + sizeChange, expectedAdjustedStart, sizeChange + " " +  testName);
          Assert.assertTrue(sizeChange >= 0, "sizeChange >= 0. " + sizeChange);
      }
       Assert.assertEquals (codec.encode(new Cigar(result)), expectedCigar, testName);
    }

    @DataProvider(name = "clipData")
    private Object[][] getCigarClippingTestData() {
        // numClippedBases = (readLength - clipPosition) +1
        return new Object[][]{
            {"Test 1:simple + strand", 100, "50M", false, 43, "42M8S", 100},
            {"Test 1s:+ strand already clipped", 100, "42M8S", false, 43, "42M8S", 100},
            {"Test 2:simple - strand", 100, "50M", true, 41, "10S40M", 110},
            {"Test 3:boundary + strand", 100, "42M3D8M", false, 43, "42M8S", 100},
            {"Test 3s:boundary + strand", 100, "42M3D8S", false, 43, "42M8S", 100},
            {"Test 3x:stutter + strand", 100, "42M2D1D8M", false, 43, "42M8S", 100},
            {"Test 3y:stutter + strand", 100, "42M1D2D8M", false, 43, "42M8S", 100},
            {"Test 3a:boundary + strand", 100, "42M1D8M", false, 43, "42M8S", 100},
            {"Test 4:boundary - strand", 98, "10M2D40M", true, 41, "10S40M", 110},
            {"Test 5:deletion + strand", 100, "44M1D6M", false, 43, "42M8S", 110},
            {"Test 6:deletion - strand", 98, "6M2D44M", true, 41, "10S40M", 110},

            {"Test 7:insertion + strand", 100, "42M3I5M", false, 43, "42M8S", 100},
            {"Test 8:insertion - strand", 102, "8M2I40M", true, 41, "10S40M", 110},
            {"Test 9:insertion within + strand", 100, "44M2I4M", false, 43, "42M8S", 100},
            {"Test 9x:insertion stutter within + strand", 100, "44M2I2I2M", false, 43, "42M8S", 100},
            {"Test 10:insertion within - strand", 100, "3M3I44M", true, 41, "10S40M", 107},
            {"Test 11:insertion straddling + strand", 100, "40M4I6M", false, 43, "40M10S", 100},
            {"Test 11s:insertion straddling + strand", 100, "40M4I6S", false, 43, "40M10S", 100},
            {"Test 11a:insertion straddling + strand", 100, "40M2I8M", false, 43, "40M10S", 100},
            {"Test 12:insertion straddling - strand", 104, "4M4I42M", true, 41, "10S40M", 110},
            {"Test 12a:insertion straddling - strand", 102, "8M2I40M", true, 41, "10S40M", 110},

            {"Test 13:deletion before clip + strand", 100, "10M5D35M", false, 38, "10M5D27M8S", 100},
            {"Test 14:deletion before clip - strand", 100, "35M5D10M", true, 36, "10S25M5D10M", 110},
            {"Test 15:insertion before clip + strand", 100, "10M5I35M", false, 43, "10M5I27M8S", 100},
            {"Test 16:insertion before clip - strand", 100, "16M5I29M", true, 41, "10S6M5I29M", 110},

            {"Test 17:second, earlier clip", 100, "48M2S", false, 43, "42M8S", 100},
            {"Test 17s:second, earlier clip", 100, "2S48M", true, 43, "8S42M", 106},
            {"Test 18:second, later clip", 100, "42M8S", false, 48, "42M8S", 100},
            {"Test 18s:second, later clip", 100, "8S42M", true, 48, "8S42M", 100},
        };
    }
}
