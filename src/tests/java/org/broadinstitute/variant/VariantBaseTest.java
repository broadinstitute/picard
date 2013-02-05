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

package org.broadinstitute.variant;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.testng.Assert;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Base class for test classes within org.broadinstitute.variant
 */
public class VariantBaseTest {

    public static final String variantTestDataRoot = new File("testdata/variant/").getAbsolutePath() + "/";

    /**
     * Creates a temp file that will be deleted on exit after tests are complete.
     * @param name Prefix of the file.
     * @param extension Extension to concat to the end of the file.
     * @return A file in the temporary directory starting with name, ending with extension, which will be deleted after the program exits.
     */
    public static File createTempFile(String name, String extension) {
        try {
            File file = File.createTempFile(name, extension);
            file.deleteOnExit();
            return file;
        } catch (IOException ex) {
            throw new RuntimeException("Cannot create temp file: " + ex.getMessage(), ex);
        }
    }

    private static final double DEFAULT_FLOAT_TOLERANCE = 1e-1;

    public static final void assertEqualsDoubleSmart(final Object actual, final Double expected) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double)(Double)actual, (double)expected);
    }

    public static final void assertEqualsDoubleSmart(final Object actual, final Double expected, final double tolerance) {
        Assert.assertTrue(actual instanceof Double, "Not a double");
        assertEqualsDoubleSmart((double)(Double)actual, (double)expected, tolerance);
    }

    public static final void assertEqualsDoubleSmart(final double actual, final double expected) {
        assertEqualsDoubleSmart(actual, expected, DEFAULT_FLOAT_TOLERANCE);
    }

    public static final <T> void assertEqualsSet(final Set<T> actual, final Set<T> expected, final String info) {
        final Set<T> actualSet = new HashSet<T>(actual);
        final Set<T> expectedSet = new HashSet<T>(expected);
        Assert.assertTrue(actualSet.equals(expectedSet), info); // note this is necessary due to testng bug for set comps
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance) {
        assertEqualsDoubleSmart(actual, expected, tolerance, null);
    }

    public static void assertEqualsDoubleSmart(final double actual, final double expected, final double tolerance, final String message) {
        if ( Double.isNaN(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isNaN(actual), "expected is nan, actual is not");
        else if ( Double.isInfinite(expected) ) // NaN == NaN => false unfortunately
            Assert.assertTrue(Double.isInfinite(actual), "expected is infinite, actual is not");
        else {
            final double delta = Math.abs(actual - expected);
            final double ratio = Math.abs(actual / expected - 1.0);
            Assert.assertTrue(delta < tolerance || ratio < tolerance, "expected = " + expected + " actual = " + actual
                    + " not within tolerance " + tolerance
                    + (message == null ? "" : "message: " + message));
        }
    }

    public static SAMSequenceDictionary createArtificialSequenceDictionary() {
        final int[] contigLengths = { 249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
                                      141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
                                      81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16569 };
        List<SAMSequenceRecord> contigs = new ArrayList<SAMSequenceRecord>();

        for ( int contig = 1; contig <= 22; contig++ ) {
            contigs.add(new SAMSequenceRecord(Integer.toString(contig), contigLengths[contig - 1]));
        }

        int position = 22;
        for ( String contigName : Arrays.asList("X", "Y", "MT") ) {
            contigs.add(new SAMSequenceRecord(contigName, contigLengths[position]));
            position++;
        }

        return new SAMSequenceDictionary(contigs);
    }

}
