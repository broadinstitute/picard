/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
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

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author alecw@broadinstitute.org
 */
public class CreateSequenceDictionaryTest extends CommandLineProgramTest {
    public static File TEST_DATA_DIR = new File("testdata/picard");
    public static File BASIC_FASTA = new File(TEST_DATA_DIR + "/sam", "basic.fasta");
    public static File EQUIVALENCE_TEST_FASTA = new File(TEST_DATA_DIR + "/reference", "test.fasta");
    public static File DUPLICATE_FASTA = new File(TEST_DATA_DIR + "/sam", "duplicate_sequence_names.fasta");

    public String getCommandLineProgramName() {
        return CreateSequenceDictionary.class.getSimpleName();
    }

    @Test
    public void testBasic() throws Exception {
        final File outputDict = File.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        outputDict.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + BASIC_FASTA,
                "OUTPUT=" + outputDict,
                "TRUNCATE_NAMES_AT_WHITESPACE=false"
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);
    }

    @Test
    public void testDefaultOutputFile() throws Exception {
        final File expectedDict = new File(TEST_DATA_DIR + "/sam", "basic.dict");
        expectedDict.deleteOnExit();
        Assert.assertFalse(expectedDict.exists());
        final String[] argv = {
                "REFERENCE=" + BASIC_FASTA,
                "TRUNCATE_NAMES_AT_WHITESPACE=false"
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);
        Assert.assertTrue(expectedDict.exists());
    }

    @Test
    public void testForEquivalence() throws Exception {
        final File outputDict = File.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        final String[] argv = {
                "REFERENCE=" + EQUIVALENCE_TEST_FASTA,
                "OUTPUT=" + outputDict,
                "TRUNCATE_NAMES_AT_WHITESPACE=false"
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);

        List<String> currentDict = new BufferedReader(new FileReader(outputDict))
                .lines()
                //remove info about location fasta file
                .map(s -> s.replaceAll("UR:.*", ""))
                .collect(Collectors.toList());

        List<String> expectedDict = new BufferedReader(new FileReader(TEST_DATA_DIR + "/reference/csd_dict.dict"))
                .lines()
                //remove info about location fasta file
                .map(s -> s.replaceAll("UR:.*", ""))
                .collect(Collectors.toList());

        Assert.assertEquals(currentDict, expectedDict);
    }

    /**
     * Should throw an exception because with TRUNCATE_NAMES_AT_WHITESPACE, sequence names are not unique.
     */
    @Test(expectedExceptions = PicardException.class)
    public void testNonUniqueSequenceName() throws Exception {
        final File outputDict = File.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + DUPLICATE_FASTA,
                "OUTPUT=" + outputDict,
                "TRUNCATE_NAMES_AT_WHITESPACE=true"
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);
        Assert.fail("Exception should have been thrown.");
    }
}
