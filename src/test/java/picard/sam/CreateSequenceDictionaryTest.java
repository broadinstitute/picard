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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import picard.cmdline.CommandLineProgramTest;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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
    
    
    @Test
    public void testAltNames() throws Exception {
        final File altFile = File.createTempFile("CreateSequenceDictionaryTest", ".alt");
        final PrintWriter pw = new PrintWriter(altFile);
        pw.println("chr1\t1");
        pw.println("chr1\t01");
        pw.println("chr1\tk1");
        pw.println("chrMT\tM");
        pw.flush();
        pw.close();
        altFile.deleteOnExit();
        
        final File outputDict = File.createTempFile("CreateSequenceDictionaryTest.", ".dict");
        outputDict.delete();
        outputDict.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + BASIC_FASTA,
                "AN=" + altFile,
                "OUTPUT=" + outputDict,
                "TRUNCATE_NAMES_AT_WHITESPACE=true"
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);
        final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(outputDict);
        Assert.assertNotNull(dict, "dictionary is null");
       
        // check chr1
        SAMSequenceRecord ssr = dict.getSequence("chr1");
        Assert.assertNotNull(ssr, "chr1 missing in dictionary");
        String an = ssr.getAttribute("AN");
        Assert.assertNotNull(ssr, "AN Missing");
        Set<String> anSet = new HashSet<>(Arrays.asList(an.split("[,]")));
        Assert.assertTrue(anSet.contains("1"));
        Assert.assertTrue(anSet.contains("01"));
        Assert.assertTrue(anSet.contains("k1"));
        Assert.assertFalse(anSet.contains("M"));
        
        // check chr2
        ssr = dict.getSequence("chr2");
        Assert.assertNotNull(ssr, "chr2 missing in dictionary");
        an = ssr.getAttribute("AN");
        Assert.assertNull(an, "AN Present");
        
        // check chrM
        ssr = dict.getSequence("chrM");
        Assert.assertNull(ssr, "chrM present in dictionary");
    }

}
