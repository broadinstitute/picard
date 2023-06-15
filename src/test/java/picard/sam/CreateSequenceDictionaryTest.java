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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.PicardException;
import picard.cmdline.CommandLineProgramTest;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * @author alecw@broadinstitute.org
 */
public class CreateSequenceDictionaryTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard");
    private static final File BASIC_FASTA = new File(TEST_DATA_DIR + "/sam", "basic.fasta");
    private static final File EQUIVALENCE_TEST_FASTA = new File(TEST_DATA_DIR + "/reference", "test.fasta");
    private static final File DUPLICATE_FASTA = new File(TEST_DATA_DIR + "/sam", "duplicate_sequence_names.fasta");

    private static final File SPECIAL_FASTA = new File(TEST_DATA_DIR + "/sam", "special.fasta");
    private static final File SPECIAL_ALT = new File(TEST_DATA_DIR + "/sam", "special.alt");

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
    public void testDefaultOutputFile()  {
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

        final List<String> currentDict = new BufferedReader(new FileReader(outputDict))
                .lines()
                //remove info about location fasta file
                .map(s -> s.replaceAll("UR:.*", ""))
                .collect(Collectors.toList());

        final List<String> expectedDict = new BufferedReader(new FileReader(TEST_DATA_DIR + "/reference/csd_dict.dict"))
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
        final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(outputDict.toPath());
        Assert.assertNotNull(dict, "dictionary is null");
       
        // check chr1

        final SAMSequenceRecord ssr1 = dict.getSequence("chr1");
        Assert.assertNotNull(ssr1, "chr1 missing in dictionary");
        final String an1 = ssr1.getAttribute("AN");
        Assert.assertNotNull(ssr1, "AN Missing");
        Set<String> anSet = new HashSet<>(Arrays.asList(an1.split("[,]")));

        Assert.assertTrue(anSet.contains("1"));
        Assert.assertTrue(anSet.contains("01"));
        Assert.assertTrue(anSet.contains("k1"));
        Assert.assertFalse(anSet.contains("M"));
        
        // check chr2
        SAMSequenceRecord ssr2 = dict.getSequence("chr2");
        Assert.assertNotNull(ssr2, "chr2 missing in dictionary");
        final String an2 = ssr2.getAttribute("AN");
        Assert.assertNull(an2, "AN Present");
        
        // check chrM
        final SAMSequenceRecord ssrM = dict.getSequence("chrM");
        Assert.assertNull(ssrM, "chrM present in dictionary");
    }

    @Test
    public void testPunctuationInReferenceNames() {
        final File expectedDict = new File(TEST_DATA_DIR + "/sam", "special.dict");
        expectedDict.deleteOnExit();
        Assert.assertFalse(expectedDict.exists());

        final String[] argv = {
                "REFERENCE=" + SPECIAL_FASTA,
                "AN=" + SPECIAL_ALT,
                "TRUNCATE_NAMES_AT_WHITESPACE=false"
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);
        Assert.assertTrue(expectedDict.exists());

        final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(expectedDict);
        Assert.assertNotNull(dict, "dictionary is null");

        // check chr1.01
        {
            final SAMSequenceRecord ssr = dict.getSequence("chr1.01");
            Assert.assertNotNull(ssr, "chr1 missing in dictionary");
            final String an = ssr.getAttribute("AN");
            Assert.assertNotNull(an, "AN Missing");
            final Set<String> anSet = new HashSet<>(Arrays.asList(an.split("[,]")));
            Assert.assertTrue(anSet.contains("1.01"));
        }

        // check chr2/hello
        {
            final SAMSequenceRecord ssr = dict.getSequence("chr2/hello");
            Assert.assertNotNull(ssr, "chr2/hello missing in dictionary");
            final String an = ssr.getAttribute("AN");
            Assert.assertNotNull(an, "AN Missing");
            final Set<String> anSet = new HashSet<>(Arrays.asList(an.split("[,]")));
            Assert.assertTrue(anSet.contains("2/hello"));
        }
    }

    // Test that dictionary MD5s are invariant across dos/unix line endings
    @Test
    public void testSequenceMD5sRespectLineEndings() throws IOException {
        // NOTE: the .fasta test files here that use windows cr/lf line endings have to be listed in
        // the repo .gitattributes file in order to prevent git from converting them to use local
        // platform line endings on checkout
        final Path picardDictPath = makeTempDictionary(
                Paths.get(TEST_DATA_DIR + "/reference/test.fasta"),
                "CreateSequenceDictionaryTest_unix.");

        final Path picardDosDictPath = makeTempDictionary(
                Paths.get(TEST_DATA_DIR + "/reference/test_dos_line_endings.fasta"),
                "CreateSequenceDictionaryTest_dos.");

        final Path picardDosGZDictPath = makeTempDictionary(
                Paths.get(TEST_DATA_DIR + "/reference/test_dos_line_endings.fasta.gz"),
                "CreateSequenceDictionaryTest_dos_gz.");
        final SAMSequenceDictionary unixDict = SAMSequenceDictionaryExtractor.extractDictionary(picardDictPath);
        final SAMSequenceDictionary dosDict = SAMSequenceDictionaryExtractor.extractDictionary(picardDosDictPath);
        final SAMSequenceDictionary dosGZDict = SAMSequenceDictionaryExtractor.extractDictionary(picardDosGZDictPath);

        // test files generated with samtools 1.17

        // converted by samtools from test.fasta
        final SAMSequenceDictionary samtoolsDict = SAMSequenceDictionaryExtractor.extractDictionary(
                Paths.get(TEST_DATA_DIR + "/reference/samtools_test.dict"));
        // converted by samtools from test_windows_line_endings.fasta
        final SAMSequenceDictionary samtoolsDosDict = SAMSequenceDictionaryExtractor.extractDictionary(
                Paths.get(TEST_DATA_DIR + "/reference/samtools_test_dos_line_endings.dict"));
        // from test_windows_line_endings.fasta.gz
        final SAMSequenceDictionary samtoolsDosGZDict = SAMSequenceDictionaryExtractor.extractDictionary(
                Paths.get(TEST_DATA_DIR + "/reference/samtools_test_dos_line_endings_gz.dict"));

        final List<SAMSequenceRecord> picardRecords = unixDict.getSequences();
        final List<SAMSequenceRecord> picardDosRecords = dosDict.getSequences();
        final List<SAMSequenceRecord> picardDosGZRecords = dosGZDict.getSequences();
        final List<SAMSequenceRecord> samtoolsRecords = samtoolsDict.getSequences();
        final List<SAMSequenceRecord> samtoolsDosRecords = samtoolsDosDict.getSequences();
        final List<SAMSequenceRecord> samtoolsDosGZRecords = samtoolsDosGZDict.getSequences();

        Assert.assertEquals(picardRecords.size(), picardDosRecords.size());
        Assert.assertEquals(picardRecords.size(), picardDosGZRecords.size());
        Assert.assertEquals(picardRecords.size(), samtoolsRecords.size());
        Assert.assertEquals(picardRecords.size(), samtoolsDosRecords.size());
        Assert.assertEquals(picardRecords.size(), samtoolsDosGZRecords.size());

        for (int i = 0; i < picardRecords.size(); i++) {
            Assert.assertEquals(picardRecords.get(i).getMd5(), picardDosRecords.get(i).getMd5());
            Assert.assertEquals(picardRecords.get(i).getMd5(), picardDosGZRecords.get(i).getMd5());
            Assert.assertEquals(picardRecords.get(i).getMd5(), samtoolsRecords.get(i).getMd5());
            Assert.assertEquals(picardRecords.get(i).getMd5(), samtoolsDosRecords.get(i).getMd5());
            Assert.assertEquals(picardRecords.get(i).getMd5(), samtoolsDosGZRecords.get(i).getMd5());
        }
    }

    final Path makeTempDictionary(final Path inputFasta, final String dictNamePrefix) throws IOException {
        final File tempDict = File.createTempFile(dictNamePrefix, ".dict");
        tempDict.delete();
        tempDict.deleteOnExit();
        final String[] argv = {
                "REFERENCE=" + inputFasta.toAbsolutePath(),
                "OUTPUT=" + tempDict
        };
        Assert.assertEquals(runPicardCommandLine(argv), 0);
        return tempDict.toPath();
    }
}
