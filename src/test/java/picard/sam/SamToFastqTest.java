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

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.IOUtil;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Tests for SamToFastq
 */
public class SamToFastqTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/bam2fastq/paired");
    private static final String CLIPPING_TEST_DATA = "ok/clipping_test.sam";

    public String getCommandLineProgramName() {
        return SamToFastq.class.getSimpleName();
    }

    @DataProvider(name = "okFiles")
    public Object[][] okFiles() {
        return new Object[][] {
                {"ok/sorted-pair.sam"}, // 5 sorted pairs (10 records) - mate1, mate2
                {"ok/sorted-pair-no-rg.sam"}, // 5 sorted pairs (10 records) - mate1, mate2, no read group
                {"ok/last-pair-mates-flipped.sam" }, // 4 pairs, :05 mate2, mate1
                {"ok/first-mate-bof-last-mate-eof.sam"}, // :01 mate1, 4 pairs, :01 mate2
        };
    }


    @DataProvider(name = "badFiles")
    public Object[][] badFiles() {
        return new Object[][] {
            {"bad/unpaired-mate.sam"} // mate1 without its mate2
        };
    }

    private void convertFile(final String [] args) {
        Assert.assertEquals(runPicardCommandLine(args), 0);
    }

    @Test(dataProvider = "clippingTests")
    public void testClipping(final String clippingAction, final String bases1_1, final String quals1_1, final String bases1_2, final String quals1_2,
                             final String bases2_1, final String quals2_1, final String bases2_2, final String quals2_2, final String testName) throws IOException {
        final File samFile = new File(TEST_DATA_DIR, CLIPPING_TEST_DATA) ;
        final File f1 = File.createTempFile("clippingtest1", "fastq");
        final File f2 = File.createTempFile("clippingtest2", "fastq");
        f1.deleteOnExit();
        f2.deleteOnExit();

        if (clippingAction != null) {
            convertFile(new String[]{
                "INPUT="            + samFile.getAbsolutePath(),
                "FASTQ="            + f1.getAbsolutePath(),
                "SECOND_END_FASTQ=" + f2.getAbsolutePath(),
                "CLIPPING_ACTION="  + clippingAction,
                "CLIPPING_ATTRIBUTE=" + "XT"
            });
        } else {
            convertFile(new String[]{
                "INPUT="            + samFile.getAbsolutePath(),
                "FASTQ="            + f1.getAbsolutePath(),
                "SECOND_END_FASTQ=" + f2.getAbsolutePath(),
            });
        }

        Iterator<FastqRecord> it = new FastqReader(f1).iterator();
        FastqRecord first = it.next();
        Assert.assertEquals(first.getReadString(), bases1_1, testName);
        Assert.assertEquals(first.getBaseQualityString(), quals1_1, testName);
        FastqRecord second = it.next();
        Assert.assertEquals(second.getReadString(), bases1_2, testName);
        Assert.assertEquals(second.getBaseQualityString(), quals1_2, testName);
        it = new FastqReader(f2).iterator();
        first = it.next();
        Assert.assertEquals(first.getReadString(), bases2_1, testName);
        Assert.assertEquals(first.getBaseQualityString(), quals2_1, testName);
        second = it.next();
        Assert.assertEquals(second.getReadString(), bases2_2, testName);
        Assert.assertEquals(second.getBaseQualityString(), quals2_2, testName);
    }

    @DataProvider(name = "clippingTests")
    public Object[][] clippingTests() {
        return new Object[][] {
            {null, "AAAAAAAAAA", "1111111111", "AAAAAAAAAA", "1111111111", "CCCCCCCCCC", "2222222222", "GGGGGGGGGG", "2222222222", "No clipping test"},
            {"X",  "AAAAAAA",    "1111111",    "AAAAAA",     "111111",     "CCCCCCCC",   "22222222",   "GGGGGG",     "222222",     "Cut clipped bases test"},
            {"N",  "AAAAAAANNN", "1111111111", "AAAAAANNNN", "1111111111", "CCCCCCCCNN", "2222222222", "GGGGGGNNNN", "2222222222", "Mask clipped bases test"},
            {"2",  "AAAAAAAAAA", "1111111###", "AAAAAAAAAA", "111111####", "CCCCCCCCCC", "22222222##", "GGGGGGGGGG", "222222####", "Change clipped qualities test"}
        };
    }

    @Test(dataProvider = "okFiles")
    public void testOkFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");

        convertFile(new String[]{
              "INPUT=" + samFile.getAbsolutePath(),
              "FASTQ=" + pair1File.getAbsolutePath(),
              "SECOND_END_FASTQ=" + pair2File.getAbsolutePath()
        });

        verifyFastq(pair1File, pair2File, samFile);
    }

    @Test(dataProvider =  "okFiles")
    public void testOkInterleavedFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pairFile = newTempFastqFile("pair");

        convertFile(new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "FASTQ=" + pairFile.getAbsolutePath(),
                "INTERLEAVE=true"
        });

        final Set<String> outputHeaderSet = createFastqReadHeaderSet(pairFile);
        // Create map of mate pairs from SAM records
        final Map<String,MatePair> map = createSamMatePairsMap(samFile) ;
        Assert.assertEquals(map.size() * 2, outputHeaderSet.size());

        // Ensure that each mate of each pair in SAM file is in the correct fastq pair file
        for (final Map.Entry<String,MatePair> entry : map.entrySet() ) {
            final MatePair mpair = entry.getValue();
            Assert.assertNotNull(mpair.mate1); // ensure we have two mates
            Assert.assertNotNull(mpair.mate2);
            Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
            final String readName = mpair.mate1.getReadName() ;
            Assert.assertTrue(outputHeaderSet.contains(readName+"/1")); // ensure mate is in correct file
            Assert.assertTrue(outputHeaderSet.contains(readName+"/2"));
        }
    }


    @Test (dataProvider = "badFiles", expectedExceptions= SAMFormatException.class)
    public void testBadFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pair1 = File.createTempFile("tt-pair1.", ".fastq");
        final File pair2 = File.createTempFile("tt-pair2.", ".fastq");
        pair1.deleteOnExit();
        pair2.deleteOnExit();
        convertFile(new String[]{
              "INPUT=" + samFile.getAbsolutePath(),
              "FASTQ=" + pair1.getAbsolutePath(),
              "SECOND_END_FASTQ=" + pair2.getAbsolutePath()
        });
    }

    @DataProvider(name = "okGroupedFiles")
    public Object[][] okGroupedFiles() {
        return new Object[][] {
            {"ok/grouped-last-pair-mates-flipped.sam", new String[]{"rg1","rg2"}},
        };
    }


    @DataProvider(name = "badGroupedFiles")
    public Object[][] badGroupedFiles() {
        return new Object[][] {
            {"bad/grouped-unpaired-mate.sam"}
        };
    }

    @DataProvider(name = "missingRgFiles")
    public Object[][] missingRgFiles() {
        return new Object[][] {
                {"bad/missing-rg-info.sam"}
        };
    }

    @Test(dataProvider = "okGroupedFiles")
    public void testOkGroupedFiles(final String samFilename, final String [] groupFiles) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final Map<String, Set<String>> outputSets = new HashMap<>(groupFiles.length);

        final String tmpDir = IOUtil.getDefaultTmpDir().getAbsolutePath() + "/";
        final String [] args = new String[]{
              "INPUT=" + samFile.getAbsolutePath(),
              "OUTPUT_PER_RG=true",
              "OUTPUT_DIR=" + tmpDir,
        };
        runPicardCommandLine(args);

        File f1;
        File f2;
        String fname1;
        String fname2;
        String keyName1;
        String keyName2;
        Set<String> outputHeaderSet1;
        Set<String> outputHeaderSet2;
        for(final String groupPUName : groupFiles)
        {
            keyName1 = groupPUName + "_1";
            keyName2 = groupPUName + "_2";
            fname1 = tmpDir + "/" + keyName1 + ".fastq";
            fname2 = tmpDir + "/" + keyName2 + ".fastq";
            f1 = new File(fname1);
            f2 = new File(fname2);
            f1.deleteOnExit();
            f2.deleteOnExit();
            IOUtil.assertFileIsReadable(f1);
            IOUtil.assertFileIsReadable(f2);

            // Check that paired fastq files are same size and store them for later comparison
            outputHeaderSet1 = createFastqReadHeaderSet(f1);
            outputHeaderSet2 = createFastqReadHeaderSet(f2);
            outputSets.put(keyName1 , outputHeaderSet1);
            outputSets.put(keyName2, outputHeaderSet2);
            Assert.assertEquals(outputHeaderSet1.size(), outputHeaderSet2.size());
        }

        // Create map of read groups and mate pairs from SAM records
        final Map<String, Map<String,MatePair>> map = createPUPairsMap(samFile);

        for(final Map.Entry<String, Map<String, MatePair>> groupEntry : map.entrySet()) {
            // Ensure that for each group, each mate of each pair in the SAM file is in the correct fastq pair file
            for (final Map.Entry<String,MatePair> entry : groupEntry.getValue().entrySet() ) {
                final MatePair mpair = entry.getValue();
                outputHeaderSet1 = outputSets.get(groupEntry.getKey() + "_1");
                outputHeaderSet2 = outputSets.get(groupEntry.getKey() + "_2");

                Assert.assertNotNull(mpair.mate1); // ensure we have two mates
                Assert.assertNotNull(mpair.mate2);
                Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
                final String readName = mpair.mate1.getReadName() ;
                Assert.assertTrue(outputHeaderSet1.contains(readName+"/1")); // ensure mate is in correct file
                Assert.assertTrue(outputHeaderSet2.contains(readName+"/2"));
            }
        }
    }

    @Test (dataProvider = "badGroupedFiles", expectedExceptions=SAMException.class)
    public void testBadGroupedFileOutputPerRg(final String samFilename) throws IOException {
        convertFile(new String[]{
                "INPUT=" + TEST_DATA_DIR + "/" + samFilename,
                "OUTPUT_DIR=" + IOUtil.getDefaultTmpDir().getAbsolutePath() + "/",
                "OUTPUT_PER_RG=true"
        });
    }

    @Test (dataProvider = "badGroupedFiles", expectedExceptions=SAMFormatException.class)
    public void testBadGroupedFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR, samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");

        convertFile(new String[]{
                "INPUT=" + TEST_DATA_DIR + "/" + samFilename,
                "FASTQ=" + pair1File.getAbsolutePath(),
                "SECOND_END_FASTQ=" + pair2File.getAbsolutePath()
        });
    }

    @Test (dataProvider = "missingRgFiles", expectedExceptions=PicardException.class)
    public void testMissingRgFileOutputPerRg(final String samFilename) throws IOException {
        convertFile(new String[]{
                "INPUT=" + TEST_DATA_DIR + "/" + samFilename,
                "OUTPUT_DIR=" + IOUtil.getDefaultTmpDir().getAbsolutePath() + "/",
                "OUTPUT_PER_RG=true"
        });
    }

    @Test (dataProvider = "missingRgFiles")
    public void testMissingRgFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR, samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");
        pair1File.deleteOnExit();
        pair2File.deleteOnExit();
        convertFile(new String[]{
                "INPUT=" + samFile.getAbsolutePath(),
                "FASTQ=" + pair1File.getAbsolutePath(),
                "SECOND_END_FASTQ=" + pair2File.getAbsolutePath()
        });
        verifyFastq(pair1File, pair2File, samFile);
    }

    @Test(dataProvider = "trimmedData")
    public void testTrimming(final String samFilename, final int read1Trim,
                             final int read1MaxBases, final int expectedRead1Length, final int read2Trim,
                             final int read2MaxBases, final int expectedRead2Length) throws IOException {

        final File samFile = new File(TEST_DATA_DIR, samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");
        pair1File.deleteOnExit();
        pair2File.deleteOnExit();

        convertFile(new String[]{
              "INPUT=" + samFile.getAbsolutePath(),
              "FASTQ=" + pair1File.getAbsolutePath(),
              "SECOND_END_FASTQ=" + pair2File.getAbsolutePath(),
              "READ1_TRIM=" + read1Trim,
              "READ1_MAX_BASES_TO_WRITE=" + read1MaxBases,
              "READ2_TRIM=" + read2Trim,
              "READ2_MAX_BASES_TO_WRITE=" + read2MaxBases
        });

        for (final FastqRecord first : new FastqReader(pair1File)) {
            Assert.assertEquals(first.getReadString().length(), expectedRead1Length, "Incorrect read length");
            Assert.assertEquals(first.getBaseQualityString().length(), expectedRead1Length, "Incorrect quality string length");
        }
        for (final FastqRecord second : new FastqReader(pair2File)) {
            Assert.assertEquals(second.getReadString().length(), expectedRead2Length, "Incorrect read length");
            Assert.assertEquals(second.getBaseQualityString().length(), expectedRead2Length, "Incorrect quality string length");
        }
    }

    @DataProvider(name = "trimmedData")
    public Object[][] trimmedData() {
        return new Object[][] {
            // There are 13 bases in each of these reads
            {"ok/sorted-pair.sam", 6, 7, 7, 5, 8, 8}, // exact matches for everything
            {"ok/sorted-pair.sam", 7, 7, 6, 3, 8, 8}  // fewer or more bases
        };
    }

    private Set<String> createFastqReadHeaderSet(final File file) {
        final Set<String> set = new HashSet<String>();
        final FastqReader freader = new FastqReader(file);
        while (freader.hasNext()) {
            final FastqRecord frec = freader.next();
            set.add(frec.getReadHeader());
        }
        return set ;
    }

    private Map<String,MatePair> createSamMatePairsMap(final File samFile) throws IOException {
        IOUtil.assertFileIsReadable(samFile);
        final SamReader reader = SamReaderFactory.makeDefault().open(samFile);

        final Map<String,MatePair> map = new LinkedHashMap<String,MatePair>();
        for (final SAMRecord record : reader ) {
            MatePair mpair = map.get(record.getReadName());
            if (mpair == null) {
                 mpair = new MatePair();
                 map.put(record.getReadName(), mpair);
            }
            mpair.add(record);
        }
        reader.close();
        return map;
    }


    private Map<String, Map<String, MatePair>> createPUPairsMap(final File samFile) throws IOException {
        IOUtil.assertFileIsReadable(samFile);
        final SamReader reader = SamReaderFactory.makeDefault().open(samFile);
        final Map<String, Map<String, MatePair>> map = new LinkedHashMap<String, Map<String,MatePair>>();

        Map<String,MatePair> curFileMap;
        for (final SAMRecord record : reader ) {
            final String platformUnit = record.getReadGroup().getPlatformUnit();
            curFileMap = map.get(platformUnit);
            if(curFileMap == null)
            {
                curFileMap = new LinkedHashMap<String, MatePair>();
                map.put(platformUnit, curFileMap);
            }

            MatePair mpair = curFileMap.get(record.getReadName());
            if (mpair == null) {
                 mpair = new MatePair();
                 curFileMap.put(record.getReadName(), mpair);
            }
            mpair.add(record);
        }
        reader.close();
        return map;
    }

    private void verifyFastq(final File pair1File, final File pair2File, final File samFile) throws IOException {
        // Check that paired fastq files are same size
        final Set<String> outputHeaderSet1 = createFastqReadHeaderSet(pair1File);
        final Set<String> outputHeaderSet2 = createFastqReadHeaderSet(pair2File);
        Assert.assertEquals(outputHeaderSet1.size(), outputHeaderSet2.size());

        // Create map of mate pairs from SAM records
        final Map<String,MatePair> map = createSamMatePairsMap(samFile) ;
        Assert.assertEquals(map.size(), outputHeaderSet2.size());

        // Ensure that each mate of each pair in SAM file is in the correct fastq pair file
        for (final Map.Entry<String,MatePair> entry : map.entrySet() ) {
            final MatePair mpair = entry.getValue();
            Assert.assertNotNull(mpair.mate1); // ensure we have two mates
            Assert.assertNotNull(mpair.mate2);
            Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
            final String readName = mpair.mate1.getReadName() ;
            Assert.assertTrue(outputHeaderSet1.contains(readName+"/1")); // ensure mate is in correct file
            Assert.assertTrue(outputHeaderSet2.contains(readName+"/2"));
        }
    }

    class MatePair {
        SAMRecord mate1 ;
        SAMRecord mate2 ;
        void add(final SAMRecord record) {
            if (!record.getReadPairedFlag()) throw new PicardException("Record "+record.getReadName()+" is not paired");
            if (record.getFirstOfPairFlag()) { 
                if (mate1 != null) throw new PicardException("Mate 1 already set for record: "+record.getReadName());
                mate1 = record ;
            }
            else if (record.getSecondOfPairFlag()) { 
                if (mate2 != null) throw new PicardException("Mate 2 already set for record: "+record.getReadName());
                mate2 = record ;
            }
            else throw new PicardException("Neither FirstOfPairFlag or SecondOfPairFlag is set for a paired record");
        }
    }

    private File newTempFastqFile(final String filename) throws IOException {
        if(filename == null) return null;
        final File file = File.createTempFile(filename,".fastq");
        file.deleteOnExit();
        return file; 
    }
}
