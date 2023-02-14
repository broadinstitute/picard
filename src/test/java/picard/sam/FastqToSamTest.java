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
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.util.FastqQualityFormat;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.cmdline.CommandLineProgramTest;
import picard.PicardException;

import java.io.File;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 * Tests for FastqToBam
 */
public class FastqToSamTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File("testdata/picard/sam/fastq2bam");

    public String getCommandLineProgramName() {
        return FastqToSam.class.getSimpleName();
    }

    // fastq files with legal values for each fastq version
    @DataProvider(name = "okVersionFiles")
    public Object[][] okVersionFiles() {
        return new Object[][] {
            {"fastq-sanger/5k-v1-Rhodobacter_LW1.sam.fastq.gz",      FastqQualityFormat.Standard },
            {"fastq-sanger/5k-30BB2AAXX.3.aligned.sam.fastq",     FastqQualityFormat.Standard },
            {"fastq-sanger/sanger_full_range_as_sanger-63.fastq", FastqQualityFormat.Standard }, // all sanger chars

            {"fastq-solexa/s_1_sequence.txt", FastqQualityFormat.Solexa},
            {"fastq-solexa/solexa_full_range_as_solexa.fastq", FastqQualityFormat.Solexa}, // all solexa chars

            {"fastq-illumina/s_1_2_sequence.txt", FastqQualityFormat.Illumina},
            {"fastq-illumina/s_1_2_sequence.txt", FastqQualityFormat.Illumina},
            {"fastq-illumina/s_1_sequence.txt", FastqQualityFormat.Illumina},
        };
    }

    // Illegal values fastq files for each fastq version
    @DataProvider(name = "badVersionFiles")
    public Object[][] badVersionFiles() {
        return new Object[][] {
            {"fastq-sanger/sanger_full_range_as_sanger-63.fastq", FastqQualityFormat.Illumina},
            {"fastq-solexa/s_1_sequence.txt", FastqQualityFormat.Illumina},
        };
    }

    // Illegal fastq format, i.e. doesn't contain correct four lines per record
    @DataProvider(name = "badFormatFiles")
    public Object[][] badFormatFiles() {
        return new Object[][] {
            {"bad-format/bad-qual-header.txt"},
            {"bad-format/bad-seq-header.txt"},
            {"bad-format/extra-line.txt"},
            {"bad-format/too-many-quals.txt"},
            {"bad-format/1lines.txt"},
            {"bad-format/2lines.txt"},
            {"bad-format/3lines.txt"},
        };
    }

    // permissive fastq format, i.e. contain blank lines at various places
    @DataProvider(name = "permissiveFormatFiles")
    public Object[][] permissiveFormatFiles() {
        return new Object[][] {
                {"permissive-format/pair1.txt",          "permissive-format/pair2.txt", FastqQualityFormat.Standard },
                {"permissive-format/s_1_1_sequence.txt",    "permissive-format/s_1_2_sequence.txt", FastqQualityFormat.Illumina},
                {"permissive-format/pair1.txt", null, FastqQualityFormat.Standard},
                {"permissive-format/pair2.txt", null, FastqQualityFormat.Standard},
                {"permissive-format/s_1_1_sequence.txt", null, FastqQualityFormat.Illumina},
                {"permissive-format/s_1_2_sequence.txt", null, FastqQualityFormat.Illumina},
                {"permissive-format/s_1_sequence.txt", null, FastqQualityFormat.Illumina},

        };
    }


    // OK paired fastq files
    @DataProvider(name = "okPairedFiles")
    public Object[][] okPairedFiles() {
        return new Object[][] {
            {"ok-paired/pair1.txt",          "ok-paired/pair2.txt", FastqQualityFormat.Standard },
            {"fastq-illumina/s_1_1_sequence.txt", "fastq-illumina/s_1_2_sequence.txt", FastqQualityFormat.Illumina}
        };
    }

    // Inconsistent paired fastq files
    @DataProvider(name = "badPairedFiles")
    public Object[][] badPairedFiles() {
        return new Object[][] {
            {"ok-paired/pair1.txt",          "bad-paired/pair2-one-more-record.txt" },
            {"bad-paired/pair1-one-more-record.txt", "ok-paired/pair2.txt" },
            {"ok-paired/pair1.txt",          "bad-paired/pair2-badnum.txt" },
            {"bad-paired/pair1-badnum.txt",   "ok-paired/pair2.txt" },
            {"bad-paired/pair1-nonum.txt",    "ok-paired/pair2.txt" },
            {"bad-paired/pair1-onetoken.txt", "ok-paired/pair2.txt" },
        };
    }

    @Test(dataProvider = "permissiveFormatFiles")
    public void testPermissiveOk(final String filename1, final String filename2, final FastqQualityFormat version) throws IOException {
        convertFile(filename1, filename2, version, true);
    }

    @Test(dataProvider = "permissiveFormatFiles",expectedExceptions = SAMException.class)
    public void testPermissiveFail(final String filename1, final String filename2, final FastqQualityFormat version) throws IOException {
        convertFile(filename1,filename2,version,false);
    }

    @Test(dataProvider = "okVersionFiles")
    public void testFastqVersionOk(final String fastqVersionFilename, final FastqQualityFormat version) throws IOException {
        convertFile(fastqVersionFilename, version);
    }

    @Test(dataProvider = "badVersionFiles", expectedExceptions= SAMException.class)
    public void testFastqVersionBad(final String fastqVersionFilename, final FastqQualityFormat version) throws IOException {
        convertFile(fastqVersionFilename, version);
    }

    @Test(dataProvider = "badFormatFiles", expectedExceptions= SAMException.class)
    public void testBadFile(final String filename) throws IOException {
        convertFile(filename, null, FastqQualityFormat.Standard);
    }

    @Test(dataProvider = "badPairedFiles", expectedExceptions= PicardException.class)
    public void testPairedBad(final String filename1, final String filename2) throws IOException {
        convertFile(filename1, filename2, FastqQualityFormat.Standard);
    }

    @Test(dataProvider = "okPairedFiles")
    public void testPairedOk(final String filename1, final String filename2, final FastqQualityFormat version) throws IOException {
        convertFile(filename1, filename2, version);
    }

    private File convertFile(final String filename, final FastqQualityFormat version) throws IOException {
        return convertFile(filename, null, version);
    }

    private File convertFile(final String fastqFilename1, final String fastqFilename2, final FastqQualityFormat version) throws IOException{
        return convertFile(fastqFilename1, fastqFilename2, version,false);
    }

    private File convertFile(final String fastqFilename1, final String fastqFilename2, final FastqQualityFormat version,final boolean permissiveFormat) throws IOException {
        return convertFile(fastqFilename1, fastqFilename2, version, permissiveFormat, false);
    }

    private File convertFile(final String fastqFilename1,
                             final String fastqFilename2, 
                             final FastqQualityFormat version,
                             final boolean permissiveFormat,
                             final boolean useSequentialFastqs) throws IOException {
        final File fastq1 = new File(TEST_DATA_DIR, fastqFilename1);
        final File fastq2 = (fastqFilename2 != null) ? new File(TEST_DATA_DIR, fastqFilename2) : null;
        final File samFile = newTempSamFile(fastq1.getName());

        final List<String> args = new ArrayList<>();

        args.add("FASTQ=" + fastq1.getAbsolutePath());
        args.add("OUTPUT=" + samFile.getAbsolutePath());
        args.add("QUALITY_FORMAT=" + version);
        args.add("READ_GROUP_NAME=rg");
        args.add("SAMPLE_NAME=s1");

        if (fastqFilename2 != null) args.add("FASTQ2=" + fastq2.getAbsolutePath());
        if (permissiveFormat) args.add("ALLOW_AND_IGNORE_EMPTY_LINES=true");
        if (useSequentialFastqs) args.add("USE_SEQUENTIAL_FASTQS=true");

        Assert.assertEquals(runPicardCommandLine(args), 0);
        return samFile ;
    
    }

    private static File newTempSamFile(final String filename) throws IOException {
        final File file = File.createTempFile(filename,".sam");
        file.deleteOnExit();
        return file; 
    }

    private static File newTempFile(final String filename) throws IOException {
        final File file = File.createTempFile(filename,".tmp");
        file.deleteOnExit();
        return file; 
    }

//  Test for legal syntax for pair read names for FastqToSam.getBaseName()
//  We create a dummy file to test the getBaseName() method since it expects 
//  an existing file.

    // TODO - Should switch over to using invocation via new PicardCommandLine() - BUT the tests using this are accessing class members directly.
    private static final FastqToSam fastqToSam = new FastqToSam();
    private static FastqReader freader1 ;
    private static FastqReader freader2 ;

    @BeforeClass
    public static void beforeClass() throws IOException {
        final File dummyFile = newTempFile("dummy");
        freader1 = new FastqReader(dummyFile);
        freader2 = new FastqReader(dummyFile);
    }

    @DataProvider(name = "okPairNames")
    public Object[][] okPairNames() {
        return new Object[][] {
            {"aa/1", "aa/2" },
            {"aa", "aa" },
            {"aa/bb", "aa/bb" },
            {"aa/bb/", "aa/bb/" },
            {"aa/bb/1", "aa/bb/2" },
            {"aa/bb/cc/dd/ee/ff/1", "aa/bb/cc/dd/ee/ff/2" },
            {"////1", "////2" },
            {"/", "/" },
            {"////", "////" },
            {"/aa", "/aa" },
            {"aa/", "aa/" },
            {"ab/c", "ab/c"}
        };
    }

    @DataProvider(name = "badPairNames")
    public Object[][] badPairNames() {
        return new Object[][] {
            {"", "" },
            {"aa/1", "bb/2" },
            {"aa"  , "bb" },
            {"aa/1", "aa" },
            {"aa",   "aa/2" },
            {"aa/1", "aa/1" },
            {"aa/2", "aa/2" },
        };
    }

    @Test(dataProvider = "okPairNames")
    public void readPairNameOk(final String name1, final String name2) throws IOException {
        fastqToSam.getBaseName(name1, name2, freader1, freader2);
    }

    @Test(dataProvider = "badPairNames", expectedExceptions= PicardException.class) 
    public void readPairNameBad(final String name1, final String name2) throws IOException {
        fastqToSam.getBaseName(name1, name2, freader1, freader2);
    }

    // TODO: convert other tests to use this too
    private void convertFileAndVerifyRecordCount(final int expectedCount,
                       final String fastqFilename1,
                       final String fastqFilename2,
                       final FastqQualityFormat version,
                       final boolean permissiveFormat,
                       final boolean useSequentialFastqs) throws IOException {
        final File samFile = convertFile(fastqFilename1, fastqFilename2, version, permissiveFormat, useSequentialFastqs);
        final SamReader samReader = SamReaderFactory.makeDefault().open(samFile);
        final SAMRecordIterator iterator = samReader.iterator();
        int actualCount = 0;
        while (iterator.hasNext()) {
            iterator.next();
            actualCount++;
        }
        samReader.close();
        Assert.assertEquals(actualCount, expectedCount);
    }
    
    @Test
    public void testSequentialFiles() throws IOException {
        final String singleEnd = "sequential-files/single_end_R1_001.fastq";
        final String pairedEnd1 = "sequential-files/paired_end_R1_001.fastq";
        final String pairedEnd2 = "sequential-files/paired_end_R2_001.fastq";
        
        Assert.assertEquals(FastqToSam.getSequentialFileList(Paths.get(TEST_DATA_DIR.getPath(), singleEnd)).size(), 2);
        Assert.assertEquals(FastqToSam.getSequentialFileList(Paths.get(TEST_DATA_DIR.getPath(),  pairedEnd1)).size(), 2);
        Assert.assertEquals(FastqToSam.getSequentialFileList(Paths.get(TEST_DATA_DIR.getPath(), pairedEnd2)).size(), 2);

        convertFileAndVerifyRecordCount(1, singleEnd, null, FastqQualityFormat.Illumina, true, false);
        convertFileAndVerifyRecordCount(2, singleEnd, null, FastqQualityFormat.Illumina, true, true);
        convertFileAndVerifyRecordCount(2, pairedEnd1, pairedEnd2, FastqQualityFormat.Illumina, true, false);
        convertFileAndVerifyRecordCount(4, pairedEnd1, pairedEnd2, FastqQualityFormat.Illumina, true, true);
    }
}

