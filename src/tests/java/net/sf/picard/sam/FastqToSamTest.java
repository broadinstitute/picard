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
package net.sf.picard.sam;

import net.sf.picard.PicardException;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.util.FastqQualityFormat;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;

/**
 * Tests for FastqToBam
 */
public class FastqToSamTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam/fastq2bam");

    // fastq files with legal values for each fastq version
    @DataProvider(name = "okVersionFiles")
    public Object[][] okVersionFiles() {
        return new Object[][] {
            {"fastq-sanger/5k-v1-Rhodobacter_LW1.sam.fastq",      FastqQualityFormat.Standard },
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


    @Test(dataProvider = "okVersionFiles")
    public void testFastqVersionOk(final String fastqVersionFilename, final FastqQualityFormat version) throws IOException {
        final File fastqVersionSamFile = convertFile(fastqVersionFilename, version);
    }

    @Test(dataProvider = "badVersionFiles", expectedExceptions= IllegalArgumentException.class)
    public void testFastqVersionBad(final String fastqVersionFilename, final FastqQualityFormat version) throws IOException {
        final File fastqVersionSamFile = convertFile(fastqVersionFilename, version);
    }

    @Test(dataProvider = "badFormatFiles", expectedExceptions= PicardException.class) 
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

    private File convertFile(final String fastqFilename1, final String fastqFilename2, final FastqQualityFormat version) throws IOException {
        final File fastqFile1 = new File(TEST_DATA_DIR, fastqFilename1);
        final File samFile = newTempSamFile(fastqFile1.getName());

        final FastqToSam program = new FastqToSam();
        program.FASTQ = fastqFile1;
        if (fastqFilename2 != null) program.FASTQ2 = new File(TEST_DATA_DIR, fastqFilename2);
        program.OUTPUT = samFile;
        program.QUALITY_FORMAT = version;
        program.READ_GROUP_NAME = "rg" ;
        program.SAMPLE_NAME = "s1" ;
        Assert.assertEquals(program.doWork(), 0);
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

    private static final FastqToSam fastqToSam = new FastqToSam();
    private static File dummyFile ;
    private static FastqReader freader1 ;
    private static FastqReader freader2 ;

    @BeforeClass
    public static void beforeClass() throws IOException {
        dummyFile = newTempFile("dummy");
        freader1 = new FastqReader(dummyFile);
        freader2 = new FastqReader(dummyFile);
    }

    @DataProvider(name = "okPairNames")
    public Object[][] okPairNames() {
        return new Object[][] {
            {"aa/1", "aa/2" },
            {"aa/bb/1", "aa/bb/2" },
            {"aa/bb/cc/dd/ee/ff/1", "aa/bb/cc/dd/ee/ff/2" },
            {"////1", "////2" },
        };
    }

    @DataProvider(name = "badPairNames")
    public Object[][] badPairNames() {
        return new Object[][] {
            {"", "" },
            {"/", "/" },
            {"////", "////" },
            {"aa", "aa" },
            {"/aa", "/aa" },
            {"aa/", "aa/" },
            {"aa/1", "aa/1" },
            {"aa/2", "aa/2" },
            {"aa/1", "bb/2" },
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
}

