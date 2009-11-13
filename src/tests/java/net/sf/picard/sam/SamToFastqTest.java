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
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Tests for BamToFastq
 */
public class SamToFastqTest {
    private static final File TEST_DATA_DIR = new File("testdata/net/sf/picard/sam/bam2fastq/paired");

    @DataProvider(name = "okFiles")
    public Object[][] okFiles() {
        return new Object[][] {
            {"ok/sorted-pair.sam"}, // 5 sorted pairs (10 records) - mate1, mate2
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

    private void convertFile(final String samFilename, final String fastqFilename1, final String fastqFilename2) throws IOException {
        convertFile(new File(TEST_DATA_DIR, samFilename), 
            newTempFastqFile(fastqFilename1), newTempFastqFile(fastqFilename2));
    }

    private void convertFile(final File samFile, final File fastqFile1, final File fastqFile2) throws IOException {
        final SamToFastq program = new SamToFastq();
        program.INPUT = samFile ;
        program.FASTQ = fastqFile1;
        program.SECOND_END_FASTQ = fastqFile2;
        Assert.assertEquals(program.doWork(), 0);
    }

    @Test(dataProvider = "okFiles")
    public void testOkFile(final String samFilename) throws IOException {
        final File samFile = new File(TEST_DATA_DIR,samFilename);
        final File pair1File = newTempFastqFile("pair1");
        final File pair2File = newTempFastqFile("pair2");

        convertFile(samFile, pair1File, pair2File);

        // Check that paired fastq files are same size
        final Set<String> set1 = createFastqReadHeaderSet(pair1File);
        final Set<String> set2 = createFastqReadHeaderSet(pair2File);
        Assert.assertEquals(set1.size(), set2.size());

        // Create map of mate pairs from SAM records
        final Map<String,MatePair> map = createSamMatePairsMap(samFile) ;
        Assert.assertEquals(map.size(), set2.size());

        // Ensure that each mate of each pair in SAM file is in the correct fastq pair file
        for (final Map.Entry<String,MatePair> entry : map.entrySet() ) {
            final MatePair mpair = entry.getValue();
            Assert.assertNotNull(mpair.mate1); // ensure we have two mates
            Assert.assertNotNull(mpair.mate2);
            Assert.assertEquals(mpair.mate1.getReadName(),mpair.mate2.getReadName());
            final String readName = mpair.mate1.getReadName() ;
            Assert.assertTrue(set1.contains(readName+"/1")); // ensure mate is in correct file
            Assert.assertTrue(set2.contains(readName+"/2"));
        }
    }

    @Test (dataProvider = "badFiles", expectedExceptions= PicardException.class)
    public void testBadFile(final String samFilename) throws IOException {
        convertFile(samFilename, "tt-pair1.fastq", "tt-pair2.fastq");
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
        IoUtil.assertFileIsReadable(samFile);
        final SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(samFile));

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
        final File file = File.createTempFile(filename,".fastq");
        file.deleteOnExit();
        return file; 
    }
}
