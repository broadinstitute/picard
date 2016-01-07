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
package picard.fingerprint;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;

/**
 */
public class HaplotypeMapTest {

    public static final File TEST_MAP =
            new File("testdata/picard/fingerprint/haplotypeMap.txt");

    @Test
    public void testHaplotypeMapReader() {
        HaplotypeMap map = new HaplotypeMap(TEST_MAP);
        Assert.assertEquals(map.getHaplotypes().size(), 23, "Wrong number of haplotypes returned.");
        Assert.assertEquals(map.getAllSnps().size(), 26, "Wrong number of snps returned.");

    }

    @Test
    public void testHaplotypeMapWriter() throws Exception  {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary sd = new SAMSequenceDictionary();
        sd.addSequence(new SAMSequenceRecord("chr1", 15000000));
        sd.addSequence(new SAMSequenceRecord("chr2", 15000000));
        sd.addSequence(new SAMSequenceRecord("chr3", 15000000));
        header.setSequenceDictionary(sd);


        HaplotypeMap newMap = new HaplotypeMap(header);
        HaplotypeBlock t1 = new HaplotypeBlock(0.151560926);
        t1.addSnp(new Snp("snp1", "chr1", 13969408, (byte)'T', (byte)'C',
                  0.151560926, null));
        t1.addSnp(new Snp("snp2", "chr1", 1234567, (byte)'A', (byte)'T', 1-0.151560926,
                Arrays.asList("SQNM_1CHIP_FingerprintAssays")));
        newMap.addHaplotype(t1);
        HaplotypeBlock t2 = new HaplotypeBlock(.02d);
        t2.addSnp(new Snp("snp3", "chr2", 1234567, (byte)'C', (byte)'G', .02, null));
        newMap.addHaplotype(t2);
        File temp = File.createTempFile("haplotypeMap", "txt");
        temp.deleteOnExit();
        newMap.writeToFile(temp);

        BufferedReader reader = new BufferedReader(new FileReader(temp));
        // Skip the header and sequence dictionary
        for (int i = 0; i < 5; i++) {
            reader.readLine();
        }

        String first[] = reader.readLine().split("\t");
        Assert.assertEquals(first[0], "chr1", "Wrong chromosome on first snp: " + first[0]);
        Assert.assertEquals(first[2], "snp2", "Wrong name on first snp: " + first[2]);
        Assert.assertEquals(first[6].trim(), "", "anchor snp should be null on first snp: " + first[6] );
        Assert.assertEquals(first[7], "SQNM_1CHIP_FingerprintAssays",
                "Incorrect fingerprint panel on first snp: " + first[7] );

        String second[] = reader.readLine().split("\t");
        Assert.assertEquals(second[0], "chr1", "Wrong chromosome on second snp: " + second[0]);
        Assert.assertEquals(second[2], "snp1", "Wrong name on second snp: " + second[2]);
        Assert.assertEquals(second[6], "snp2", "anchor snp is incorrect on second snp: " + second[6] );

        String third[] = reader.readLine().split("\t");
        Assert.assertEquals(third[0], "chr2", "Wrong chromosome on third snp: " + third[0]);
        Assert.assertEquals(third[2], "snp3", "Wrong name on third snp: " + third[2]);
        Assert.assertEquals(6, third.length, "Third snp should not have anchor snp or fingerprint " + Arrays.asList(third) );

    }

}