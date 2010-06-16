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

package net.sf.picard.reference;

import net.sf.picard.PicardException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import org.testng.annotations.Test;
import org.testng.annotations.DataProvider;
import org.testng.Assert;

/**
 * Test the fasta sequence index reader.
 */
public class FastaSequenceIndexTest {
    private static File TEST_DATA_DIR = new File("testdata/net/sf/picard/reference");

    @DataProvider(name="homosapiens")
    public Object[][] provideHomoSapiens() throws FileNotFoundException {
        final File sequenceIndexFile = new File(TEST_DATA_DIR,"Homo_sapiens_assembly18.fasta.fai");
        return new Object[][] { new Object[] { new FastaSequenceIndex(sequenceIndexFile) } };    
    }

    @DataProvider(name="specialcharacters")
    public Object[][] provideSpecialCharacters() throws FileNotFoundException {
        final File sequenceIndexFile = new File(TEST_DATA_DIR,"testing.fai");
        return new Object[][] { new Object[] { new FastaSequenceIndex(sequenceIndexFile) } };
    }

    @Test(dataProvider="homosapiens")
    public void testInitialContig(FastaSequenceIndex sequenceIndex) {
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chrM"),"Contig chrM is not present");
        FastaSequenceIndexEntry entry = sequenceIndex.getIndexEntry("chrM");
        Assert.assertEquals(entry.getContig(),"chrM","Contig chrM name is incorrect");
        Assert.assertEquals(entry.getLocation(),6L,"Contig chrM location is incorrect");
        Assert.assertEquals(entry.getSize(),16571L,"Contig chrM size is incorrect");
        Assert.assertEquals(entry.getBasesPerLine(),50,"Contig chrM bases per line is incorrect");
        Assert.assertEquals(entry.getBytesPerLine(),51,"Contig chrM bytes per line is incorrect");
    }

    @Test(dataProvider="homosapiens")
    public void testMiddleContig(FastaSequenceIndex sequenceIndex) {
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr8"),"Contig chr8 is not present");
        FastaSequenceIndexEntry entry = sequenceIndex.getIndexEntry("chr8");
        Assert.assertEquals(entry.getContig(),"chr8","Contig chr8 name is incorrect");
        Assert.assertEquals(entry.getLocation(),1419403101L,"Contig chr8 location is incorrect");
        Assert.assertEquals(entry.getSize(),146274826L,"Contig chr8 size is incorrect");
        Assert.assertEquals(entry.getBasesPerLine(),50,"Contig chr8 bases per line is incorrect");
        Assert.assertEquals(entry.getBytesPerLine(),51,"Contig chr8 bytes per line is incorrect");
    }

    @Test(dataProvider="homosapiens")
    public void testLastContig(FastaSequenceIndex sequenceIndex) {
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chrX_random"),"Contig chrX_random is not present");
        FastaSequenceIndexEntry entry = sequenceIndex.getIndexEntry("chrX_random");
        Assert.assertEquals(entry.getContig(),"chrX_random","Contig chrX_random name is incorrect");
        Assert.assertEquals(entry.getLocation(),3156698441L,"Contig chrX_random location is incorrect");
        Assert.assertEquals(entry.getSize(),1719168L,"Contig chrX_random size is incorrect");
        Assert.assertEquals(entry.getBasesPerLine(),50,"Contig chrX_random bases per line is incorrect");
        Assert.assertEquals(entry.getBytesPerLine(),51,"Contig chrX_random bytes per line is incorrect");
    }

    @Test(dataProvider="homosapiens")
    public void testAllContigsPresent(FastaSequenceIndex sequenceIndex) {
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chrM"),"Contig chrM is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr1"),"Contig chr1 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr2"),"Contig chr2 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr3"),"Contig chr3 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr4"),"Contig chr4 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr5"),"Contig chr5 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr6"),"Contig chr6 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr7"),"Contig chr7 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr8"),"Contig chr8 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr9"),"Contig chr9 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr10"),"Contig chr10 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr11"),"Contig chr11 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr12"),"Contig chr12 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr13"),"Contig chr13 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr14"),"Contig chr14 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr15"),"Contig chr15 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr16"),"Contig chr16 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr17"),"Contig chr17 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr18"),"Contig chr18 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr19"),"Contig chr19 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr20"),"Contig chr20 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr21"),"Contig chr21 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr22"),"Contig chr22 is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chrX"),"Contig chrX is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chrY"),"Contig chrY is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr1_random"),"Contig chr1_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr2_random"),"Contig chr2_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr3_random"),"Contig chr3_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr4_random"),"Contig chr4_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr5_random"),"Contig chr5_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr6_random"),"Contig chr6_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr7_random"),"Contig chr7_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr8_random"),"Contig chr8_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr9_random"),"Contig chr9_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr10_random"),"Contig chr10_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr11_random"),"Contig chr11_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr13_random"),"Contig chr13_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr15_random"),"Contig chr15_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr16_random"),"Contig chr16_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr17_random"),"Contig chr17_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr18_random"),"Contig chr18_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr19_random"),"Contig chr19_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr21_random"),"Contig chr21_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chr22_random"),"Contig chr22_random is not present");
        Assert.assertTrue(sequenceIndex.hasIndexEntry("chrX_random"),"Contig chrX_random is not present");
    }

    @Test(dataProvider="homosapiens")
    public void testHasInvalidEntry(FastaSequenceIndex sequenceIndex) {
        Assert.assertFalse(sequenceIndex.hasIndexEntry("invalid"),"Found an invalid entry");
    }

    @Test(dataProvider="homosapiens",expectedExceptions=PicardException.class)
    public void testGetInvalidEntry(FastaSequenceIndex sequenceIndex) {
        sequenceIndex.getIndexEntry("invalid");
    }

    @Test(dataProvider="homosapiens")
    public void testIteration(FastaSequenceIndex sequenceIndex) {
        Iterator<FastaSequenceIndexEntry> sequenceIndexEntries = sequenceIndex.iterator();

        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chrM","Contig chrM is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr1","Contig chr1 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr2","Contig chr2 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr3","Contig chr3 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr4","Contig chr4 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr5","Contig chr5 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr6","Contig chr6 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr7","Contig chr7 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr8","Contig chr8 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr9","Contig chr9 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr10","Contig chr10 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr11","Contig chr11 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr12","Contig chr12 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr13","Contig chr13 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr14","Contig chr14 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr15","Contig chr15 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr16","Contig chr16 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr17","Contig chr17 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr18","Contig chr18 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr19","Contig chr19 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr20","Contig chr20 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr21","Contig chr21 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr22","Contig chr22 is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chrX","Contig chrX is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chrY","Contig chrY is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr1_random","Contig chr1_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr2_random","Contig chr2_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr3_random","Contig chr3_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr4_random","Contig chr4_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr5_random","Contig chr5_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr6_random","Contig chr6_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr7_random","Contig chr7_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr8_random","Contig chr8_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr9_random","Contig chr9_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr10_random","Contig chr10_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr11_random","Contig chr11_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr13_random","Contig chr13_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr15_random","Contig chr15_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr16_random","Contig chr16_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr17_random","Contig chr17_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr18_random","Contig chr18_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr19_random","Contig chr19_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr21_random","Contig chr21_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chr22_random","Contig chr22_random is not present");
        Assert.assertEquals(sequenceIndexEntries.next().getContig(),"chrX_random","Contig chrX_random is not present");
        Assert.assertFalse(sequenceIndexEntries.hasNext(),"Iterator still has more entries");
    }

    @Test(dataProvider="specialcharacters")
    public void testSpecialCharacters(FastaSequenceIndex specialCharactersIndex) {
        /* file contents:
        chrM	16571	6	50	51
        chr1;boat	247249719	16915	50	51
        chr2:money	242951149	252211635	50	51
        chr3::;	199501827	500021813	50	51
        ;;;;;;  1234            1234            1234    1234
        file:gi|17981852|ref|NC_001807.4|    16571   2911876801      70      71
        */
        Iterator<FastaSequenceIndexEntry> sequenceIndexEntries = specialCharactersIndex.iterator();
        FastaSequenceIndexEntry ent = sequenceIndexEntries.next();
        Assert.assertEquals(ent.getContig(),"chrM","Contig chrM is not present");
        Assert.assertEquals(ent.getSize(),16571,"Contig chrM size is not correct");
        Assert.assertEquals(ent.getLocation(),6,"Contig chrM location is not correct");
        Assert.assertEquals(ent.getBasesPerLine(),50,"Contig chrM bases per line is not correct");
        Assert.assertEquals(ent.getBytesPerLine(),51,"Contig chrM bytes per line is not correct");

        ent = sequenceIndexEntries.next();
        Assert.assertEquals(ent.getContig(),"chr1;boat","Contig chr1;boat is not present");
        Assert.assertEquals(ent.getSize(),247249719,"Contig chr1;boat size is not correct");
        Assert.assertEquals(ent.getLocation(),16915,"Contig chr1;boat location is not correct");
        Assert.assertEquals(ent.getBasesPerLine(),50,"Contig chr1;boat bases per line is not correct");
        Assert.assertEquals(ent.getBytesPerLine(),51,"Contig chr1;boat bytes per line is not correct");

        ent = sequenceIndexEntries.next();
        Assert.assertEquals(ent.getContig(),"chr2:money","Contig chr2:money is not present");
        Assert.assertEquals(ent.getSize(),242951149,"Contig chr2:money size is not correct");
        Assert.assertEquals(ent.getLocation(),252211635,"Contig chr2:money location is not correct");
        Assert.assertEquals(ent.getBasesPerLine(),50,"Contig chr2:money bases per line is not correct");
        Assert.assertEquals(ent.getBytesPerLine(),51,"Contig chr2:money bytes per line is not correct");

        ent = sequenceIndexEntries.next();
        Assert.assertEquals(ent.getContig(),"chr3::;","Contig chr3::; is not present");
        Assert.assertEquals(ent.getSize(),199501827,"Contig chr3::; size is not correct");
        Assert.assertEquals(ent.getLocation(),500021813,"Contig chrM location is not correct");
        Assert.assertEquals(ent.getBasesPerLine(),50,"Contig chr3::; bases per line is not correct");
        Assert.assertEquals(ent.getBytesPerLine(),51,"Contig chr3::; bytes per line is not correct");

        ent = sequenceIndexEntries.next();
        Assert.assertEquals(ent.getContig(),";;;;;;;;","Contig ;;;;;;;; is not present");
        Assert.assertEquals(ent.getSize(),123,"Contig ;;;;;;;; size is not correct");
        Assert.assertEquals(ent.getLocation(),234,"Contig ;;;;;;;; location is not correct");
        Assert.assertEquals(ent.getBasesPerLine(),456,"Contig ;;;;;;;; bases per line is not correct");
        Assert.assertEquals(ent.getBytesPerLine(),789,"Contig ;;;;;;;; bytes per line is not correct");

        ent = sequenceIndexEntries.next();
        Assert.assertEquals(ent.getContig(),"file:gi|17981852|ref|NC_001807.4|","Contig file:gi|17981852|ref|NC_001807.4| is not present");
        Assert.assertEquals(ent.getSize(),16571,"Contig file:gi|17981852|ref|NC_001807.4| size is not correct");
        Assert.assertEquals(ent.getLocation(),2911876801L,"Contig file:gi|17981852|ref|NC_001807.4| location is not correct");
        Assert.assertEquals(ent.getBasesPerLine(),70,"Contig file:gi|17981852|ref|NC_001807.4| bases per line is not correct");
        Assert.assertEquals(ent.getBytesPerLine(),71,"Contig file:gi|17981852|ref|NC_001807.4| bytes per line is not correct");
    }
}
