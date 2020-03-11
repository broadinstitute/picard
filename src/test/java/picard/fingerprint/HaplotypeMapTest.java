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
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;

/**
 */
public class HaplotypeMapTest {

    public static final File TEST_MAP =
            new File("testdata/picard/fingerprint/haplotypeMap.txt");

    public static final File TEST_VCF_MAP =
            new File("testdata/picard/fingerprint/haplotypeMap.vcf");

    public static final File TEST_FASTA =
            new File("testdata/picard/fingerprint/reference.fasta");


    @DataProvider(name="haplotypeMapReaderData")
    public Object[][] haplotypeMapReaderData(){
        return new Object[][]{ {TEST_MAP}, {TEST_VCF_MAP}};
    }
    @Test(dataProvider = "haplotypeMapReaderData")
    public void testHaplotypeMapReader(final File file) {
        final HaplotypeMap map = new HaplotypeMap(file);
        Assert.assertEquals(map.getHaplotypes().size(), 23, "Wrong number of haplotypes returned.");
        Assert.assertEquals(map.getAllSnps().size(), 26, "Wrong number of snps returned.");
    }

    @Test
    public void testEquivalenceHapDBandVCF() {
        final HaplotypeMap mapHapDB = new HaplotypeMap(TEST_MAP);
        final HaplotypeMap mapVCF = new HaplotypeMap(TEST_VCF_MAP);

        for (final HaplotypeBlock haplotypeBlock: mapHapDB.getHaplotypes()){
            final Snp firstSnp = haplotypeBlock.getFirstSnp();

            //names in the two schemes are different, so need to use position
            HaplotypeBlock haplotypeBlockVcf = mapVCF.getHaplotype(firstSnp.getChrom(), firstSnp.getPos());
            Assert.assertNotNull(haplotypeBlockVcf);

            Assert.assertEquals(haplotypeBlock, haplotypeBlockVcf); // only checks that domain is right.

            for(final Snp snp: haplotypeBlock.getSnps()){

                final Snp otherSnp = haplotypeBlockVcf.getSnps().stream()
                        .filter(s->s.getChrom().equals(snp.getChrom()) && s.getPos()==snp.getPos())
                        .findFirst()
                        .orElseThrow(()->new AssertionError("missing snp:"+ snp.toString()));

                Assert.assertEquals(otherSnp.getAllele1(), snp.getAllele1());
                Assert.assertEquals(otherSnp.getAllele2(), snp.getAllele2());
            }
        }
    }


    @DataProvider(name="haplotypeMapForWriting")
    public Object[][] haplotypeMapForWriting() {

        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        SAMSequenceDictionary sd = new SAMSequenceDictionary();
        sd.addSequence(new SAMSequenceRecord("chr1", 101));
        sd.addSequence(new SAMSequenceRecord("chr2", 101));
        sd.addSequence(new SAMSequenceRecord("chr3", 101));
        header.setSequenceDictionary(sd);


        HaplotypeMap newMap = new HaplotypeMap(header);
        HaplotypeBlock t1 = new HaplotypeBlock(0.151560926);
        t1.addSnp(new Snp("snp2", "chr1", 83, (byte)'A', (byte)'G', .16, null));
        t1.addSnp(new Snp("snp1", "chr1", 51, (byte)'T', (byte)'C', 0.15, Collections.singletonList("SQNM_1CHIP_FingerprintAssays")));
        newMap.addHaplotype(t1);
        HaplotypeBlock t2 = new HaplotypeBlock(.02d);
        t2.addSnp(new Snp("snp3", "chr2", 24, (byte)'C', (byte)'A', .20, null));
        newMap.addHaplotype(t2);

        return new Object[][]{{newMap}};
    }

    @Test(dataProvider = "haplotypeMapForWriting")
    public void testHaplotypeMapWriteToVcf(final HaplotypeMap haplotypeMap) throws Exception {
        final File temp = File.createTempFile("haplotypeMap", ".vcf");
        temp.deleteOnExit();
        haplotypeMap.writeAsVcf(temp, TEST_FASTA);

        final VCFFileReader reader = new VCFFileReader(temp);

        Assert.assertEquals(reader.getFileHeader().getNGenotypeSamples(), 1, "VCF should have exactly one sample");
        Assert.assertEquals(reader.getFileHeader().getSampleNamesInOrder().get(0), HaplotypeMap.HET_GENOTYPE_FOR_PHASING, "VCF sole sample should be " + HaplotypeMap.HET_GENOTYPE_FOR_PHASING);

        final Iterator<VariantContext> iter = reader.iterator();
        final VariantContext first = iter.next();
        Assert.assertEquals(first.getContig(), "chr1", "Wrong chromosome on first snp: " + first);
        Assert.assertEquals(first.getID(), "snp1", "Wrong name on first snp: " + first);
        Assert.assertEquals(first.getGenotype(0).getExtendedAttribute(VCFConstants.PHASE_SET_KEY), Integer.toString(first.getStart()), "anchor snp should have PS equal to its position " + first);
        Assert.assertEquals(first.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0D), 1 - 0.15); // because it's swapped w.r.t the reference

        final VariantContext second = iter.next();
        Assert.assertEquals(second.getContig(), "chr1", "Wrong chromosome on second snp: " + second);
        Assert.assertEquals(second.getID(), "snp2", "Wrong name on second snp: " + second);
        Assert.assertEquals(second.getGenotype(0).getExtendedAttribute(VCFConstants.PHASE_SET_KEY), Integer.toString(first.getStart()), "Phase set is incorrect on second snp: " + second);
        Assert.assertEquals(second.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0D), 0.16);

        final VariantContext third = iter.next();
        Assert.assertEquals(third.getContig(), "chr2", "Wrong chromosome on third snp: " + third);
        Assert.assertEquals(third.getID(), "snp3", "Wrong name on third snp: " + third);
        Assert.assertFalse (third.getGenotype(0).hasAnyAttribute(VCFConstants.PHASE_SET_KEY), "Third snp should not have a phaseset" + third);
        Assert.assertEquals(third.getAttributeAsDouble(VCFConstants.ALLELE_FREQUENCY_KEY, 0D), 0.2);
    }


    @Test(dataProvider = "haplotypeMapForWriting")
    public void testHaplotypeMapWriter(final HaplotypeMap haplotypeMap) throws Exception  {

        File temp = File.createTempFile("haplotypeMap", "txt");
        temp.deleteOnExit();
        haplotypeMap.writeToFile(temp);

        BufferedReader reader = new BufferedReader(new FileReader(temp));
        // Skip the header and sequence dictionary
        for (int i = 0; i < 5; i++) {
            reader.readLine();
        }

        String first[] = reader.readLine().split("\t");
        Assert.assertEquals(first[0], "chr1", "Wrong chromosome on first snp: " + first[0]);
        Assert.assertEquals(first[2], "snp1", "Wrong name on first snp: " + first[2]);
        Assert.assertEquals(Double.parseDouble(first[5]), 0.15, "Wrong AlleleFrequency on first snp: " + first[5]);
        Assert.assertEquals(first[6].trim(), "", "anchor snp should be null on first snp: " + first[6] );
        Assert.assertEquals(first[7], "SQNM_1CHIP_FingerprintAssays",
                "Incorrect fingerprint panel on first snp: " + first[7] );

        String second[] = reader.readLine().split("\t");
        Assert.assertEquals(second[0], "chr1", "Wrong chromosome on second snp: " + second[0]);
        Assert.assertEquals(second[2], "snp2", "Wrong name on second snp: " + second[2]);
        Assert.assertEquals(Double.parseDouble(second[5]), 0.16, "Wrong AlleleFrequency on second snp: " + second[5]);
        Assert.assertEquals(second[6], "snp1", "anchor snp is incorrect on second snp: " + second[6] );

        String third[] = reader.readLine().split("\t");
        Assert.assertEquals(third[0], "chr2", "Wrong chromosome on third snp: " + third[0]);
        Assert.assertEquals(third[2], "snp3", "Wrong name on third snp: " + third[2]);
        Assert.assertEquals(Double.parseDouble(third[5]), 0.20, "Wrong AlleleFrequency on second snp: " + second[5]);
        Assert.assertEquals(6, third.length, "Third snp should not have anchor snp or fingerprint " + Arrays.asList(third) );
    }
}
