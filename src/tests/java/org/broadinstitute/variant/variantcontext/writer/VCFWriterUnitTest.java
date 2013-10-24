/*
* Copyright (c) 2012 The Broad Institute
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.variant.variantcontext.writer;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.util.TestUtil;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.Tribble;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderVersion;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * @author aaron
 *         <p/>
 *         Class VCFWriterUnitTest
 *         <p/>
 *         This class tests out the ability of the VCF writer to correctly write VCF files
 */
public class VCFWriterUnitTest extends VariantBaseTest {
    private Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
    private Set<String> additionalColumns = new HashSet<String>();
    private File fakeVCFFile = new File("FAKEVCFFILEFORTESTING.vcf");

    /** test, using the writer and reader, that we can output and input a VCF file without problems */
    @Test
    public void testBasicWriteAndRead() {
        VCFHeader header = createFakeHeader(metaData,additionalColumns);
        final EnumSet<Options> options = EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        VariantContextWriter writer = VariantContextWriterFactory.create(fakeVCFFile, createArtificialSequenceDictionary(), options);
        writer.writeHeader(header);
        writer.add(createVC(header));
        writer.add(createVC(header));
        writer.close();
        VCFCodec codec = new VCFCodec();
        VCFHeader headerFromFile = null;
        FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(fakeVCFFile.getAbsolutePath(), codec, false);
        headerFromFile = (VCFHeader)reader.getHeader();

        int counter = 0;

        // validate what we're reading in
        validateHeader(headerFromFile);
        
        try {
            Iterator<VariantContext> it = reader.iterator();
            while(it.hasNext()) {
                VariantContext vc = it.next();
                counter++;
            }
            Assert.assertEquals(counter, 2);
            Tribble.indexFile(fakeVCFFile).delete();
            fakeVCFFile.delete();
        }
        catch (IOException e ) {
            throw new RuntimeException(e.getMessage());
        }

    }

    /**
     * create a fake header of known quantity
     * @param metaData           the header lines
     * @param additionalColumns  the additional column names
     * @return a fake VCF header
     */
    public static VCFHeader createFakeHeader(Set<VCFHeaderLine> metaData, Set<String> additionalColumns) {
        metaData.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_0.getFormatString(), VCFHeaderVersion.VCF4_0.getVersionString()));
        metaData.add(new VCFHeaderLine("two", "2"));
        additionalColumns.add("extra1");
        additionalColumns.add("extra2");
        return new VCFHeader(metaData, additionalColumns);
    }

    /**
     * create a fake VCF record
     * @param header the VCF header
     * @return a VCFRecord
     */
    private VariantContext createVC(VCFHeader header) {

       return createVCGeneral(header,"chr1",1);
    }

    private VariantContext createVCGeneral(VCFHeader header,String chrom, int position) {
        List<Allele> alleles = new ArrayList<Allele>();
        Set<String> filters = null;
        Map<String, Object> attributes = new HashMap<String,Object>();
        GenotypesContext genotypes = GenotypesContext.create(header.getGenotypeSamples().size());

        alleles.add(Allele.create("A",true));
        alleles.add(Allele.create("ACC",false));

        attributes.put("DP","50");
        for (String name : header.getGenotypeSamples()) {
            Genotype gt = new GenotypeBuilder(name,alleles.subList(1,2)).GQ(0).attribute("BB", "1").phased(true).make();
            genotypes.add(gt);
        }
        return new VariantContextBuilder("RANDOM", chrom, position, position, alleles)
                .genotypes(genotypes).attributes(attributes).make();
    }


    /**
     * validate a VCF header
     * @param header the header to validate
     */
    public void validateHeader(VCFHeader header) {
        // check the fields
        int index = 0;
        for (VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
            Assert.assertEquals(VCFHeader.HEADER_FIELDS.values()[index], field);
            index++;
        }
        Assert.assertEquals(header.getMetaDataInSortedOrder().size(), metaData.size());
        index = 0;
        for (String key : header.getGenotypeSamples()) {
            Assert.assertTrue(additionalColumns.contains(key));
            index++;
        }
        Assert.assertEquals(index, additionalColumns.size());
    }

    @Test(enabled=true)
    public void TestWritingLargeVCF() throws FileNotFoundException, InterruptedException {

        final Set<String> Columns = new HashSet<String>();
        for (int i = 0; i < 123; i++) {

            Columns.add(String.format("SAMPLE_%d", i));
        }

        final VCFHeader header = createFakeHeader(metaData,Columns);
        final EnumSet<Options> options = EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER,Options.INDEX_ON_THE_FLY);

        final File tempDir = TestUtil.getTempDirecory("VCFWriter", "StaleIndex");

        tempDir.deleteOnExit();

        final File vcf = new File(tempDir, "test.vcf");
        final File vcfIndex = new File(tempDir, "test.vcf.idx");
        final SAMSequenceDictionary dict=createArtificialSequenceDictionary();

        for(int count=1;count<2; count++){
            final VariantContextWriter writer = VariantContextWriterFactory.create(vcf, dict, options);
            writer.writeHeader(header);

            for (int i = 1; i < 17 ; i++) { // write 17 chromosomes
                for (int j = 1; j < 10; j++) { //10 records each
                    writer.add(createVCGeneral(header, String.format("chr%d", i), j * 100));
                }
            }
            writer.close();

            Assert.assertTrue(vcf.lastModified() <= vcfIndex.lastModified());
        }
    }

}

