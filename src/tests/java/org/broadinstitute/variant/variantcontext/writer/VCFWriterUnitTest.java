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

import net.sf.picard.reference.IndexedFastaSequenceFile;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.Tribble;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderVersion;
import org.broadinstitute.variant.variantcontext.*;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;


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
        return new VariantContextBuilder("RANDOM", "chr1", 1, 1, alleles)
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

    @DataProvider(name = "VCFWriterDoubleFormatTestData")
    public Object[][] makeVCFWriterDoubleFormatTestData() {
        List<Object[]> tests = new ArrayList<Object[]>();
        tests.add(new Object[]{1.0, "1.00"});
        tests.add(new Object[]{10.1, "10.10"});
        tests.add(new Object[]{10.01, "10.01"});
        tests.add(new Object[]{10.012, "10.01"});
        tests.add(new Object[]{10.015, "10.02"});
        tests.add(new Object[]{0.0, "0.00"});
        tests.add(new Object[]{0.5, "0.500"});
        tests.add(new Object[]{0.55, "0.550"});
        tests.add(new Object[]{0.555, "0.555"});
        tests.add(new Object[]{0.5555, "0.556"});
        tests.add(new Object[]{0.1, "0.100"});
        tests.add(new Object[]{0.050, "0.050"});
        tests.add(new Object[]{0.010, "0.010"});
        tests.add(new Object[]{0.012, "0.012"});
        tests.add(new Object[]{0.0012, "1.200e-03"});
        tests.add(new Object[]{1.2e-4, "1.200e-04"});
        tests.add(new Object[]{1.21e-4, "1.210e-04"});
        tests.add(new Object[]{1.212e-5, "1.212e-05"});
        tests.add(new Object[]{1.2123e-6, "1.212e-06"});
        tests.add(new Object[]{Double.POSITIVE_INFINITY, "Infinity"});
        tests.add(new Object[]{Double.NEGATIVE_INFINITY, "-Infinity"});
        tests.add(new Object[]{Double.NaN, "NaN"});
        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "VCFWriterDoubleFormatTestData")
    public void testVCFWriterDoubleFormatTestData(final double d, final String expected) {
        Assert.assertEquals(VCFWriter.formatVCFDouble(d), expected, "Failed to pretty print double in VCFWriter");
    }
}

