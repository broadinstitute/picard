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
import org.broad.tribble.util.TabixUtils;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.*;
import org.testng.Assert;
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
    private Set<VCFHeaderLine> metaData;
    private Set<String> additionalColumns;

    /** test, using the writer and reader, that we can output and input a VCF file without problems */
    @Test(dataProvider = "vcfExtensionsDataProvider")
    public void testBasicWriteAndRead(final String extension) throws IOException {
        final File fakeVCFFile = File.createTempFile("testBasicWriteAndRead.", extension);
        fakeVCFFile.deleteOnExit();
        Tribble.indexFile(fakeVCFFile).deleteOnExit();
        metaData = new HashSet<VCFHeaderLine>();
        additionalColumns = new HashSet<String>();
        final SAMSequenceDictionary sequenceDict = createArtificialSequenceDictionary();
        final VCFHeader header = createFakeHeader(metaData, additionalColumns, sequenceDict);
        final VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(fakeVCFFile)
                .setReferenceDictionary(sequenceDict)
                .setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER, Options.INDEX_ON_THE_FLY))
                .build();
        writer.writeHeader(header);
        writer.add(createVC(header));
        writer.add(createVC(header));
        writer.close();
        final VCFCodec codec = new VCFCodec();
        final FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(fakeVCFFile.getAbsolutePath(), codec, false);
        final VCFHeader headerFromFile = (VCFHeader)reader.getHeader();

        int counter = 0;

        // validate what we're reading in
        validateHeader(headerFromFile, sequenceDict);
        
        try {
            final Iterator<VariantContext> it = reader.iterator();
            while(it.hasNext()) {
                it.next();
                counter++;
            }
            Assert.assertEquals(counter, 2);
        }
        catch (final IOException e ) {
            throw new RuntimeException(e.getMessage());
        }

    }

    /**
     * create a fake header of known quantity
     * @param metaData           the header lines
     * @param additionalColumns  the additional column names
     * @return a fake VCF header
     */
    public static VCFHeader createFakeHeader(final Set<VCFHeaderLine> metaData, final Set<String> additionalColumns,
                                             final SAMSequenceDictionary sequenceDict) {
        metaData.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_0.getFormatString(), VCFHeaderVersion.VCF4_0.getVersionString()));
        metaData.add(new VCFHeaderLine("two", "2"));
        additionalColumns.add("extra1");
        additionalColumns.add("extra2");
        final VCFHeader ret = new VCFHeader(metaData, additionalColumns);
        ret.setSequenceDictionary(sequenceDict);
        return ret;
    }

    /**
     * create a fake VCF record
     * @param header the VCF header
     * @return a VCFRecord
     */
    private VariantContext createVC(final VCFHeader header) {

       return createVCGeneral(header,"1",1);
    }

    private VariantContext createVCGeneral(final VCFHeader header, final String chrom, final int position) {
        final List<Allele> alleles = new ArrayList<Allele>();
        final Map<String, Object> attributes = new HashMap<String,Object>();
        final GenotypesContext genotypes = GenotypesContext.create(header.getGenotypeSamples().size());

        alleles.add(Allele.create("A",true));
        alleles.add(Allele.create("ACC",false));

        attributes.put("DP","50");
        for (final String name : header.getGenotypeSamples()) {
            final Genotype gt = new GenotypeBuilder(name,alleles.subList(1,2)).GQ(0).attribute("BB", "1").phased(true).make();
            genotypes.add(gt);
        }
        return new VariantContextBuilder("RANDOM", chrom, position, position, alleles)
                .genotypes(genotypes).attributes(attributes).make();
    }


    /**
     * validate a VCF header
     * @param header the header to validate
     */
    public void validateHeader(final VCFHeader header, final SAMSequenceDictionary sequenceDictionary) {
        // check the fields
        int index = 0;
        for (final VCFHeader.HEADER_FIELDS field : header.getHeaderFields()) {
            Assert.assertEquals(VCFHeader.HEADER_FIELDS.values()[index], field);
            index++;
        }
        Assert.assertEquals(header.getMetaDataInSortedOrder().size(), metaData.size() + sequenceDictionary.size());
        index = 0;
        for (final String key : header.getGenotypeSamples()) {
            Assert.assertTrue(additionalColumns.contains(key));
            index++;
        }
        Assert.assertEquals(index, additionalColumns.size());
    }

    @Test(dataProvider = "vcfExtensionsDataProvider")
    public void TestWritingLargeVCF(final String extension) throws FileNotFoundException, InterruptedException {

        final Set<VCFHeaderLine> metaData = new HashSet<VCFHeaderLine>();
        final Set<String> Columns = new HashSet<String>();
        for (int i = 0; i < 123; i++) {

            Columns.add(String.format("SAMPLE_%d", i));
        }

        final SAMSequenceDictionary dict = createArtificialSequenceDictionary();
        final VCFHeader header = createFakeHeader(metaData,Columns, dict);

        final File tempDir = TestUtil.getTempDirectory("VCFWriter", "StaleIndex");

        tempDir.deleteOnExit();

        final File vcf = new File(tempDir, "test" + extension);
        final String indexExtension;
        if (extension.equals(".vcf.gz")) {
            indexExtension = TabixUtils.STANDARD_INDEX_EXTENSION;
        } else {
            indexExtension = Tribble.STANDARD_INDEX_EXTENSION;
        }
        final File vcfIndex = new File(vcf.getAbsolutePath() + indexExtension);

        for(int count=1;count<2; count++){
            final VariantContextWriter writer =  new VariantContextWriterBuilder()
                    .setOutputFile(vcf)
                    .setReferenceDictionary(dict)
                    .setOptions(EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER, Options.INDEX_ON_THE_FLY))
                    .build();
            writer.writeHeader(header);

            for (int i = 1; i < 17 ; i++) { // write 17 chromosomes
                for (int j = 1; j < 10; j++) { //10 records each
                    writer.add(createVCGeneral(header, String.format("%d", i), j * 100));
                }
            }
            writer.close();

            Assert.assertTrue(vcf.lastModified() <= vcfIndex.lastModified());
        }
    }

    @DataProvider(name = "vcfExtensionsDataProvider")
    public Object[][]vcfExtensionsDataProvider() {
        return new Object[][] {
                // TODO: BCF doesn't work because header is not properly constructed.
                // {".bcf"},
                {".vcf"},
                {".vcf.gz"}
        };
    }
}

