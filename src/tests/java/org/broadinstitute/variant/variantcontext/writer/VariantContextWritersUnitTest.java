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


// the imports for unit testing.


import net.sf.samtools.SAMSequenceDictionary;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.bcf2.BCF2Codec;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextTestProvider;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;


public class VariantContextWritersUnitTest extends VariantBaseTest {
    private SAMSequenceDictionary dictionary;

    @BeforeSuite
    public void before() throws IOException {
        dictionary = createArtificialSequenceDictionary();
        VariantContextTestProvider.initializeTests();
    }

    @DataProvider(name = "VariantContextTest_SingleContexts")
    public Object[][] SiteVCsTest() {
        List<Object[]> tests = new ArrayList<Object[]>();
        for ( VariantContextTestProvider.VariantContextTestData testData : VariantContextTestProvider.generateSiteTests() )
            tests.add(new Object[]{testData});
        return tests.toArray(new Object[][]{});
    }

    // --------------------------------------------------------------------------------
    //
    // Test BCF2 reader / writer
    //
    // --------------------------------------------------------------------------------

    @Test(dataProvider = "VariantContextTest_SingleContexts")
    public void testBCF2WriterReader(final VariantContextTestProvider.VariantContextTestData testData) throws IOException {
        VariantContextTestProvider.testReaderWriter(new BCFIOTester(), testData);
    }

    @Test(dataProvider = "VariantContextTest_SingleContexts")
    public void testBCF2WriterReaderMissingGenotypes(final VariantContextTestProvider.VariantContextTestData testData) throws IOException {
        VariantContextTestProvider.testReaderWriterWithMissingGenotypes(new BCFIOTester(), testData);
    }

    private class BCFIOTester extends VariantContextTestProvider.VariantContextIOTest<BCF2Codec> {
        @Override
        public String getExtension() {
            return ".bcf";
        }

        @Override
        public BCF2Codec makeCodec() {
            return new BCF2Codec();
        }

        @Override
        public VariantContextWriter makeWriter(final File file, final EnumSet<Options> baseOptions) {
            return VariantContextWriterFactory.create(file, dictionary, baseOptions);
        }

        @Override
        public VariantContextTestProvider.VariantContextContainer readAllVCs(File input) throws IOException {
            final BCF2Codec codec = this.makeCodec();
            return VariantContextTestProvider.readAllVCs(input, codec);
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Test VCF reader / writer
    //
    // --------------------------------------------------------------------------------

    @Test(enabled = true, dataProvider = "VariantContextTest_SingleContexts")
    public void testVCF4WriterReader(final VariantContextTestProvider.VariantContextTestData testData) throws IOException {
        VariantContextTestProvider.testReaderWriter(new VCFIOTester(), testData);
    }

    @Test(enabled = true, dataProvider = "VariantContextTest_SingleContexts")
    public void testVCF4WriterReaderMissingGenotypes(final VariantContextTestProvider.VariantContextTestData testData) throws IOException {
        VariantContextTestProvider.testReaderWriterWithMissingGenotypes(new VCFIOTester(), testData);
    }

    private class VCFIOTester extends VariantContextTestProvider.VariantContextIOTest<VCFCodec> {
        @Override
        public String getExtension() {
            return ".vcf";
        }

        @Override
        public List<VariantContext> postprocess(final VCFHeader header, final List<VariantContext> vcsAfterIO) {
            final List<VariantContext> fullyDecoded = new ArrayList<VariantContext>(vcsAfterIO.size());

            for ( final VariantContext withStrings : vcsAfterIO )
                fullyDecoded.add(withStrings.fullyDecode(header, false));

            return fullyDecoded;
        }

        @Override
        public VCFCodec makeCodec() {
            return new VCFCodec();
        }

        @Override
        public VariantContextWriter makeWriter(final File file, final EnumSet<Options> baseOptions) {
            return VariantContextWriterFactory.create(file, dictionary, baseOptions);
        }

        @Override
        public VariantContextTestProvider.VariantContextContainer readAllVCs(File input) throws FileNotFoundException {
            final VCFCodec codec = this.makeCodec();
            return VariantContextTestProvider.readAllVCs(input, codec);
        }
    }
}