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

package org.broadinstitute.variant.vcf;

import net.sf.samtools.SAMSequenceDictionary;
import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.Tribble;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.IndexFactory;
import org.broadinstitute.variant.VariantBaseTest;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriterFactory;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.EnumSet;

/**
 * tests out the various functions in the index factory class
 */
public class IndexFactoryUnitTest extends VariantBaseTest {

    File inputFile = new File(variantTestDataRoot + "HiSeq.10000.vcf");
    File outputFile = createTempFile("onTheFlyOutputTest", ".vcf");
    File outputFileIndex = Tribble.indexFile(outputFile);

    private SAMSequenceDictionary dict;

    @BeforeTest
    public void setup() {
        dict = createArtificialSequenceDictionary();
    }

    //
    // test out scoring the indexes
    //
    @Test
    public void testOnTheFlyIndexing1() throws IOException {
        final Index indexFromInputFile = IndexFactory.createDynamicIndex(inputFile, new VCFCodec());
        if ( outputFileIndex.exists() ) {
            System.err.println("Deleting " + outputFileIndex);
            outputFileIndex.delete();
        }

        for ( int maxRecords : Arrays.asList(0, 1, 10, 100, 1000, -1)) {
            final AbstractFeatureReader source = AbstractFeatureReader.getFeatureReader(inputFile.getAbsolutePath(), new VCFCodec(), indexFromInputFile);

            int counter = 0;
            final EnumSet<Options> options = EnumSet.of(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
            VariantContextWriter writer = VariantContextWriterFactory.create(outputFile, dict, options);
            writer.writeHeader((VCFHeader)source.getHeader());
            CloseableTribbleIterator<VariantContext> it = source.iterator();
            while (it.hasNext() && (counter++ < maxRecords || maxRecords == -1) ) {
                VariantContext vc = it.next();
                writer.add(vc);
            }
            writer.close();

            // test that the input index is the same as the one created from the identical input file
            // test that the dynamic index is the same as the output index, which is equal to the input index
            //WalkerTest.assertOnDiskIndexEqualToNewlyCreatedIndex(outputFileIndex, "unittest", outputFile);
        }
    }
}
