/*
 * The MIT License
 *
 * Copyright (c) 2014 The Broad Institute
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
package org.broadinstitute.variant.variantcontext.writer;

import org.broad.tribble.AbstractFeatureReader;
import org.broad.tribble.CloseableTribbleIterator;
import org.broad.tribble.FeatureReader;
import org.broad.tribble.index.tabix.TabixIndex;
import org.broad.tribble.util.TabixUtils;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCF3Codec;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.testng.annotations.Test;

import java.io.File;
import java.util.EnumSet;

public class TabixOnTheFlyIndexCreationTest {
    private static final File SMALL_VCF = new File("testdata/tribble/tabix/trioDup.vcf.gz");
    @Test
    public void simpleTest() throws Exception {
        final VCF3Codec codec = new VCF3Codec();
        final FeatureReader<VariantContext> reader = AbstractFeatureReader.getFeatureReader(SMALL_VCF.getAbsolutePath(), codec, false);
        final VCFHeader headerFromFile = (VCFHeader)reader.getHeader();
        final File vcf = File.createTempFile("TabixOnTheFlyIndexCreationTest.", ".vcf.gz");
        final File tabix = new File(vcf.getAbsolutePath() + TabixUtils.STANDARD_INDEX_EXTENSION);
        vcf.deleteOnExit();
        tabix.deleteOnExit();
        final VariantContextWriter vcfWriter = new VariantContextWriterBuilder()
                .setOutputFile(vcf)
                .setReferenceDictionary(headerFromFile.getSequenceDictionary())
                .setOptions(EnumSet.of(Options.INDEX_ON_THE_FLY, Options.ALLOW_MISSING_FIELDS_IN_HEADER))
                .build();
        vcfWriter.writeHeader(headerFromFile);
        final CloseableTribbleIterator<VariantContext> it = reader.iterator();
        while (it.hasNext()) {
            vcfWriter.add(it.next());
        }
        it.close();
        vcfWriter.close();

        // Hard to validate, so just confirm that index can be read.
        new TabixIndex(tabix);
    }
}
