/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
package picard.vcf.processor;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author mccowan
 */
public class ThreadsafeTest {
    static final int TEN_MILLION = (int) 10e6;
    static final File VCF_WITH_MULTI_ALLELIC_VARIANT_AT_POSITION_10MILLION = new File("testdata/picard/vcf/chunking/multi_allelic_at_10M.vcf");

    @Test
    public void ensureUniqueVariantObservationsEspeciallyMultiAllelicOnesThatAppearAtChunkingBoundaries() {
        final VariantIteratorProducer.Threadsafe iteratorFactory =
                new VariantIteratorProducer.Threadsafe(
                        VcfFileSegmentGenerator.byWholeContigSubdividingWithWidth(TEN_MILLION),
                        Arrays.asList(VCF_WITH_MULTI_ALLELIC_VARIANT_AT_POSITION_10MILLION)
                );
        final Set<String> observed = new HashSet<String>();
        for (final CloseableIterator<VariantContext> i : iteratorFactory.iterators()) {
            while (i.hasNext()) {
                final VariantContext next = i.next();
                Assert.assertTrue(observed.add(next.toString()), "Second observation for " + next.toString());
            }
        }
    }

    /** This test doesn't even test the class, it just makes sure the cornercase test data is really a cornercase */
    @Test
    public void ensureTestDataActuallyHasWideVariantAtTenMillion() {
        final Joiner joiner = Joiner.on(":"); // Cheat: do a string compare
        final VCFFileReader r = new VCFFileReader(VCF_WITH_MULTI_ALLELIC_VARIANT_AT_POSITION_10MILLION);
        Assert.assertEquals(
                joiner.join(r.query("1", TEN_MILLION, TEN_MILLION)),
                joiner.join(r.query("1", TEN_MILLION + 5, TEN_MILLION + 5))
        );
        r.close();
    }
    
    @Test
    public void ensureSameVariantsReadAsSimpleVcfFileIterator() {
        final VariantIteratorProducer.Threadsafe iteratorFactory =
                new VariantIteratorProducer.Threadsafe(
                        VcfFileSegmentGenerator.byWholeContigSubdividingWithWidth(TEN_MILLION),
                        Arrays.asList(VCF_WITH_MULTI_ALLELIC_VARIANT_AT_POSITION_10MILLION)
                );
        final Set<String> observedVcs = new HashSet<String>();
        final Set<String> actual = new HashSet<String>();
        final VCFFileReader actualVcs = new VCFFileReader(VCF_WITH_MULTI_ALLELIC_VARIANT_AT_POSITION_10MILLION);
        for (final VariantContext actualVc : actualVcs) {
            actual.add(actualVc.toString());
        }

        for (final CloseableIterator<VariantContext> i : iteratorFactory.iterators()) {
            while (i.hasNext()) {
                observedVcs.add(i.next().toString());
            }
        }
        
        Assert.assertEquals(actual, observedVcs);
    }

}
