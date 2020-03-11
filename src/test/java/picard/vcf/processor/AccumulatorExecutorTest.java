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

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author mccowan
 */
public class AccumulatorExecutorTest {
    final static List<File> TEST_VCFS = Arrays.asList(
            new File("testdata/picard/vcf/CEUTrio-indels-bad-samples.vcf"),
            new File("testdata/picard/vcf/CEUTrio-indels-dissimilar-contigs.vcf"),
            new File("testdata/picard/vcf/CEUTrio-merged-indels-snps.vcf")
    );

    @Test
    public void test() throws Exception {
        // Fist, read variants via a known, functional mechanism
        final Set<String> actualVariantContextStrings = Collections.synchronizedSet(new HashSet<String>());
        for (final File testVcf : TEST_VCFS) {
            for (final VariantContext variantContext : new VCFFileReader(testVcf)) {
                actualVariantContextStrings.add(variantContext.toString());
            }
        }

        // Then ensure for a variety of thread counts that we observe the same variants
        for (int i = 1; i <= 24; i++) {
            final Set<String> observedVariantContextStrings = Collections.synchronizedSet(new HashSet<String>());
            final VariantAccumulatorExecutor executor = new VariantAccumulatorExecutor.MultiThreadedChunkBased(
                    i,
                    VariantIteratorProducer.byHundredMegabaseChunks(TEST_VCFS),
                    new VariantProcessor.AccumulatorGenerator() {
                        @Override
                        public VariantProcessor.Accumulator build() {
                            return new VariantProcessor.Accumulator() {
                                @Override
                                public void accumulate(final VariantContext vc) {
                                    observedVariantContextStrings.add(vc.toString());
                                }

                                @Override
                                public Object result() {
                                    return null;
                                }
                            };
                        }
                    }
            );
            executor.start();
            executor.awaitCompletion();
            Assert.assertTrue(actualVariantContextStrings.equals(observedVariantContextStrings));
        }
    }
}
