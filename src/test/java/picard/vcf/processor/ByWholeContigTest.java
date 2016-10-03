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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author mccowan
 */
public class ByWholeContigTest {
    public final static File TEST_VCF = new File("testdata/picard/vcf/CEUTrio-indels-bad-samples.vcf");
    
    @Test
    public void test() throws Exception {
        final VcfFileSegmentGenerator.ByWholeContig segmenter = VcfFileSegmentGenerator.ByWholeContig.getInstance();
        int chunkCount = 0;
        for (final VcfFileSegment variantContextCloseableIterator : segmenter.forVcf(TEST_VCF)) {
            chunkCount++;
        }
        Assert.assertEquals(chunkCount, 84);
    }
}
