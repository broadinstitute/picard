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

import com.google.common.collect.Iterables;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

/**
 * @author mccowan
 */
public class VcfFileSegmentGeneratorTest {
    final static Log LOG = Log.getInstance(VcfFileSegmentGeneratorTest.class);
    
    final File VCF_WITH_LOGS_OF_GAPS =  new File("testdata/picard/vcf/chunking/multi_allelic_at_10M.vcf");
    final int TEN_MILLION = (int) 10e6;

    @Test
    public void ensureOverlapExclusionTest() {
        final OverlapDetector<Interval> oneTinyIntervalDetector = new OverlapDetector<Interval>(0, 0);
        final Interval theInterval = new Interval("1", 5, 10);
        oneTinyIntervalDetector.addLhs(theInterval, theInterval);
        final VcfFileSegmentGenerator noFilter = VcfFileSegmentGenerator.byWholeContigSubdividingWithWidth(TEN_MILLION);
        Assert.assertEquals(Iterables.size(noFilter.forVcf(VCF_WITH_LOGS_OF_GAPS)), 382); // The number of subdivisions of 10 million of this vcf
        
        final VcfFileSegmentGenerator allFiltered = VcfFileSegmentGenerator.excludingNonOverlaps(noFilter, oneTinyIntervalDetector);
        Assert.assertEquals(Iterables.size(allFiltered.forVcf(VCF_WITH_LOGS_OF_GAPS)), 1);
    }
}
